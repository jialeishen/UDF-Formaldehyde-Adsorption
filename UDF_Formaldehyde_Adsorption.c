/**************************************************************************
		UDF of formaldehyde adsorption by room surfaces
@author:Jialei Shen
@contact: www.jialeishen.com
@latest:2022.05.07
This is an UDF file to simulate indoor formaldehyde distribution under the 
condition with several sorptive surfaces indoors (walls and ceiling). A 
sink term of formaldehyde is added to the species transport governing 
equation in Fluent.

Descriptions of the simulation case:
 - Room dimensions: 10' x 13' x 8' high
 - Supply air diffuser: 12" x 2"
 - Return air opening: 12" x 6"
 - Supply and return airflow rate: 54 cfm
 - A constant course of formaldehyde is injected to the supply air through
 0.1 ACH (qi = 5.4 cfm) outdoor flow at 100 ug/m3. And to keep airflow balance 
 in the room, the same amount of airflow is exhausted from the return air
 (at the concentration of formaldehyde in the return/room air)
 - The 4 walls and ceiling surface are sorptive materials with the following 
 characteristics:
 	s = ka * Ca - kd * M	(Eq. 1)
 where, 
 	s - net sorptive rate per unit surface area, ug/(m2s)
	ka - adsorption rate constant, m/s (assuming ka = 0.25 m/s)
	Ca - local concentration right above the surface (near surface)
	kd - desorption rate constant, /s (assuming 0.01 /s)
	M - adsorbed mass per unit surface area, ug/m2 (M = 0 at t = 0)
 - At t = 0, the initial chamber concentration is 0
 - Isothermal condition assumed, T = 23 degC, RH = 50%
 - Airflow and concentration distribution are computed thorugh the CFD model 
 with the sink surfaces as boundary conditions
 - Assuming zero flux for heat and moist
 - Compared to a single-zone (well-mixed) model for pollutants:
 	V * dC/dt = (Ci - CR) * qi - s * A	(Eq. 2)
 where,
 	V - chamber volume, m3
	Ci - supply air concentration (injected air), Ci = 100 ug/m3
	CR - return air concentration, CR = C for well-mixed indoor air
	qi - incoming airflow rate (outdoor airflow rate), qi = 5.4 cfm
	s - sorption rate, ug/(m2s)
	A - total sorptive surface area, m2
**************************************************************************/

#include "udf.h"
#include "sg.h"

#define rho 1.225          //density of air (kg/m3)
#define ka 0.25            //adsorption rate constant (m/s)
#define kd 0.01            //desorption rate constant (1/s)
#define dt 1.              //time step size (s)

real c_outlet;

/************formaldehyde sink term (adsorption & desorption)*************/
/************for the near surface cells***********************************/

DEFINE_SOURCE(formaldehyde_sink_udf,c,t,dS,eqn) // define a sink term
{
	face_t f;
	Thread *tf;
	int n;
	real NV_VEC(A);
	real xc[ND_ND], xf[ND_ND],y0[ND_ND]; // define variables storing the coordinates of cells and surfaces
	real source;
	real xxx,yyy,zzz,xx,yy,zz; 
	C_CENTROID(xc,c,t); //get the centroid coordinates of all (volume) cells
	c_face_loop (c,t,n) // a loop for all face cells
	{
		f=C_FACE(c,t,n); // get all face cells
		tf=C_FACE_THREAD(c,t,n); // get the thread of all face cells
		F_CENTROID(xf,f,tf); // get the centroid coordinates of face cells, and store them to xf
		// round the coordinates to 4 decimal place (in case the obtained coordinates were close but never equal to the desired coordinates, e.g. x = 1.000000000001 instead of x = 1.0)
		xxx=ROUND(xf[0]*10000.0);
		xx=xxx/10000.0;
		yyy=ROUND(xf[1]*10000.0);
		yy=yyy/10000.0;
		zzz=ROUND(xf[2]*10000.0);
		zz=zzz/10000.0;
		// IF the face cells are boundary faces (instead of interior faces), AND (&&) are not on the floor surface (zz != 0.0), AND (&&) are not faces on the diffuser and return surfaces (yy != ..., xx != ...) 
		// which means only the 4 walls and ceiling surface will be set up as sorptive surfaces (sink terms)
		if (THREAD_TYPE(tf)==THREAD_F_WALL && zz!=0.0 && zz!=0.15 && yy!=0.15 && yy!=-0.15 && xx!=0.05)  
		{
			NV_VV(y0,=,xc,-,xf); // vector operation: y0 = xc - xf, y0 is the distance from the centroid of the first cell on the boundary to the centroid of the boundary surface (aka the half height of the boundary cell)
			F_AREA(A,f,tf); // get the area of the boundary surface area (of walls and ceiling)
			// define a user defined memory (UDM, index = 1 here) to store the net sorption rate per unit surface area (s in Eq. 1 and 2)
			// the variable is determined through Eq. 1: s = ka * Ca - kd * M
			// Ca: C_YI(c,t,0)*rho*NV_MAG(A)/C_VOLUME(c,t) --> local concentration of the pollutant (index = 0) at the boundary cell. C_YI(c,t,0) returns the volumetric concentration of pollutant 0, it needs to be converted into a mass concentration as in the equation
			// M: C_UDMI(c,t,0) --> adsorbed mass per unit surface area by the sorptive surfaces, is stored by another UDM (index = 0), which is determined by another UDF below
			C_UDMI(c,t,1)=ka*C_YI(c,t,0)*rho*NV_MAG(A)/C_VOLUME(c,t)-kd*C_UDMI(c,t,0);
			// set up the source term (negative value means a sink term)
			source=-1*C_UDMI(c,t,1); // required by Fluent
			dS[eqn]=-1*ka*rho*NV_MAG(A)/C_VOLUME(c,t); // required by Fluent
		}
		else // for other surfaces, no sorptive/sink term
		{
			source=0.;
			dS[eqn]=0.;
		}
	}
	return source;
}

/********calculate adsorbed mass by each cell on sorptive surfaces********/
/********executed at the end of each time step****************************/
DEFINE_EXECUTE_AT_END(execute_at_end) // an execute_at_end function needs to be defined
{
	Domain *d;
	Thread *t;
	cell_t c;
	d = Get_Domain(1);
	thread_loop_c(t,d) // loop the entire cells
	{
		if (FLUID_THREAD_P(t))
		{
			begin_c_loop(c,t)
			{
				// store the adsorbed mass to a UDM (index = 0)
				// M = s * dt
				C_UDMI(c,t,0)+=C_UDMI(c,t,1)*dt; 
			}
			end_c_loop(c,t)
		}
	}
}

/***************calculate adsorbed mass of each surface*******************/
/***************executed at the end of each time step*********************/
DEFINE_EXECUTE_AT_END(adsorption) // an execute_at_end function needs to be defined
{
	Domain *d;
	Thread *t;
	Thread *tf;
	
	FILE *fp_ad; 
	
	int n;
	
	cell_t c;
	face_t f;
	
	real NV_VEC(A);
	real xc[ND_ND], xf[ND_ND];
	real xxx,yyy,zzz,xx,yy,zz;
	real S_ceiling=0.; // net sorption rate at each time step, s in Eq. 1 & 2
	real S_wall1=0.;
	real S_wall2=0.;
	real S_wall3=0.;
	real S_wall4=0.;
	real M_ceiling=0.; // accumulated/adsorbed mass by the sorptive surfaces at each time step, M in Eq. 1 & 2
	real M_wall1=0.;
	real M_wall2=0.;
	real M_wall3=0.;
	real M_wall4=0.;
	real A_ceiling=0.; // surface area
	real A_wall1=0.;
	real A_wall2=0.;
	real A_wall3=0.;
	real A_wall4=0.;
	
	d = Get_Domain(1);
	
	fp_ad=fopen("adsorption.txt","a"); // write the desired data to the file
	
	thread_loop_c(t,d)
	{
		if (FLUID_THREAD_P(t))
		{
			begin_c_loop(c,t)
			{
				c_face_loop (c,t,n) // loop all boundary face cells
				{
					C_CENTROID(xc,c,t); 
					f=C_FACE(c,t,n);
					tf=C_FACE_THREAD(c,t,n);
					F_CENTROID(xf,f,tf);
					xxx=ROUND(xf[0]*10000.0); // rounding the coordinates
					xx=xxx/10000.0;
					yyy=ROUND(xf[1]*10000.0);
					yy=yyy/10000.0;
					zzz=ROUND(xf[2]*10000.0);
					zz=zzz/10000.0;
					if (THREAD_TYPE(tf)==THREAD_F_WALL && zz==2.44) // ceiling
					{
						F_AREA(A,f,tf);
						S_ceiling+=C_UDMI(c,t,1)*C_VOLUME(c,t); // calculate the s of the ceiling 
						M_ceiling+=C_UDMI(c,t,0)*C_VOLUME(c,t); // calculate the M of the ceiling
						A_ceiling+=NV_MAG(A); // calculate the ceiling area
					}
					else if (THREAD_TYPE(tf)==THREAD_F_WALL && xx==0.0) // wall 1
					{
						F_AREA(A,f,tf);
						S_wall1+=C_UDMI(c,t,1)*C_VOLUME(c,t);
						M_wall1+=C_UDMI(c,t,0)*C_VOLUME(c,t);
						A_wall1+=NV_MAG(A);
					}
					else if (THREAD_TYPE(tf)==THREAD_F_WALL && yy==-1.525) // wall 2
					{
						F_AREA(A,f,tf);
						S_wall2+=C_UDMI(c,t,1)*C_VOLUME(c,t);
						M_wall2+=C_UDMI(c,t,0)*C_VOLUME(c,t);
						A_wall2+=NV_MAG(A);
					}
					else if (THREAD_TYPE(tf)==THREAD_F_WALL && xx==3.96) // wall 3
					{
						F_AREA(A,f,tf);
						S_wall3+=C_UDMI(c,t,1)*C_VOLUME(c,t);
						M_wall3+=C_UDMI(c,t,0)*C_VOLUME(c,t);
						A_wall3+=NV_MAG(A);
					}
					else if (THREAD_TYPE(tf)==THREAD_F_WALL && yy==1.525) // wall 4
					{
						F_AREA(A,f,tf);
						S_wall4+=C_UDMI(c,t,1)*C_VOLUME(c,t);
						M_wall4+=C_UDMI(c,t,0)*C_VOLUME(c,t);
						A_wall4+=NV_MAG(A);
					}
				}
			}
			end_c_loop(c,t)
		}
	}
	S_ceiling=S_ceiling/A_ceiling; // per unit surface area
	S_wall1=S_wall1/A_wall1;
	S_wall2=S_wall2/A_wall2;
	S_wall3=S_wall3/A_wall3;
	S_wall4=S_wall4/A_wall4;
	M_ceiling=M_ceiling/A_ceiling; // per unit surface area
	M_wall1=M_wall1/A_wall1;
	M_wall2=M_wall2/A_wall2;
	M_wall3=M_wall3/A_wall3;
	M_wall4=M_wall4/A_wall4;
	
	// write to the file
	fprintf(fp_ad,"%g %g %g %g %g %g %g %g %g %g\n",S_ceiling,S_wall1,S_wall2,S_wall3,S_wall4,M_ceiling,M_wall1,M_wall2,M_wall3,M_wall4);
	fclose(fp_ad);
}

/***************calculate return air concentration************************/
/***************executed at the end of each time step*********************/
DEFINE_EXECUTE_AT_END(outflow_concentration) // an execute_at_end function needs to be defined
{
	Domain *d;
	Thread *t;
	face_t f;
	
	FILE *fp_out;
	
	real NV_VEC(A);
	real c_sum=0.;
	real A_outlet=0.;
	real x[ND_ND];
	real xxx,yyy,zzz,xx,yy,zz;

	d = Get_Domain(1);
	
	fp_out=fopen("out_concentration.txt","a"); // write data to the txt file
	
	thread_loop_f(t,d)
	{
		begin_f_loop(f,t)
		{
			F_CENTROID(x,f,t);
			xxx=ROUND(x[0]*1000.0);
			yyy=ROUND(x[1]*1000.0);
			zzz=ROUND(x[2]*1000.0);
			xx=xxx/1000.0;
			yy=yyy/1000.0;
			zz=zzz/1000.0;
			if(xx==0.05 && yy>=-0.15 && yy<=0.15 && zz>=0 && zz<=0.15)
			{
				F_AREA(A,f,t); 
				c_sum+=F_YI(f,t,0)*NV_MAG(A); // average pollutant concentration of the return air
				A_outlet+=NV_MAG(A);
			}
		}
		end_f_loop(f,t)
	}
	c_outlet=rho*c_sum/A_outlet; // c_outlet is a global variable that can be used by all functions, it will be used to define the inlet concentration as the AHU has a certain ratio of recirculated air
	fprintf(fp_out,"%g\n",c_outlet);
	fclose(fp_out);
}

/***************calculate supply air concentration************************/
/***************executed at the end of each time step*********************/
DEFINE_EXECUTE_AT_END(inflow_concentration) // an execute_at_end function needs to be defined
{
	Domain *d;
	Thread *t;
	face_t f;
	
	FILE *fp_in;

	real NV_VEC(A);
	real c_sum=0.;
	real A_inlet=0.;
	real c_in;
	real x[ND_ND];
	real xxx,yyy,zzz,xx,yy,zz;
	
	d = Get_Domain(1);
	
	fp_in=fopen("in_concentration.txt","a"); // write data to the txt file
	
	thread_loop_f(t,d)
	{
		begin_f_loop(f,t)
		{
			F_CENTROID(x,f,t);
			xxx=ROUND(x[0]*1000.0);
			yyy=ROUND(x[1]*1000.0);
			zzz=ROUND(x[2]*1000.0);
			xx=xxx/1000.0;
			yy=yyy/1000.0;
			zz=zzz/1000.0;
			if(xx>=0.0 && xx<=0.05 && yy>=-0.15 && yy<=0.15 && zz==0.15)
			{
				F_AREA(A,f,t);
				c_sum+=F_YI(f,t,0)*NV_MAG(A); // average pollutant concentration of the supply air
				A_inlet+=NV_MAG(A);
			}
		}
		end_f_loop(f,t)
	}
	c_in=rho*c_sum/A_inlet;
	fprintf(fp_in,"%g\n",c_in);
	fclose(fp_in);
}

/***************define supply air concentration***************************/
DEFINE_PROFILE(c_inlet,t,i) // a profile function for defining the inlet concentration
{
	face_t f;
	begin_f_loop(f,t)
	{
		// inlet pollutant concentration = return air concentration * recirculated ratio + outdoor air concentration * outdoor ratio
		F_PROFILE(f,t,i)=((500.0e-9*8.667+(86.67-8.667)*c_outlet)/86.67)/rho;
	}
	end_f_loop(f,t)
}
