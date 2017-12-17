/**************************************************************************
			UDF of formaldehyde adsorption on room surfaces
@author:Jialei Shen
@e-mail:jialeishen@smail.nju.edu.cn
@latest:2017.12.17
This is an UDF file to simulate indoor formaldehyde adsorption on room
surfaces. A sink term of formaldehyde is added to the species transport
governing equation in Fluent.
**************************************************************************/

#include "udf.h"
#include "sg.h"

#define rho 1.225          //density of air (kg/m3)
#define ka 0.25            //adsorption rate constant (m/s)
#define kd 0.01            //desorption rate constant (1/s)
#define dt 1.              //time step size (s)

real c_outlet;

/************formaldehyde sink term (adsorption & desorption)*************/

DEFINE_SOURCE(formaldehyde_sink_udf,c,t,dS,eqn)
{
	face_t f;
	Thread *tf;
	int n;
	real NV_VEC(A);
	real xc[ND_ND], xf[ND_ND],y0[ND_ND];
	real source;
	real xxx,yyy,zzz,xx,yy,zz;
	C_CENTROID(xc,c,t);
	c_face_loop (c,t,n)
	{
		f=C_FACE(c,t,n);
		tf=C_FACE_THREAD(c,t,n);
		F_CENTROID(xf,f,tf);
		xxx=ROUND(xf[0]*10000.0);
		xx=xxx/10000.0;
		yyy=ROUND(xf[1]*10000.0);
		yy=yyy/10000.0;
		zzz=ROUND(xf[2]*10000.0);
		zz=zzz/10000.0;
		if (THREAD_TYPE(tf)==THREAD_F_WALL && zz!=0.0 && zz!=0.15 && yy!=0.15 && yy!=-0.15 && xx!=0.05)
		{
			NV_VV(y0,=,xc,-,xf);
			F_AREA(A,f,tf);
			C_UDMI(c,t,1)=ka*C_YI(c,t,0)*rho*NV_MAG(A)/C_VOLUME(c,t)-kd*C_UDMI(c,t,0);
			source=-1*C_UDMI(c,t,1);
			dS[eqn]=-1*ka*rho*NV_MAG(A)/C_VOLUME(c,t);
		}
		else
		{
			source=0.;
			dS[eqn]=0.;
		}
	}

	return source;
}

DEFINE_EXECUTE_AT_END(execute_at_end)
{
	Domain *d;
	Thread *t;
	cell_t c;
	d = Get_Domain(1);
	thread_loop_c(t,d)
	{
		if (FLUID_THREAD_P(t))
		{
			begin_c_loop(c,t)
			{
				C_UDMI(c,t,0)+=C_UDMI(c,t,1)*dt;
			}
			end_c_loop(c,t)
		}
	}
}

DEFINE_EXECUTE_AT_END(adsorption)
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
	real S_ceiling=0.;
	real S_wall1=0.;
	real S_wall2=0.;
	real S_wall3=0.;
	real S_wall4=0.;
	real M_ceiling=0.;
	real M_wall1=0.;
	real M_wall2=0.;
	real M_wall3=0.;
	real M_wall4=0.;
	real A_ceiling=0.;
	real A_wall1=0.;
	real A_wall2=0.;
	real A_wall3=0.;
	real A_wall4=0.;
	
	d = Get_Domain(1);
	
	fp_ad=fopen("adsorption.txt","a");
	
	thread_loop_c(t,d)
	{
		if (FLUID_THREAD_P(t))
		{
			begin_c_loop(c,t)
			{
				c_face_loop (c,t,n)
				{
					C_CENTROID(xc,c,t);
					f=C_FACE(c,t,n);
					tf=C_FACE_THREAD(c,t,n);
					F_CENTROID(xf,f,tf);
					xxx=ROUND(xf[0]*10000.0);
					xx=xxx/10000.0;
					yyy=ROUND(xf[1]*10000.0);
					yy=yyy/10000.0;
					zzz=ROUND(xf[2]*10000.0);
					zz=zzz/10000.0;
					if (THREAD_TYPE(tf)==THREAD_F_WALL && zz==2.44)
					{
						F_AREA(A,f,tf);
						S_ceiling+=C_UDMI(c,t,1)*C_VOLUME(c,t);
						M_ceiling+=C_UDMI(c,t,0)*C_VOLUME(c,t);
						A_ceiling+=NV_MAG(A);
					}
					else if (THREAD_TYPE(tf)==THREAD_F_WALL && xx==0.0)
					{
						F_AREA(A,f,tf);
						S_wall1+=C_UDMI(c,t,1)*C_VOLUME(c,t);
						M_wall1+=C_UDMI(c,t,0)*C_VOLUME(c,t);
						A_wall1+=NV_MAG(A);
					}
					else if (THREAD_TYPE(tf)==THREAD_F_WALL && yy==-1.525)
					{
						F_AREA(A,f,tf);
						S_wall2+=C_UDMI(c,t,1)*C_VOLUME(c,t);
						M_wall2+=C_UDMI(c,t,0)*C_VOLUME(c,t);
						A_wall2+=NV_MAG(A);
					}
					else if (THREAD_TYPE(tf)==THREAD_F_WALL && xx==3.96)
					{
						F_AREA(A,f,tf);
						S_wall3+=C_UDMI(c,t,1)*C_VOLUME(c,t);
						M_wall3+=C_UDMI(c,t,0)*C_VOLUME(c,t);
						A_wall3+=NV_MAG(A);
					}
					else if (THREAD_TYPE(tf)==THREAD_F_WALL && yy==1.525)
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
	S_ceiling=S_ceiling/A_ceiling;
	S_wall1=S_wall1/A_wall1;
	S_wall2=S_wall2/A_wall2;
	S_wall3=S_wall3/A_wall3;
	S_wall4=S_wall4/A_wall4;
	M_ceiling=M_ceiling/A_ceiling;
	M_wall1=M_wall1/A_wall1;
	M_wall2=M_wall2/A_wall2;
	M_wall3=M_wall3/A_wall3;
	M_wall4=M_wall4/A_wall4;
	
	fprintf(fp_ad,"%g %g %g %g %g %g %g %g %g %g\n",S_ceiling,S_wall1,S_wall2,S_wall3,S_wall4,M_ceiling,M_wall1,M_wall2,M_wall3,M_wall4);
	fclose(fp_ad);
}

DEFINE_EXECUTE_AT_END(outflow_concentration)
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
	
	fp_out=fopen("out_concentration.txt","a");
	
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
				c_sum+=F_YI(f,t,0)*NV_MAG(A);
				A_outlet+=NV_MAG(A);
			}
		}
		end_f_loop(f,t)
	}
	c_outlet=rho*c_sum/A_outlet;
	fprintf(fp_out,"%g\n",c_outlet);
	fclose(fp_out);
}

DEFINE_EXECUTE_AT_END(inflow_concentration)
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
	
	fp_in=fopen("in_concentration.txt","a");
	
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
				c_sum+=F_YI(f,t,0)*NV_MAG(A);
				A_inlet+=NV_MAG(A);
			}
		}
		end_f_loop(f,t)
	}
	c_in=rho*c_sum/A_inlet;
	fprintf(fp_in,"%g\n",c_in);
	fclose(fp_in);
}

DEFINE_PROFILE(c_inlet,t,i)
{
	face_t f;
	begin_f_loop(f,t)
	{
		F_PROFILE(f,t,i)=((500.0e-9*8.667+(86.67-8.667)*c_outlet)/86.67)/rho;
	}
	end_f_loop(f,t)
}
