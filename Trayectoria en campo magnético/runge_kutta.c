#include<stdio.h>
#include<stdlib.h>
#include<math.h>

double qB0_m = 2874; // qB0/m [1/s]
double Rt = 6378.1E3; // m
double c = 3E8; // m/s^2
double m = 938.272; // MeV/c^2

double vx_prime(double x,double y,double z,double vx,double vy,double vz,double v0);
double vy_prime(double x,double y,double z,double vx,double vy,double vz,double v0);
double vz_prime(double x,double y,double z,double vx,double vy,double vz,double v0);
void make_step(double *x,double *y,double  *z,double *vx,double *vy,double *vz,double *x_old,double  *y_old,double *z_old,double *vx_old,double *vy_old,double  *vz_old,double *t,double  dt, double v0);

int main(int argc, char **argv){
  int points =10000;
  double KE = atof(argv[1]); // KE is MeV
  double alpha = atof(argv[2])*M_PI/180; // alpha is given in degrees
  double x=2*Rt,y=0,z=0;
  double x_old, y_old,z_old,vx_old,vy_old,vz_old;
  double v0 = c*sqrt(1- 1/pow( KE/(2*m) + 1 ,2) ); // Note is constant as the only force is due to the magnetic field. Hence gamma is constant
  double vx = 0, vy = v0*sin(alpha), vz = v0*cos(alpha);
  double t=0;
  double dt = 100.0/points;
  char filename[50];
  sprintf(filename, "trayectoria_%f_%f_runge_kutta.dat",KE,alpha);
  FILE *out= fopen(filename,"w");
  fprintf(out,"%f %f %f %f\n",t,x,y,z);
  t=t+dt;
  x_old=x,y_old=y,z_old=z,vx_old=vx,vy_old=vy,vz_old=vz;
  vx = vx_old + dt*vx_prime(x_old,y_old,z_old,vx_old,vy_old,vz_old,v0);
  vy = vy_old + dt*vy_prime(x_old,y_old,z_old,vx_old,vy_old,vz_old,v0);
  vz = vz_old + dt*vz_prime(x_old,y_old,z_old,vx_old,vy_old,vz_old,v0);
  x=x_old+dt*vx/Rt;
  y=y_old+dt*vy/Rt;
  z=z_old+dt*vz/Rt;
  fprintf(out,"%f %f %f %f\n",t,x,y,z);
  x_old=x,y_old=y,z_old=z,vx_old=vx,vy_old=vy,vz_old=vz;

  int i;

  for(i=1; i<points; i++){
    make_step(&x,&y,&z,&vx,&vy,&vz,&x_old,&y_old,&z_old,&vx_old,&vy_old,&vz_old,&t, dt, v0);
    fprintf(out,"%f %f %f %f\n",t,x,y,z);
  }
  fclose(out);
  return 0;

}

double vx_prime(double x,double y,double z,double vx,double vy,double vz,double v0){
  return pow(Rt,3)*sqrt(1+pow(v0/c,2))*qB0_m*(-vy*(2*pow(z,2)-pow(x,2)-pow(y,2))+3*y*z*vz)/pow(pow(x,2)+pow(y,2)+pow(z,2),5.0/2.0);
}

double vy_prime(double x,double y,double z,double vx,double vy,double vz,double v0){
  return pow(Rt,3)*sqrt(1+pow(v0/c,2))/pow(pow(x,2)+pow(y,2)+pow(z,2),5.0/2.0)*qB0_m*(vx*(2*pow(z,2)-pow(x,2)-pow(y,2))-3*x*z*vz);
}

double vz_prime(double x,double y,double z,double vx,double vy,double vz,double v0){
  return pow(Rt,3)*sqrt(1+pow(v0/c,2))*qB0_m*(-3*vx*y*z+3*x*z*vy)/pow(pow(x,2)+pow(y,2)+pow(z,2),5.0/2.0);
}

void make_step(double *x, double *y,double  *z,double *vx,double *vy,double *vz,double *x_old,double  *y_old, double *z_old,double *vx_old,double *vy_old,double  *vz_old,double *t,double  dt, double v0){

  double ax1 = vx_prime(*x_old,*y_old,*z_old,*vx_old,*vy_old,*vz_old,v0);
  double ay1 = vy_prime(*x_old,*y_old,*z_old,*vx_old,*vy_old,*vz_old,v0);
  double az1 = vz_prime(*x_old,*y_old,*z_old,*vx_old,*vy_old,*vz_old,v0);

  double vx1 = *vx_old+dt/2.0*ax1;
  double vy1 = *vy_old+dt/2.0*ay1;
  double vz1 = *vz_old+dt/2.0*az1;
  double x1 = *x_old+dt/2.0*vx1;
  double y1 = *y_old+dt/2.0*vy1;
  double z1 = *z_old+dt/2.0*vz1;
  double ax2 = vx_prime(x1,y1,z1,vx1,vy1,vz1,v0);
  double ay2 = vy_prime(x1,y1,z1,vx1,vy1,vz1,v0); 
  double az2 = vz_prime(x1,y1,z1,vx1,vy1,vz1,v0);

  double vx2 = *vx_old+dt/2.0*ax2;
  double vy2 = *vy_old+dt/2.0*ay2;
  double vz2 = *vz_old+dt/2.0*az2;
  double x2 = *x_old+dt/2.0*vx2;
  double y2 = *y_old+dt/2.0*vy2;
  double z2 = *z_old+dt/2.0*vz2;
  double ax3 = vx_prime(x2,y2,z2,vx2,vy2,vz2,v0);
  double ay3 = vy_prime(x2,y2,z2,vx2,vy2,vz2,v0);
  double az3 = vz_prime(x2,y2,z2,vx2,vy2,vz2,v0);

  double vx3 = *vx_old+dt*ax3;
  double vy3 = *vy_old+dt*ay3;
  double vz3 = *vz_old+dt*az3;
  double x3 = *x_old+dt*vx3;
  double y3 = *y_old+dt*vy3;
  double z3 = *z_old+dt*vz3;
  double ax4 = vx_prime(x3,y3,z3,vx3,vy3,vz3,v0);
  double ay4 = vy_prime(x3,y3,z3,vx3,vy3,vz3,v0);
  double az4 = vz_prime(x3,y3,z3,vx3,vy3,vz3,v0);

  double ax=1.0/6.0*(ax1+2.0*ax2+2.0*ax3+ax4);
  double ay=1.0/6.0*(ay1+2.0*ay2+2.0*ay3+ay4);
  double az=1.0/6.0*(az1+2.0*az2+2.0*az3+az4);

  double new_vx = 1.0/6.0* (*vx + 2.0*vx1+ 2.0*vx2 +vx3);
  double new_vy = 1.0/6.0* (*vx + 2.0*vx1+ 2.0*vx2 +vx3);
  double new_vz = 1.0/6.0* (*vx + 2.0*vx1+ 2.0*vx2 +vx3);

  *y = *y_old + new_vy*dt;
  *x = *x_old + new_vx*dt;
  *z = *z_old+new_vz*dt;
  *vx = *vx_old+ax*dt;
  *vy = *vy_old+ay*dt;
  *vz = *vz_old+az*dt;

  *x_old = *x;
  *y_old = *y;
  *z_old = *z;
  *vx_old = *vx;
  *vy_old = *vz;
  *vz_old = *vz;
  *t=*t+dt;
}
