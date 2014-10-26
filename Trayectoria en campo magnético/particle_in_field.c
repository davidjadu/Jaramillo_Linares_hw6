#include<stdio.h>
#include<stdlib.h>
#include<math.h>

double qB0_m = 2874; // qB0/m [1/s]
double Rt = 6378.1E3; // m
double c = 3E8; // m/s^2
double m = 0.938272; // GeV/c^2

double vx_prime(double x,double y,double z,double vx,double vy,double vz,double v0);
double vy_prime(double x,double y,double z,double vx,double vy,double vz,double v0);
double vz_prime(double x,double y,double z,double vx,double vy,double vz,double v0);

int main(int argc, char **argv){
  int points =10000;
  double KE = atof(argv[1]); // KE is GeV
  double alpha = atof(argv[2]);
  double x=2*Rt,y=0,z=0;
  double x_old, y_old,z_old,vx_old,vy_old,vz_old,x_old2, y_old2,z_old2,vx_old2,vy_old2,vz_old2;
  double v0 = c*sqrt(1- 1/pow( KE/(2*m) + 1 ,2) ); // Note is constant as the only force is due to the magnetic field. Hence gamma is constant
  double vx = 0, vy = v0*sin(alpha), vz = v0*cos(alpha);
  double t=0;
  double dt = 100.0/points;
  printf("%f\n",v0);
  char filename[50];
  sprintf(filename, "trayectoria_%f_%f.dat",KE,alpha);
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
  x_old2=x_old,y_old2=y_old,z_old2=z_old,vx_old2=vx_old,vy_old2=vy_old,vz_old2=vz_old;
  x_old=x,y_old=y,z_old=z,vx_old=vx,vy_old=vy,vz_old=vz;

  int i;

  //Leap frog
  for(i=1; i<points; i++){
    vx=vx_old2 + 2*dt*vx_prime(x_old,y_old,z_old,vx_old,vy_old,vz_old,v0);
    vy=vy_old2 + 2*dt*vy_prime(x_old,y_old,z_old,vx_old,vy_old,vz_old,v0);
    vz=vz_old2 + 2*dt*vz_prime(x_old,y_old,z_old,vx_old,vy_old,vz_old,v0);
    x = x_old2 + 2*dt*vx;
    y = y_old2 + 2*dt*vy;
    z = z_old2 + 2*dt*vz;
    t=t+dt;
    fprintf(out,"%f %f %f %f\n",t,x,y,z);
    x_old2=x_old,y_old2=y_old,z_old2=z_old,vx_old2=vx_old,vy_old2=vy_old,vz_old2=vz_old;
    x_old=x,y_old=y,z_old=z , vx_old=vx ,vy_old=vy,vz_old=vz;
  }
  fclose(out);
  return 0;

}

double vx_prime(double x,double y,double z,double vx,double vy,double vz,double v0){
  return  sqrt(1+pow(v0/c,2))*qB0_m*(-vy*(2*pow(z,2)-pow(x,2)-pow(y,2))+3*y*z*vz)/pow(pow(x,2)+pow(y,2)+pow(z,2),5.0/2.0);
}

double vy_prime(double x,double y,double z,double vx,double vy,double vz,double v0){
  return  sqrt(1+pow(v0/c,2))*qB0_m*(vx*(2*pow(z,2)-pow(x,2)-pow(y,2))-3*x*z*vz)/pow(pow(x,2)+pow(y,2)+pow(z,2),5.0/2.0);
}

double vz_prime(double x,double y,double z,double vx,double vy,double vz,double v0){
  return sqrt(1+pow(v0/c,2))*qB0_m*(-3*vx*y*z+3*x*z*vy)/pow(pow(x,2)+pow(y,2)+pow(z,2),5.0/2.0);
}
