#include<stdio.h>
#include<stdlib.h>

float A = 20;
float B = 1;
float C = 30;
float D = 1;

void make_step(float *x, float *y, float *x_old, float *y_old, float *t, float dt);
float y_prime(float x,float y);
float x_prime(float x,float y);

int main(int argc, char **argv){
  float x = atof(argv[1]);
  float y = atof(argv[2]);
  int points =10000;
  float t=0;
  float dt=1.0/points;
  float x_old=x;
  float y_old=y;

  int i;
  char filename[100];
  sprintf(filename, "Data/poblaciones_%d_%d.dat", (int)x,(int)y);
  FILE *out = fopen(filename,"w");
  fprintf(out, "%f %f %f\n", t,x,y);

  for (i=0; i<points;i++){
    make_step(&x, &y,&x_old,&y_old,&t,dt);
    fprintf(out, "%f %f %f \n", t,x,y);
  }

  return 0;
}

float x_prime(float x,float y){
  return A*x - B*x*y;
}

float y_prime(float x,float y){
  return -C*y + D*x*y;
}

void make_step(float *x, float *y, float *x_old, float *y_old, float *t, float dt){
  float kx1 = x_prime(*x_old,*y_old);
  float ky1 = y_prime(*x_old,*y_old);
  
  float x1 = *x_old+dt/2.0*kx1;
  float y1 = *y_old+dt/2.0*ky1;
  float kx2 = x_prime(x1,y1);
  float ky2 = y_prime(x1,y1);

  float x2 = *x_old+dt/2.0*kx2;
  float y2 = *y_old+dt/2.0*ky2;
  float kx3 = x_prime(x2,y2);
  float ky3 = y_prime(x2,y2);

  float x3 = *x_old+dt*kx3;
  float y3 = *y_old+dt*ky3;
  float kx4 = x_prime(x3,y3);
  float ky4 = y_prime(x3,y3);

  float kx=1.0/6.0*(kx1+2.0*kx2+2.0*kx3+kx4);
  float ky=1.0/6.0*(ky1+2.0*ky2+2.0*ky3+ky4);

  *t=*t+dt;
  *y= *y_old + ky*dt;
  *x= *x_old + kx*dt;
  *x_old=*x;
  *y_old=*y;
}
