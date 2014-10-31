#include<stdio.h>
#include<stdlib.h>
#include<math.h>

float qB0_m = 2874; // qB0/m [1/s]
float c = 47;

float *v_prime(float *r, float *v, float v0);
float *copy(float *a);
float norm(float *a);
float *scalar_mult(float lambda, float *a);
float *cross_product(float *a, float *b);
float *add_vector(float *a, float *b);
float *qB_mr5vector(float *r, float *v);

int main(){
  float *a = malloc(3*sizeof(float));
  a[0]=2, a[1]=2.41, a[2]=4.18;
  float lambda = 2;
  float *b =  scalar_mult(lambda,a);
  b[0]=0,b[1]=1, b[2]=1;
  float *cross = cross_product(a,b);
  float *sum = add_vector(a,b);
  float *B = qB_mr5vector(a, b);
  printf("a =( %f, %f, %f ) \n", a[0], a[1], a[2]);
  printf("b =( %f, %f, %f ) \n", b[0], b[1], b[2]);
  printf("norm_a = %f\n", norm(a));
  printf("B= ( %f, %f, %f ) \n", B[0], B[1], B[2]);
  float *copy_B = copy(B);
  printf("copyB= ( %f, %f, %f ) \n",copy_B[0], copy_B[1], copy_B[2]);
  v_prime(a,b,norm(b));
  return 0;
}

float *qB_mr5vector(float *r, float *v){
  float *B=malloc(3*sizeof(float));
  printf("qbo_m = %f\n",qB0_m);
  B[0]=-qB0_m/pow(norm(r),5)*3*r[0]*r[2];
  B[1]=-qB0_m*3*r[1]*r[2]/pow(norm(r),5);
  B[2]=-qB0_m/pow(norm(r),5)*(2*pow(r[2],2)-pow(r[1],2)-pow(r[0],2));
  return B;
}


float *scalar_mult(float lambda, float *a){
  float *new = malloc(3*sizeof(float));
  int i;
  for(i=0; i<3; i++)
    new[i]=lambda*a[i];
  return new;
}

float *cross_product(float *a, float *b){
  float *c = malloc(3*sizeof(float));
  c[0]=a[1]*b[2]-a[2]*b[1];
  c[1]=a[2]*b[0]-a[0]*b[2];
  c[2]=a[0]*b[1]-a[1]*b[0];
  return c;
}

float *add_vector(float *a, float *b){
  float *c = malloc(3*sizeof(float));
  int i;
  for(i=0; i<3; i++)
    c[i]=a[i]+b[i];
  return c;
}

float norm(float *a){
  float norm_squared = 0;
  int i;
  for(i=0; i<3; i++)
    norm_squared+=pow(a[i],2);
  return sqrt(norm_squared);
}

float *v_prime(float *r, float *v, float v0){
  float *a= scalar_mult(sqrt(1-pow(v0/c,2)),cross_product(qB_mr5vector(r,v),v));
  printf("a in v_prime= ( %f, %f, %f ) \n", a[0], a[1], a[2]);
  return a;
}

float *copy(float *a){
  float *c = malloc(3*sizeof(float));
  int i; 
  for (i=0; i<3; i++)
    c[i]=a[i];
  return c;
}
