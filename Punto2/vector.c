#include<stdio.h>
#include<stdlib.h>
#include<math.h>

float qB0_m = 2874; // qB0/m [1/s]
float Rt = 1; // m
float c = 47; // m/s^2
float m = 938.272; // MeV/c^2

void make_step(float *r, float *v, float *r_old, float *v_old, float *t, float v0, float dt);
float *scalar_mult(float lambda, float *a);
float *v_prime(float *r, float *v, float v0);
float *r_prime(float *r, float *v, float v0);
float *qB_mr5vector(float *r);
float norm(float *a);
float *cross_product(float *a, float *b);
float *copy(float *a);
float *add_vector(float *a, float *b);

int main(int argc, char **argv){
  float KE = atof(argv[1]); // KE is MeV
  int points = 1E7;
  float alpha_deg = atof(argv[2]);
  float alpha = alpha_deg*M_PI/180;
  float *r = malloc(3*sizeof(float));
  r[0]=2*Rt, r[1]=0, r[2]=0;
  float *r_old = malloc(3*sizeof(float));
  float v0 = c*sqrt(1- 1/pow( KE/(2*m) + 1 ,2) ); // Note is constant as the only force is due to the magnetic field. Hence gamma is constant
  printf("v0 =%f \n", v0);
  float *v =  malloc(3*sizeof(float));
  v[0] = 0, v[1] = v0*sin(alpha), v[2] = v0*cos(alpha);
  printf("v = ( %f, %f, %f )\n",v[0],v[1],v[2]);
  float *v_old = malloc(3*sizeof(float));
  float t=0;
  float dt = 100.0/points;
  char filename[50];
  sprintf(filename, "trayectoria_%f_%f_runge_kutta.dat",KE,alpha_deg);
  FILE *out= fopen(filename,"w");
  fprintf(out,"%f %f %f %f\n",t,r[0],r[1],r[2]);
  /*
printf("r = ( %f, %f, %f )\n",r[0],r[1],r[2]);
printf("v = ( %f, %f, %f )\n",v[0],v[1],v[2]);*/
  t=t+dt;
  r_old = copy(r);
  r = add_vector(scalar_mult(dt,v),r_old);
  v_old = copy(v);
  float *a = v_prime(r_old,v_old,v0);
  v = add_vector(scalar_mult(dt,a), v_old);
  fprintf(out,"%f %f %f %f\n",t,r[0],r[1],r[2]);
  
  /*printf("r = ( %f, %f, %f )\n",r[0],r[1],r[2]);
    printf("a = ( %f, %f, %f )\n",a[0],a[1],a[2]);*/

  int i;

  FILE *out2 =fopen("speeds","w");

  for(i=1; i<points; i++){

    float *a1 = v_prime(r_old,v_old,v0);

    float *v1 = add_vector(scalar_mult(dt/2.0,a1),v_old);
    float *r1 = add_vector(scalar_mult(dt/2.0,v),r_old);
    float *a2 = v_prime(r1,v1,v0);

    float *v2 = add_vector(scalar_mult(dt/2.0,a2),v_old);
    float *r2 = add_vector(scalar_mult(dt/2.0,v1),r_old);
    float *a3 = v_prime(r2,v2,v0);

    float *v3 = add_vector(scalar_mult(dt,a3),v_old);
    float *r3 = add_vector(scalar_mult(dt,v2),r_old);
    float *a4 = v_prime(r3,v3,v0);

    float *a = scalar_mult(1/6.0, add_vector(scalar_mult(2.0,add_vector(a2,a3)),add_vector(a1, a4)));
    float *newv = scalar_mult(1/6.0, add_vector(scalar_mult(2.0,add_vector(v1,v2)),add_vector(v, v3)));
    /*
      printf("a = (%f, %f, %f)\n", a[0], a[1], a[2]);
      printf("v = (%f, %f, %f)\n", newv[0], newv[1], newv[2]);*/

    t=t+dt;
    r_old = copy(r);
    v_old = copy(v);
    r = add_vector(scalar_mult(dt, newv), r_old);
    v = add_vector(scalar_mult(dt, a), v_old);
    if (i%1000==0)
      fprintf(out2, "%f \n", norm(v));
    fprintf(out,"%f %f %f %f\n",t,r[0],r[1],r[2]);
  }
  fclose(out);
  return 0;

}

float *v_prime(float *r, float *v, float v0){
  float *a= scalar_mult(sqrt(1-pow(v0/c,2)), cross_product(v,qB_mr5vector(r)));
  return a;
}

//copy a vector
float *copy(float *a){
  float *c = malloc(3*sizeof(float));
  int i; 
  for (i=0; i<3; i++)
    c[i]=a[i];
  return c;
}

// q*B/(mr^5) in mks
float *qB_mr5vector(float *r){
  float *B=malloc(3*sizeof(float));
  B[0]=-qB0_m/pow(norm(r),5)*3*r[0]*r[2];
  B[1]=-qB0_m/pow(norm(r),5)*3*r[1]*r[2];
  B[2]=-qB0_m/pow(norm(r),5)*(2*pow(r[2],2)-pow(r[1],2)-pow(r[0],2));
  return B;
}

// Vector norm
float norm(float *a){
  float norm_squared = 0;
  int i;
  for(i=0; i<3; i++)
    norm_squared+=pow(a[i],2);
  return sqrt(norm_squared);
}

// Vector cross product
float *cross_product(float *a, float *b){
  float *c = malloc(3*sizeof(float));
  c[0]=a[1]*b[2]-a[2]*b[1];
  c[1]=a[2]*b[0]-a[0]*b[2];
  c[2]=a[0]*b[1]-a[1]*b[0];
  return c;
}

// Vector addition
float *add_vector(float *a, float *b){
  float *c = malloc(3*sizeof(float));
  int i;
  for(i=0; i<3; i++)
    c[i]=a[i]+b[i];
  return c;
}

//Scalar times vector multiplication
float *scalar_mult(float lambda, float *a){
  float *new = malloc(3*sizeof(float));
  int i;
  for(i=0; i<3; i++)
    new[i]=lambda*a[i];
  return new;
}


void make_step(float *r, float *v, float *r_old, float *v_old, float *t, float v0, float dt){

  float *a1 = v_prime(r_old,v_old,v0);

  float *v1 = add_vector(scalar_mult(dt/2.0,a1),v_old);
  float *r1 = add_vector(scalar_mult(dt/2.0,v),r_old);
  float *a2 = v_prime(r1,v1,v0);

  float *v2 = add_vector(scalar_mult(dt/2.0,a2),v_old);
  float *r2 = add_vector(scalar_mult(dt/2.0,v1),r_old);
  float *a3 = v_prime(r2,v2,v0);

  float *v3 = add_vector(scalar_mult(dt,a3),v_old);
  float *r3 = add_vector(scalar_mult(dt,v2),r_old);
  float *a4 = v_prime(r3,v3,v0);

  float *a = scalar_mult(1/6.0, add_vector(scalar_mult(2.0,add_vector(a2,a3)),add_vector(a1, a4)));
  float *newv = scalar_mult(1/6.0, add_vector(scalar_mult(2.0,add_vector(v1,v2)),add_vector(v, v3)));

  printf("a = ( %f, %f, %f )\n",a[0],a[1],a[2]);

  *t=*t+dt;
  r_old = copy(r);
  v_old = copy(v);
  r = add_vector(scalar_mult(dt, newv), r_old);
  v = add_vector(scalar_mult(dt, a), v_old);

  printf("r = ( %f, %f, %f ) \n", r[0], r[1], r[2]);
}
