#include<stdio.h>
#include<stdlib.h>

int main(){
  char *filename = "makefile";
  FILE *out=fopen(filename, "w");
  char actual1[50];
  char actual2[50];
  fprintf(out,"all: ");
  int i;
  for (i=30; i>0; i--){
    sprintf(actual2, "Data/poblaciones_%d_%d.pdf", i, 20);
    fprintf(out, "%s ", actual2);
  }
  fprintf(out,"\n\n");
  for (i=30; i>0; i--){
    sprintf(actual1, "Data/poblaciones_%d_%d.dat", i, 20);
    sprintf(actual2, "Data/poblaciones_%d_%d.pdf", i, 20);
    fprintf(out,"%s: %s %s\n", actual2, actual1, "plot_poblaciones.py");
    fprintf(out,"\tpython %s %s\n","plot_poblaciones.py", actual1);
    fprintf(out,"\n");
    fprintf(out,"%s: %s\n", actual1, "a.out");
    fprintf(out,"\t./a.out %d %d\n",i,20);
    fprintf(out,"\n");
  }
  fprintf(out,"a.out: runge_kutta.c \n");
  fprintf(out,"\tcc runge_kutta.c \n\n");
  fprintf(out, "clear:\n\trm -f Data/*\n\trm -f a.out\n");
  fclose(out);
}
