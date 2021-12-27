/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"

//Ecriture de la matrice de Poisson 1D en stockage general band en priorité ligne
void set_GB_operator_rowMajor_poisson1D(double* AB, int *lab, int *la){

  //TODO
  int ii, jj, kk;

  for(jj=1; jj<(*la);jj++){
    AB[(*la)+jj] = -1.0;
    AB[2*(*la)+jj]= 2.0;
    AB[3*(*la)+jj] = -1.0;
    AB[jj] = 0.0;
  }
    AB[2*(*la)] = 2.0;
    AB[(*la)] = 0.0;
  AB[3*(*la)] = -1.0;
  AB[4*(*la)-1] = 0.0;
}
//Ecriture de la matrice d'identité en priorité ligne
void set_GB_operator_rowMajor_poisson1D_Id(double* AB, int *lab, int *la){
	int ii, jj, kk;

	for(jj=1; jj<(*la);jj++){
		AB[(*la)+jj] = 0.0;
		AB[2*(*la)+jj]= 1.0;
		AB[3*(*la)+jj] = 0.0;
		AB[jj] = 0.0;
	}
		AB[2*(*la)] = 1.0;
		AB[(*la)] = 0.0;
	AB[3*(*la)] = 0.0;
	AB[4*(*la)-1] = 0.0;
  //TOD0
}
//Ecriture de la matrice de Poisson 1D en stockage general band en priorité colonne
void set_GB_operator_colMajor_poisson1D(double* AB, int *lab, int *la, int *kv){
  int ii, jj, kk;
  for (jj=0;jj<(*la);jj++){
    kk = jj*(*lab);
    if (*kv>=0){
      for (ii=0;ii< *kv;ii++){
	AB[kk+ii]=0.0;
      }
    }
    AB[kk+ *kv]=-1.0;
    AB[kk+ *kv+1]=2.0;
    AB[kk+ *kv+2]=-1.0;
  }
  AB[0]=0.0;
  if (*kv == 1) {AB[1]=0;}

  AB[(*lab)*(*la)-1]=0.0;
}
// Ecriture de la matrice identite en stockage general band en priorité colonne
void set_GB_operator_colMajor_poisson1D_Id(double* AB, int *lab, int *la, int *kv){
  int ii, jj, kk;
  for (jj=0;jj<(*la);jj++){
    kk = jj*(*lab);
    if (*kv>=0){
      for (ii=0;ii< *kv;ii++){
	AB[kk+ii]=0.0;
      }
    }
    AB[kk+ *kv]=0.0;
    AB[kk+ *kv+1]=1.0;
    AB[kk+ *kv+2]=0.0;
  }
  AB[1]=0.0;
  AB[(*lab)*(*la)-1]=0.0;
}
//On initialise le systeme avec les conditions aux bord RHS = BC0 , .... ,BC1
void set_dense_RHS_DBC_1D(double* RHS, int* la, double* BC0, double* BC1){
  int jj;
  RHS[0]= *BC0;
  RHS[(*la)-1]= *BC1;
  for (jj=1;jj<(*la)-1;jj++){
    RHS[jj]=0.0;
  }
}
//On met en place la solution analytique
void set_analytical_solution_DBC_1D(double* EX_SOL, double* X, int* la, double* BC0, double* BC1){
  int jj;
  double h, DELTA_T;
  DELTA_T=(*BC1)-(*BC0);
  for (jj=0;jj<(*la);jj++){
    EX_SOL[jj] = (*BC0) + X[jj]*DELTA_T;
  }
}
//On intialise le vecteur X = 0,1 , ... , la
void set_grid_points_1D(double* x, int* la){
  int jj;
  double h;
  h=1.0/(1.0*((*la)+1));
  for (jj=0;jj<(*la);jj++){
    x[jj]=(jj+1)*h;
  }
}
//Ecriture de la matrice de l'operateur de Poisson en 1D en stockage general band priorité ligne
void write_GB_operator_rowMajor_poisson1D(double* AB, int* lab, int* la, char* filename){
  FILE * file;
  int ii,jj;
  file = fopen(filename, "w");
  //Numbering from 1 to la
  if (file != NULL){
    for (ii=0;ii<(*lab);ii++){
      for (jj=0;jj<(*la);jj++){
	fprintf(file,"%lf\t",AB[ii*(*la)+jj]);
      }
      fprintf(file,"\n");
    }
    fclose(file);
  }
  else{
    perror(filename);
  }
}
//Ecriture de la matrice de l'operateur de Poisson en 1D en stockage general band priorité colonne
void write_GB_operator_colMajor_poisson1D(double* AB, int* lab, int* la, char* filename){
  //TODO
}
// Ecriture de vecteur dans un fichier .dat
void write_vec(double* vec, int* la, char* filename){
  int jj;
  FILE * file;
  file = fopen(filename, "w");
  // Numbering from 1 to la
  if (file != NULL){
    for (jj=0;jj<(*la);jj++){
      fprintf(file,"%lf\n",vec[jj]);
    }
    fclose(file);
  }
  else{
    perror(filename);
  }
}

void write_xy(double* vec, double* x, int* la, char* filename){
  int jj;
  FILE * file;
  file = fopen(filename, "w");
  // Numbering from 1 to la
  if (file != NULL){
    for (jj=0;jj<(*la);jj++){
      fprintf(file,"%lf\t%lf\n",x[jj],vec[jj]);
    }
    fclose(file);
  }
  else{
    perror(filename);
  }
}
//Calcul des valeurs propre de la matrice de Poisson
void eig_poisson1D(double* eigval, int *la){
  int ii;
  double scal;
  for (ii=0; ii< *la; ii++){
    scal=(1.0*ii+1.0)*M_PI_2*(1.0/(*la+1));
    eigval[ii]=sin(scal);
    eigval[ii]=4*eigval[ii]*eigval[ii];
  }
}
// Calcule de la valeur propres maximal
double eigmax_poisson1D(int *la){
  double eigmax;
  eigmax=sin(*la *M_PI_2*(1.0/(*la+1)));
  eigmax=4*eigmax*eigmax;
  return eigmax;
}
//Calcul de la valeur propre minimal
double eigmin_poisson1D(int *la){
  double eigmin;
  eigmin=sin(M_PI_2*(1.0/(*la+1)));
  eigmin=4*eigmin*eigmin;
  return eigmin;
}

double richardson_alpha_opt(int *la){
  //TODO
  
}

void richardson_alpha(double *AB, double *RHS, double *X, double *alpha_rich, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit){
  //TODO
}
