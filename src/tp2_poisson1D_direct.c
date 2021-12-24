/******************************************/
/* tp2_poisson1D_direct.c                 */
/* This file contains the main function   */
/* to solve the Poisson 1D problem        */
/******************************************/
#include "lib_poisson1D.h"
int main(int argc,char *argv[])
/* ** argc: Number of arguments */
/* ** argv: Values of arguments */
{
  int ierr;
  int jj;
  int nbpoints, la;
  int ku, kl, kv, lab;
  int *ipiv; //tableau de travaile
  int info; //Argument de diagnositque
  int NRHS;  //nombre de right hand side
  double T0, T1; //Scalaire représantant les températures
  double *RHS, *EX_SOL, *X; //Vecteur
  double *AB; //Matrice creuse en stockage generale band represntant l'operateur de poisson

  double temp, relres;

  NRHS=1;
  nbpoints=102;  //nombre de point total
  la=nbpoints-2; //nombre de point total -2
  T0=-5.0; // 5°C pour x0
  T1=5.0; // 5°C pour xfinal

  printf("--------- Poisson 1D ---------\n\n");
  RHS=(double *) malloc(sizeof(double)*la); // Right Hand Side: terme a droite de l'equation AX=Y
  EX_SOL=(double *) malloc(sizeof(double)*la); //Solution analytique
  X=(double *) malloc(sizeof(double)*la);  // le vecteur X

  set_grid_points_1D(X, &la); // On place l'axe X
  set_dense_RHS_DBC_1D(RHS,&la,&T0,&T1); //On place le vecteur de température initial
  set_analytical_solution_DBC_1D(EX_SOL, X, &la, &T0, &T1); // On calcul le vecteur representant la solution analytique
  write_vec(RHS, &la, "RHS.dat");
  write_vec(EX_SOL, &la, "EX_SOL.dat");
  write_vec(X, &la, "X_grid.dat");

  kv=1; //nombre de diagonal principal
  ku=1; //nombre de sous-diagonal supérieur
  kl=1; //nombre de sous-diagonal inférieur
  lab=kv+kl+ku+1;
	//On initialise le tableau AB qui sert de stockage GB (General Band)
  AB = (double *) malloc(sizeof(double)*lab*la);

  info=0;

  /* working array for pivot used by LU Factorization */
  ipiv = (int *) calloc(la, sizeof(int));

  int row = 0; //

  if (row == 1){ // LAPACK_ROW_MAJOR
    set_GB_operator_rowMajor_poisson1D(AB, &lab, &la);
    write_GB_operator_rowMajor_poisson1D(AB, &lab, &la, "AB_row.dat");

    info = LAPACKE_dgbsv(LAPACK_ROW_MAJOR,la, kl, ku, NRHS, AB, la, ipiv, RHS, NRHS);
    	printf("LAPACK_ROW_MAJOR\n");
  }
  else { // LAPACK_COL_MAJOR
    set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);
    write_GB_operator_colMajor_poisson1D(AB, &lab, &la, "AB_col.dat");

    info = LAPACKE_dgbsv(LAPACK_COL_MAJOR,la, kl, ku, NRHS, AB, lab, ipiv, RHS, la);
        printf("LAPACK_COL_MAJOR\n");
  }


  printf("\n INFO DGBSV = %d\n",info);

  write_xy(RHS, X, &la, "SOL.dat");

  /* Relative residual */
  temp = cblas_ddot(la, RHS, 1, RHS,1);
  temp = sqrt(temp);
  cblas_daxpy(la, -1.0, RHS, 1, EX_SOL, 1);
  relres = cblas_ddot(la, EX_SOL, 1, EX_SOL,1);
  relres = sqrt(relres);
  relres = relres / temp;

  printf("\nThe relative residual error is relres = %e\n",relres);

  free(RHS);
  free(EX_SOL);
  free(X);
  free(AB);
  free(ipiv);

  printf("\n\n--------- End -----------\n");
}
