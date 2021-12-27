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
  double *Identity;
  double *save;

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
  save=(double *) malloc(sizeof(double)*la);  // le vecteur de sauvegarde

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
  Identity = (double *) malloc(sizeof(double)*lab*la);
  info=0;

  /* working array for pivot used by LU Factorization */
  ipiv = (int *) calloc(la, sizeof(int));

  int row = 0; //
// test de dbgsv faux ku = -1
  if (row == 1){ // LAPACK_ROW_MAJOR
    set_GB_operator_rowMajor_poisson1D(AB, &lab, &la);
    write_GB_operator_rowMajor_poisson1D(AB, &lab, &la, "AB_row.dat");
    // Resout le syteme d'equation A*X= B, ou B  = RHS
    info = LAPACKE_dgbsv(LAPACK_ROW_MAJOR,la, kl, -1, NRHS, AB, la, ipiv, RHS, NRHS); // Si info = 0 , RHS => X
    	printf("LAPACK_ROW_MAJOR\n");
        printf(" INFO DGBSV = %d avec ku =-1\n ",info);
  }
  else { // LAPACK_COL_MAJOR
    set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);
    write_GB_operator_colMajor_poisson1D(AB, &lab, &la, "AB_col.dat");
    // Resout le systeme A*X=B ou B = RHS
    info = LAPACKE_dgbsv(LAPACK_COL_MAJOR,la, kl, -1, NRHS, AB, lab, ipiv, RHS, la);// Si info =0 alors RHS=>X
        printf("LAPACK_COL_MAJOR\n");
          printf(" INFO DGBSV = %d avec ku =-1\n ",info);
  }
  // test faux kl = -1
  if (row == 1){ // LAPACK_ROW_MAJOR
    set_GB_operator_rowMajor_poisson1D(AB, &lab, &la);
    write_GB_operator_rowMajor_poisson1D(AB, &lab, &la, "AB_row.dat");
    // Resout le syteme d'equation A*X= B, ou B  = RHS
    info = LAPACKE_dgbsv(LAPACK_ROW_MAJOR,la, -1, ku, NRHS, AB, la, ipiv, RHS, NRHS); // Si info = 0 , RHS => X
      printf("LAPACK_ROW_MAJOR\n") ;
        printf(" INFO DGBSV = %d avec kl =-1\n ",info);
  }
  else { // LAPACK_COL_MAJOR
    set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);
    write_GB_operator_colMajor_poisson1D(AB, &lab, &la, "AB_col.dat");
    // Resout le systeme A*X=B ou B = RHS
    info = LAPACKE_dgbsv(LAPACK_COL_MAJOR,-1, kl, ku, NRHS, AB, lab, ipiv, RHS, la);// Si info =0 alors RHS=>X
        printf("LAPACK_COL_MAJOR\n");
        printf(" INFO DGBSV = %d avec kl =-1\n ",info);
  }
  // test faux la = -1
  if (row == 1){ // LAPACK_ROW_MAJOR
    set_GB_operator_rowMajor_poisson1D(AB, &lab, &la);
    write_GB_operator_rowMajor_poisson1D(AB, &lab, &la, "AB_row.dat");
    // Resout le syteme d'equation A*X= B, ou B  = RHS
    info = LAPACKE_dgbsv(LAPACK_ROW_MAJOR,-1,kl, ku, NRHS, AB, la, ipiv, RHS, NRHS); // Si info = 0 , RHS => X
      printf("LAPACK_ROW_MAJOR\n") ;
        printf("\n INFO DGBSV = %d avec la  =-1\n ",info);
  }
  else { // LAPACK_COL_MAJOR
    set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);
    write_GB_operator_colMajor_poisson1D(AB, &lab, &la, "AB_col.dat");
    // Resout le systeme A*X=B ou B = RHS
    info = LAPACKE_dgbsv(LAPACK_COL_MAJOR,-1, kl, ku, NRHS, AB, lab, ipiv, RHS, la);// Si info =0 alors RHS=>X
        printf("LAPACK_COL_MAJOR\n");
        printf("\n INFO DGBSV = %d avec la =-1\n ",info);
  }

  // Vrai resultat
  if (row == 1){ // LAPACK_ROW_MAJOR
    set_GB_operator_rowMajor_poisson1D(AB, &lab, &la);
    write_GB_operator_rowMajor_poisson1D(AB, &lab, &la, "AB_row.dat");
    // Resout le syteme d'equation A*X= B, ou B  = RHS
    info = LAPACKE_dgbsv(LAPACK_ROW_MAJOR,la, kl, ku, NRHS, AB, la, ipiv, RHS, NRHS); // Si info = 0 , RHS => X
    	printf("LAPACK_ROW_MAJOR\n");
  }
  else { // LAPACK_COL_MAJOR
    set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);
    write_GB_operator_colMajor_poisson1D(AB, &lab, &la, "AB_col.dat");
    // Resout le systeme A*X=B ou B = RHS
    info = LAPACKE_dgbsv(LAPACK_COL_MAJOR,la, kl, ku, NRHS, AB, lab, ipiv, RHS, la);// Si info =0 alors RHS=>X
        printf("LAPACK_COL_MAJOR\n");
  }


  printf("\n INFO DGBSV = %d\n",info);

  write_xy(RHS, X, &la, "SOL.dat");

  // Relative residual //
  temp = cblas_ddot(la, RHS, 1, RHS,1);
  temp = sqrt(temp);
  cblas_daxpy(la, -1.0, RHS, 1, EX_SOL, 1); // On ecrase EX_SOL par RHS-EX_SOL
  relres = cblas_ddot(la, EX_SOL, 1, EX_SOL,1);
  relres = sqrt(relres);
  relres = relres / temp;

  printf("\nThe relative residual error is relres = %e\n",relres);

  if(row==1){
    set_dense_RHS_DBC_1D(RHS,&la,&T0,&T1); //On place le vecteur de température initial save
    set_analytical_solution_DBC_1D(EX_SOL, X, &la, &T0, &T1); // On calcul le vecteur representant la solution analytique
    info = LAPACKE_dgbsv(LAPACK_ROW_MAJOR,la, kl, ku, NRHS, AB, la, ipiv, RHS, NRHS);
    temp = cblas_ddot(la, save, 1, save,1);
    temp = sqrt(temp);
  //  set_GB_operator_rowMajor_poisson1D(AB, &lab, &la);
  cblas_dgbmv(CblasRowMajor,CblasNoTrans,la,lab,kl,ku,1.0,AB,la,RHS,1,-1.0,save,1); //on fait save = A*RHS-save
  //  dgbmv('N',la,lab,kl,ku,1.0,AB,la,EX_SOL,1,-1.0,save,1);
  relres = cblas_ddot(la, save, 1,save,1);
  relres = sqrt(relres);
  relres = relres / temp;
  printf("\n INFO DGBmV row major\n");
  }else{
    set_dense_RHS_DBC_1D(save,&la,&T0,&T1); //On place le vecteur de température initial
    set_analytical_solution_DBC_1D(EX_SOL, X, &la, &T0, &T1); // On calcul le vecteur representant la solution analytique
    temp = cblas_ddot(la, save, 1, save,1);
    temp = sqrt(temp);
  //  set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);
  cblas_dgbmv(CblasColMajor,CblasNoTrans,la,lab,kl,ku,1.0,AB,la,RHS,1,-1.0,save,1);
  //  dgbmv('N',la,lab,kl,ku,1.0,AB,la,EX_SOL,1,-1.0,save,1);
  relres = cblas_ddot(la, save, 1,save,1);
  relres = sqrt(relres);
  relres = relres / temp;
  printf("\n INFO DGBmV col major \n");
  }

printf("\nThe relative residual error is relres = %e\n",relres);
  free(RHS);
  free(EX_SOL);
  free(X);
  free(AB);
  free(ipiv);
  free(save);
  printf("\n\n--------- End -----------\n");
}
