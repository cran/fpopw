

#include "colibri.h"
#include<R_ext/Arith.h>

// this function is visible by R
extern "C" {
void colibri_op_R_c (double *profil, int *nbi, double *lambda_, double *mini, double *maxi, int *origine,
double *cout_n){
    colibri_op_c (profil, nbi, lambda_, mini, maxi, origine, cout_n);
  }

void colibri_op_weights_R_c (double *profil, double *weights, int *nbi, double *lambda_, double *mini, double *maxi, int *origine,
double *cout_n){
    colibri_op_weight_c (profil, weights, nbi, lambda_, mini, maxi, origine, cout_n);
  }

void colibri_sn_R_c (double *profil, int *nbi, int *Kmax_, double *mini, double *maxi, int *origine,
double *cout_n, double *allCost){
	colibri_sn_c (profil, nbi, Kmax_, mini, maxi, origine, cout_n, allCost);
  }

void colibri_sn_weights_R_c (double *profil, double *weights, int *nbi, int *Kmax_, double *mini, double *maxi, int *origine,
double *cout_n, double *allCost){
	colibri_sn_weight_c (profil, weights, nbi, Kmax_, mini, maxi, origine, cout_n, allCost);
  }

void colibri_sn_weights_nomemory_R_c (double *profil, double *weights, int *nbi, int *Kmax_, double *mini, double *maxi, double *cout_n){
	colibri_sn_weight_nomemory_c (profil, weights, nbi, Kmax_, mini, maxi, cout_n);
  }

}

