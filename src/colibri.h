#include "liste.h"


void colibri_op_c (const double *profil, const int *nbi, const double *lambda_, const double *mini, const double *maxi, int *origine,
double *cout_n);

void colibri_op_weight_c (const double *profil, double *weights, const int *nbi, const double *lambda_, const double *mini, const double *maxi, int *origine,
double *cout_n);

void colibri_sn_c (const double *profil, const int *nbi, const int *Kmaxi, const double *mini, const double *maxi, int *origine,
double *cout_n, double *allCost);

void colibri_sn_weight_c (const double *profil, double *weights, const int *nbi, const int *Kmaxi, const double *mini, const double *maxi, int *origine,
double *cout_n, double *allCost);

void colibri_sn_weight_nomemory_c (const double *profil, double *weights, const int *nbi, const int *Kmaxi, const double *mini, const double *maxi, double *cout_n);
