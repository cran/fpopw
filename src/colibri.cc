
#include "colibri.h"
// op +fp version for the L2 loss
// profil: represent a vector (double) of length nbi
// lambda: is the value (double) of the penalty (recall that the penalty i K lambda where K is the number of segment)
// mini and maxi: are the minimum and maximum value of mu during the function pruning. They should be equal to the min and max value of the vector profil)
// origine : will be updated by the function and contain at index i the position of the last breakpoint of the best segmentation arriving at position 
// cout_n will be updated by the function and contain at index i the cost of the best segmentation arriving at position i
void colibri_op_c (const double *profil, const int *nbi, const double *lambda_, const double *mini, const double *maxi, int *origine,
double *cout_n)
{
	int nb=*nbi;
	double lambda = *lambda_;
	double min=*mini;
	double max=*maxi;

	int minPosition=-1;
	double minCurrent=-10.0;
	
        /* Initialisation of object once and for all */
	//Polynome2 * p1;
	Liste * l1;  
	//Polynome2 * pTest;

	Polynome2 **stock= new Polynome2* [nb]; 
	for ( int t =0; t < nb; t++ ) stock[t]= new Polynome2();
		
         
	/* Parametrization of the first candidate segmentation (i.e 1 segment from ]0, 1] the min cost is (0 + lambda)  */
	stock[0]->reset(1.0, -2*profil[0], lambda,  -10);
	stock[0]->setStatus(2);
	
	l1 = new Liste(max, min, stock[0]);

	l1->computeMinOrMax(&minCurrent, &minPosition);
	cout_n[0] = minCurrent + lambda;
	origine[0] = minPosition;
        /* For any new data point t do 1), 2) and 3) */
	for ( int t =1; t < nb; t++ ){
        	 /* Slide 1 and Prune */
		 l1->computeRoots(cout_n[t-1]);
		 stock[t]->reset(0.0, 0.0, cout_n[t-1],  t);
		 l1->resetAllBorders(stock[t]);
		 l1->checkForDoublon();
		 l1->add(1.0, -2*profil[t], 0.0);

		 /* Compute Min */
		 l1->computeMinOrMax(&minCurrent, &minPosition);
		 cout_n[t]=minCurrent + lambda;
		 origine[t] = minPosition;	
	}
	  
	/* free stock */
	for ( int t =0; t < nb; t++ ) delete(stock[t]);
	delete[] stock;  
}

void colibri_op_weight_c (const double *profil, double *weights, const int *nbi, const double *lambda_, const double *mini, const double *maxi, int *origine,
double *cout_n)
{
	int nb=*nbi;
	double lambda = *lambda_;
	double min=*mini;
	double max=*maxi;

	int minPosition=-1;
	double minCurrent=-10.0;
	
        /* Initialisation of object once and for all */
	//Polynome2 * p1;
	Liste * l1;  
	//Polynome2 * pTest;

	Polynome2 **stock= new Polynome2* [nb]; 
	for ( int t =0; t < nb; t++ ) stock[t]= new Polynome2();
		
         
	/* Parametrization of the first candidate segmentation (i.e 1 segment from ]0, 1] the min cost is (0 + lambda)  */
	stock[0]->reset(weights[0], -2*profil[0]*weights[0], lambda,  -10);
	stock[0]->setStatus(2);
	
	l1 = new Liste(max, min, stock[0]);

	l1->computeMinOrMax(&minCurrent, &minPosition);
	cout_n[0] = minCurrent + lambda;
	origine[0] = minPosition;
        /* For any new data point t do 1), 2) and 3) */
	for ( int t =1; t < nb; t++ ){
        	 /* Slide 1 and Prune */
		 l1->computeRoots(cout_n[t-1]);
		 stock[t]->reset(0.0, 0.0, cout_n[t-1],  t);
		 l1->resetAllBorders(stock[t]);
		 l1->checkForDoublon();
		 l1->add(weights[t], -2*profil[t]*weights[t], 0.0);

		 /* Compute Min */
		 l1->computeMinOrMax(&minCurrent, &minPosition);
		 cout_n[t]=minCurrent + lambda;
		 origine[t] = minPosition;	
	}
	  
	/* free stock */
	for ( int t =0; t < nb; t++ ) delete(stock[t]);
	delete [] stock;  
}

////////////////////
////////////////////

void colibri_sn_c (const double *profil, const int *nbi, const int *Kmaxi, const double *mini, const double *maxi, int *origine,
double *cout_n, double *allCost)
{
	int nb=*nbi;
	int Kmax=*Kmaxi;
	double min=*mini;
	double max=*maxi;
	double *minCostBefore = new double[nb];
	double *minCostCurrent = new double[nb];
	double *tmp; //1
	int minPosition;
	double minCurrent;
	//int * origine = (int *) malloc(nb * sizeof(int));
	int i = 0;
    int i2 = 0;
	double somme = 0;
	double sommeC = 0;
	int turn = 1;
	//char c = 13;

    /* Initialisation Cout en 1 segment */
    while(i < nb)
	{
		somme = somme + profil[i];
		sommeC = sommeC + profil[i]*profil[i];
		minCostBefore[i] = sommeC - pow(somme, 2) / (i+1);
		origine[i]=0;
		allCost[i] = minCostBefore[i];
		i++;
	}
	/* Save */
    cout_n[0] = minCostBefore[nb-1];


    /* Initialisation Polynome Cost */
	//Polynome2 * p1;
	Liste * l1;  
	//Polynome2 * pTest;

	Polynome2 **stock= new Polynome2* [nb]; 

    i=0;
	while(i < nb)
	{
		stock[i]=new Polynome2();
		i++;	
	}


    /* Boucle turn 1 -> Kmax -1 */
	while( turn < Kmax)
	{
	  /* Print turn / Kmax */
	  /*fprintf(stderr, "%c Turn :   %d  / %d  ", c, turn, Kmax);*/
	  /* initalisation */
	  i= turn;
      i2= turn+ turn*nb;
	  stock[i]->reset(1.0, -2*profil[i], profil[i]*profil[i]+minCostBefore[turn -1],  turn);
	  stock[i]->setStatus(2);
	  l1 = new Liste(max, min, stock[i]);
	  /* Min */
	  l1->computeMinOrMax(&minCurrent, &minPosition);
	  minCostCurrent[i]=minCurrent;
	  origine[i2] = i;
	  allCost[i2] = minCurrent;

      /* iterate */
      i++;
      i2++;
	  while(i < nb)
		{
		 /* Slide 1 and Prune */
		 l1->computeRoots(minCostBefore[i-1]);
		 stock[i]->reset(0.0, 0.0, minCostBefore[i-1],  i);
		 l1->resetAllBorders(stock[i]);
		 l1->checkForDoublon();
		 l1->add(1.0, -2*profil[i], profil[i]*profil[i]);

		 /* Compute Min */
		 l1->computeMinOrMax(&minCurrent, &minPosition);
		 minCostCurrent[i]=minCurrent;
		 origine[i2] = minPosition;
		 allCost[i2] = minCurrent;
		 /* iterate */
		 i++;	
         i2++;
	  	}

	  /* Save */
      cout_n[turn] = minCostCurrent[nb-1];
	  
	  /* */
	  tmp=minCostCurrent;
	  minCostCurrent=minCostBefore;
	  minCostBefore=tmp;
	
	
	  //delete(l1);
	  /* iterate */
	  turn++;

	}
	
	/* Free All */
	/* free stock */
	i=0;
	while(i < nb)
	{
	    delete(stock[i]);	
		i++;
	}
	delete [] stock;  
	delete [] minCostBefore;
	delete [] minCostCurrent;
	//delete(origine);
	//std::cout << std::endl;

	/* Create matrix with Breakpoints positions for 0, ..., Kmax Breakpoints */
	//traceback(fileOutInt, OutPath, nb, Kmax);

    //return 0;
}
///////////////////////////////////////


void colibri_sn_weight_c (const double *profil, double *weights, const int *nbi, const int *Kmaxi, const double *mini, const double *maxi, int *origine,
double *cout_n, double *allCost)
{
	int nb=*nbi;
	int Kmax=*Kmaxi;
	double min=*mini;
	double max=*maxi;
	double *minCostBefore = new double[nb];
	double *minCostCurrent = new double[nb];
	double *tmp; //1
	int minPosition;
	double minCurrent;
	//int * origine = (int *) malloc(nb * sizeof(int));
	int i = 0;
    int i2 = 0;
	double somme = 0;
	double sommeC = 0;
	double sommeW = 0;
	int turn = 1;
	//char c = 13;

    /* Initialisation Cout en 1 segment */
    while(i < nb)
	{
		somme = somme + weights[i]*profil[i];
		sommeC = sommeC + weights[i]*profil[i]*profil[i];
		sommeW = sommeW + weights[i];
// w_i X_i^2 -2 (\sum w_i) X_b^2 + (\sum w_i) X_b^2
		minCostBefore[i] = sommeC - pow(somme, 2) / sommeW;
		origine[i]=0;
		allCost[i] = minCostBefore[i];
		i++;
	}
	/* Save */
    cout_n[0] = minCostBefore[nb-1];


    /* Initialisation Polynome Cost */
	//Polynome2 * p1;
	Liste * l1;  
	//Polynome2 * pTest;

	Polynome2 **stock= new Polynome2* [nb]; 

    i=0;
	while(i < nb)
	{
		stock[i]=new Polynome2();
		i++;	
	}


    /* Boucle turn 1 -> Kmax -1 */
	while( turn < Kmax)
	{
	  /* Print turn / Kmax */
	  /*fprintf(stderr, "%c Turn :   %d  / %d  ", c, turn, Kmax);*/
	  /* initalisation */
	  i= turn;
      i2= turn+ turn*nb;
	  stock[i]->reset(weights[i], -2*weights[i]*profil[i], weights[i]*profil[i]*profil[i]+ minCostBefore[turn -1],  turn);
	  stock[i]->setStatus(2);
	  l1 = new Liste(max, min, stock[i]);
	  /* Min */
	  l1->computeMinOrMax(&minCurrent, &minPosition);
	  minCostCurrent[i]=minCurrent;
	  origine[i2] = i;
	  allCost[i2] = minCurrent;

      /* iterate */
      i++;
      i2++;
	  while(i < nb)
		{
		 /* Slide 1 and Prune */
		 l1->computeRoots(minCostBefore[i-1]);
		 stock[i]->reset(0.0, 0.0, minCostBefore[i-1],  i);
		 l1->resetAllBorders(stock[i]);
		 l1->checkForDoublon();
		 l1->add(weights[i], -2*weights[i]*profil[i], weights[i]*profil[i]*profil[i]);

		 /* Compute Min */
		 l1->computeMinOrMax(&minCurrent, &minPosition);
		 minCostCurrent[i]=minCurrent;
		 origine[i2] = minPosition;
		 allCost[i2] = minCurrent;
		
		 /* iterate */
		 i++;	
         i2++;
	  	}

	  /* Save */
      cout_n[turn] = minCostCurrent[nb-1];
	  
	  /* */
	  tmp=minCostCurrent;
	  minCostCurrent=minCostBefore;
	  minCostBefore=tmp;
	
	
	  //delete(l1);
	  /* iterate */
	  turn++;

	}
	
	/* Free All */
	/* free stock */
	i=0;
	while(i < nb)
	{
	    delete(stock[i]);	
		i++;
	}
	delete [] stock;  
	delete [] minCostBefore;
	delete [] minCostCurrent;
	//delete(origine);
	//std::cout << std::endl;

	/* Create matrix with Breakpoints positions for 0, ..., Kmax Breakpoints */
	//traceback(fileOutInt, OutPath, nb, Kmax);

    //return 0;
}
///////////////////////////////////////



void colibri_sn_weight_nomemory_c (const double *profil, double *weights, const int *nbi, const int *Kmaxi, const double *mini, const double *maxi, double *cout_n) //, double *allCost)
{
	int nb=*nbi;
	int Kmax=*Kmaxi;
	double min=*mini;
	double max=*maxi;
	double *minCostBefore = new double[nb];
	double *minCostCurrent = new double[nb];
	double *tmp; //1
	int minPosition;
	double minCurrent;
	//int * origine = (int *) malloc(nb * sizeof(int));
	int i = 0;
    int i2 = 0;
	double somme = 0;
	double sommeC = 0;
	double sommeW = 0;
	int turn = 1;
	//char c = 13;

    /* Initialisation Cout en 1 segment */
    while(i < nb)
	{
		somme = somme + weights[i]*profil[i];
		sommeC = sommeC + weights[i]*profil[i]*profil[i];
		sommeW = sommeW + weights[i];
// w_i X_i^2 -2 (\sum w_i) X_b^2 + (\sum w_i) X_b^2
		minCostBefore[i] = sommeC - pow(somme, 2) / sommeW;
		//origine[i]=0;
		//allCost[i] = minCostBefore[i];
		i++;
	}
	/* Save */
    cout_n[0] = minCostBefore[nb-1];


    /* Initialisation Polynome Cost */
	//Polynome2 * p1;
	Liste * l1;  
	//Polynome2 * pTest;

	Polynome2 **stock= new Polynome2* [nb]; 

    i=0;
	while(i < nb)
	{
		stock[i]=new Polynome2();
		i++;	
	}


    /* Boucle turn 1 -> Kmax -1 */
	while( turn < Kmax)
	{
	  /* Print turn / Kmax */
	  /*fprintf(stderr, "%c Turn :   %d  / %d  ", c, turn, Kmax);*/
	  /* initalisation */
	  i= turn;
          i2= turn+ turn*nb;
	  stock[i]->reset(weights[i], -2*weights[i]*profil[i], weights[i]*profil[i]*profil[i]+ minCostBefore[turn -1],  turn);
	  stock[i]->setStatus(2);
	  l1 = new Liste(max, min, stock[i]);
	  /* Min */
	  l1->computeMinOrMax(&minCurrent, &minPosition);
	  minCostCurrent[i]=minCurrent;
	  //origine[i2] = i;
	  //allCost[i2] = minCurrent;

      /* iterate */
      i++;
      i2++;
	  while(i < nb)
		{
		 /* Slide 1 and Prune */
		 l1->computeRoots(minCostBefore[i-1]);
		 stock[i]->reset(0.0, 0.0, minCostBefore[i-1],  i);
		 l1->resetAllBorders(stock[i]);
		 l1->checkForDoublon();
		 l1->add(weights[i], -2*weights[i]*profil[i], weights[i]*profil[i]*profil[i]);

		 /* Compute Min */
		 l1->computeMinOrMax(&minCurrent, &minPosition);
		 minCostCurrent[i]=minCurrent;
		 //origine[i2] = minPosition;
		 //allCost[i2] = minCurrent;
		
		 /* iterate */
		 i++;	
         i2++;
	  	}

	  /* Save */
      cout_n[turn] = minCostCurrent[nb-1];
	  
	  /* */
	  tmp=minCostCurrent;
	  minCostCurrent=minCostBefore;
	  minCostBefore=tmp;
	
	
	  //delete(l1);
	  /* iterate */
	  turn++;

	}
	
	/* Free All */
	/* free stock */
	i=0;
	while(i < nb)
	{
	    delete(stock[i]);	
		i++;
	}
	delete [] stock;  
	delete [] minCostBefore;
	delete [] minCostCurrent;
	//delete(origine);
	//std::cout << std::endl;

	/* Create matrix with Breakpoints positions for 0, ..., Kmax Breakpoints */
	//traceback(fileOutInt, OutPath, nb, Kmax);

    //return 0;
}
///////////////////////////////////////

