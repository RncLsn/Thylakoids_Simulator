#ifndef __CALCOLO_UTIL__
#define __CALCOLO_UTIL__

typedef struct params {
	double cost_Hamaker, cost_Debye_quadra, sigma, D_me, teta;
} params;


double funzione_equilibrio (double x, void *parametri);

// funzioni del tempo per la luce in ingresso
double funzione_crescete_luce(double A_in, double t);
double funzione_decrescente_luce(double A_in, double t);
double funzione_periodica_luce(double A_in, double t);

#endif
