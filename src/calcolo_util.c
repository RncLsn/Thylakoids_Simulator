#include "calcolo_util.h"
#include <math.h>
#include <assert.h>



// funzione di equilibrio
double funzione_equilibrio (double x, void *parametri) {
		params *p  = (params *) parametri;
		double cost_Hamaker = p->cost_Hamaker;
		double cost_Debye_quadra = p-> cost_Debye_quadra;
		double sigma = p-> sigma;
		double D_me = p-> D_me;
		double teta = p-> teta;
			

		/* PRESSIONE DI VAN DER WAALS */

		double waals_termine_con_segno = 1.0/(x*x*x)-2.0/((x+D_me)*(x+D_me)*(x+D_me))+1.0/((x+2*D_me)*(x+2*D_me)*(x+2*D_me));
		
		double waals_termine_modulo = (waals_termine_con_segno > 0) ? waals_termine_con_segno : -waals_termine_con_segno;
		
		double pressione_waals = cost_Hamaker*waals_termine_modulo/((6*M_PI*10));
		
		double sigma_quadro = sigma * sigma;
		
		#ifdef DEBUG_EQ
		printf("sqrt(debye)=%e ham=%e\n",sqrt(cost_Debye_quadra), cost_Hamaker);
		printf("sigma_quadro=%e\n",sigma_quadro);
		#endif

		/* PRESSIONE ELETTROSTATICA */

		double esponenziale = exp(
		-sqrt(
		cost_Debye_quadra)*
		x);
		
		assert(esponenziale != NAN);
		
		#ifdef DEBUG_EQ
		printf("argomento_exp=%e, exp=%e\n",-sqrt(cost_Debye_quadra)*x,esponenziale);
		#endif

		double pressione_elettro = teta*sqrt(cost_Debye_quadra)*sigma_quadro*esponenziale ; // nuova formula
		
		#ifdef DEBUG_EQ
		printf("pressione_vdw=%e, pressione_elettro=%e, val_funzione_equlibrio=%e\n",pressione_waals,pressione_elettro, pressione_elettro -pressione_waals);
		#endif
		
		/* CONDIZIONE DI EQUILIBRIO */
		
		return pressione_elettro- pressione_waals;
}

double funzione_crescente_luce(double A_in, double t) {
	return A_in+A_in*t;
}

double funzione_decrescente_luce(double A_in, double t) {
	double risultato = A_in-A_in*t;
	return risultato > 0 ? risultato : 0;
}

double funzione_periodica_luce(double A_in, double t) {
	double coseno = cos(t);
	return A_in*(coseno > 0 ? coseno : -coseno);
}
