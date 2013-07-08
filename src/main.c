#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_sum.h>
#include <gsl/gsl_const_mksa.h>
#include <gsl/gsl_const_num.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>
#include <unistd.h>

#include "strutture_dati.h"
#include "io_util.h"
#include "calcolo_util.h"


// variabili globali
static io_util_dati_input in;
static int cost_Hamaker_calcolata;
static double cost_Hamaker;


// prototipi funzioni
void distribuzione_luce(grano* g);
void variazioni_ioni(grano* g);
void variazioni_dimensioni(grano* g);
void salva_risultati_relazione_luce_dimensioni(double A_1, rete* r);
rete* costruisci_rete();
void distruggi_rete(rete* r);
double intensita(gsl_complex A, gsl_complex B);
void produzione_energetica(rete* r);
double calcola_Hamaker();
void calcola();
int calcola_equilibrio(double cost_Hamaker, double cost_Debye_quadra, double* risultato);
void Hamaker_sommatoria(double* risultato) ;
void distribuzione_luce(grano* g);
double dim_grano(grano* gra);
void salva_risultati_relazione_tempo_dimensioni(double t, rete* r);
void salva_risultati(double x, double y, char* relazione);

int main(int argc, char* argv[]){
	
	io_util_dati_input* in_ptr =io_util_input("./input");
	in = *in_ptr;
	free(in_ptr);
	
	calcola();
	
	return 0;
}

void calcola() {
					
	rete* r;
	int n = in.N_iter;
	int i;
	double t = 0;
	
	double A_0= in.A_1;
	
	r = costruisci_rete();
	cost_Hamaker_calcolata = 0;
	
	while(n--){
		printf("\n***************\n\niterazione %d\n\n***************\n\n", in.N_iter - n);
		//~ #if(1)
		//~ in.A_1 = funzione_decrescente_luce(A_0, t);
		//~ in.A_1 = funzione_crescente_luce(A_0, t);
		in.A_1= funzione_periodica_luce(A_0, t);
		//~ printf("t=%e A_1=%e\n",t,in.A_1);
		//~ #endif
		//in.A_1 = (A_0*cos(t) > 0) ? A_0*cos(t) : -A_0*cos(t);
		
		t+=0.1;
		
		// aggiunta la i
		for(i = 0; i < r->N_gr; i++) {	
			distribuzione_luce(r->grani_figli+i);
			variazioni_ioni(r->grani_figli+i);
			variazioni_dimensioni(r->grani_figli+i);
		}
		produzione_energetica(r);
		//~ salva_risultati_relazione_luce_dimensioni(in.A_1, r);
		salva_risultati_relazione_tempo_dimensioni(t, r);
		salva_risultati(t, in.A_1, "ingresso");
	}
	
	distruggi_rete(r);
}

void salva_risultati(double x, double y, char* relazione) {
	
	io_util_coppia* coppia = malloc(sizeof(io_util_coppia));
	assert(coppia != NULL);
	coppia->x = x;
	coppia->y = y;
	
	io_util_output(coppia, 1, relazione);
	
	free(coppia);
}

void salva_risultati_relazione_tempo_dimensioni(double t, rete* r) {
	double dim = dim_grano(r->grani_figli+0);
	
	io_util_coppia* coppia = malloc(sizeof(io_util_coppia));
	assert(coppia != NULL);
	coppia->x = t;
	coppia->y = dim;
	
	io_util_output(coppia, 1, "tempo-dimensione");
	
	free(coppia);
}

void distribuzione_luce(grano* g){
	int i;
	double eta, h, fi, C, S;
	tilacoide* tilacoide_figlio;
	
	gsl_matrix* L_1 = gsl_matrix_alloc(2, 2);
	gsl_matrix* L_2 = gsl_matrix_alloc(2, 2);
	gsl_matrix* L_3 = gsl_matrix_alloc(2, 2);
	gsl_matrix* L =   gsl_matrix_alloc(2, 2);
	
	// calcolo della matrice L
	// itero sui tilacoidi
	
	// primo stroma
	// carico dati nella matrice L_1
	tilacoide_figlio = g->tilacoidi_figli;//[0];
	eta = in.n_st * in.k_onda;
	h = tilacoide_figlio->st_succ->dim;
	fi = in.n_me * in.k_onda * h;
	C = cos(fi);
	S = sin(fi);
	
	gsl_matrix_set(L_1, 0, 0, C);
	gsl_matrix_set(L_1, 0, 1, -S/eta);
	gsl_matrix_set(L_1, 1, 0, eta*S);
	gsl_matrix_set(L_1, 1, 1, C);
	
	// itero dal secondo strato in poi
	for(i = 0; i < g->N; i++) {
		tilacoide_figlio = g->tilacoidi_figli +i;// [i];
		
		// MEMEBRANA
		eta = in.n_me * in.k_onda;
		h = in.D_me;
		fi = in.n_me * in.k_onda * h;
		C = cos(fi);
		S = sin(fi);
		
		gsl_matrix_set(L_2, 0, 0, C);
		gsl_matrix_set(L_2, 0, 1, -S/eta);
		gsl_matrix_set(L_2, 1, 0, eta*S);
		gsl_matrix_set(L_2, 1, 1, C);
		
		// prodotto L_3 <- L_1 x L_2
		gsl_matrix_set(L_3, 0, 0, gsl_matrix_get(L_1, 0, 0) * gsl_matrix_get(L_2, 0, 0) + gsl_matrix_get(L_1, 0, 1) * gsl_matrix_get(L_2, 1, 0));
		gsl_matrix_set(L_3, 0, 1, gsl_matrix_get(L_1, 0, 0) * gsl_matrix_get(L_2, 0, 1) + gsl_matrix_get(L_1, 0, 1) * gsl_matrix_get(L_2, 1, 1));
		gsl_matrix_set(L_3, 1, 0, gsl_matrix_get(L_1, 1, 0) * gsl_matrix_get(L_2, 0, 0) + gsl_matrix_get(L_1, 1, 1) * gsl_matrix_get(L_2, 1, 0));
		gsl_matrix_set(L_3, 1, 1, gsl_matrix_get(L_1, 1, 0) * gsl_matrix_get(L_2, 0, 1) + gsl_matrix_get(L_1, 0, 1) * gsl_matrix_get(L_2, 1, 1));
		
		// L_1 <- L_3
		gsl_matrix_memcpy (L_1, L_3);
		
		// LUMEN
		eta = in.n_lu * in.k_onda;
		h = tilacoide_figlio->lu->dim;
		fi = in.n_lu * in.k_onda * h;
		C = cos(fi);
		S = sin(fi);
		
		gsl_matrix_set(L_2, 0, 0, C);
		gsl_matrix_set(L_2, 0, 1, -S/eta);
		gsl_matrix_set(L_2, 1, 0, eta*S);
		gsl_matrix_set(L_2, 1, 1, C);
		
		// prodotto L_3 <- L_1 x L_2
		gsl_matrix_set(L_3, 0, 0, gsl_matrix_get(L_1, 0, 0) * gsl_matrix_get(L_2, 0, 0) + gsl_matrix_get(L_1, 0, 1) * gsl_matrix_get(L_2, 1, 0));
		gsl_matrix_set(L_3, 0, 1, gsl_matrix_get(L_1, 0, 0) * gsl_matrix_get(L_2, 0, 1) + gsl_matrix_get(L_1, 0, 1) * gsl_matrix_get(L_2, 1, 1));
		gsl_matrix_set(L_3, 1, 0, gsl_matrix_get(L_1, 1, 0) * gsl_matrix_get(L_2, 0, 0) + gsl_matrix_get(L_1, 1, 1) * gsl_matrix_get(L_2, 1, 0));
		gsl_matrix_set(L_3, 1, 1, gsl_matrix_get(L_1, 1, 0) * gsl_matrix_get(L_2, 0, 1) + gsl_matrix_get(L_1, 0, 1) * gsl_matrix_get(L_2, 1, 1));
		
		// L_1 <- L_3
		gsl_matrix_memcpy (L_1, L_3);
		
		// MEMEBRANA
		eta = in.n_me * in.k_onda;
		h = in.D_me;
		fi = in.n_me * in.k_onda * h;
		C = cos(fi);
		S = sin(fi);
		
		gsl_matrix_set(L_2, 0, 0, C);
		gsl_matrix_set(L_2, 0, 1, -S/eta);
		gsl_matrix_set(L_2, 1, 0, eta*S);
		gsl_matrix_set(L_2, 1, 1, C);
		
		// prodotto L_3 <- L_1 x L_2
		gsl_matrix_set(L_3, 0, 0, gsl_matrix_get(L_1, 0, 0) * gsl_matrix_get(L_2, 0, 0) + gsl_matrix_get(L_1, 0, 1) * gsl_matrix_get(L_2, 1, 0));
		gsl_matrix_set(L_3, 0, 1, gsl_matrix_get(L_1, 0, 0) * gsl_matrix_get(L_2, 0, 1) + gsl_matrix_get(L_1, 0, 1) * gsl_matrix_get(L_2, 1, 1));
		gsl_matrix_set(L_3, 1, 0, gsl_matrix_get(L_1, 1, 0) * gsl_matrix_get(L_2, 0, 0) + gsl_matrix_get(L_1, 1, 1) * gsl_matrix_get(L_2, 1, 0));
		gsl_matrix_set(L_3, 1, 1, gsl_matrix_get(L_1, 1, 0) * gsl_matrix_get(L_2, 0, 1) + gsl_matrix_get(L_1, 0, 1) * gsl_matrix_get(L_2, 1, 1));
		
		// L_1 <- L_3
		gsl_matrix_memcpy (L_1, L_3);
		
		// stroma
		eta = in.n_st * in.k_onda;
		h = tilacoide_figlio->st_succ->dim;
		fi = in.n_me * in.k_onda * h;
		C = cos(fi);
		S = sin(fi);
		
		gsl_matrix_set(L_2, 0, 0, C);
		gsl_matrix_set(L_2, 0, 1, -S/eta);
		gsl_matrix_set(L_2, 1, 0, eta*S);
		gsl_matrix_set(L_2, 1, 1, C);
		
		// prodotto L_3 <- L_1 x L_2
		gsl_matrix_set(L_3, 0, 0, gsl_matrix_get(L_1, 0, 0) * gsl_matrix_get(L_2, 0, 0) + gsl_matrix_get(L_1, 0, 1) * gsl_matrix_get(L_2, 1, 0));
		gsl_matrix_set(L_3, 0, 1, gsl_matrix_get(L_1, 0, 0) * gsl_matrix_get(L_2, 0, 1) + gsl_matrix_get(L_1, 0, 1) * gsl_matrix_get(L_2, 1, 1));
		gsl_matrix_set(L_3, 1, 0, gsl_matrix_get(L_1, 1, 0) * gsl_matrix_get(L_2, 0, 0) + gsl_matrix_get(L_1, 1, 1) * gsl_matrix_get(L_2, 1, 0));
		gsl_matrix_set(L_3, 1, 1, gsl_matrix_get(L_1, 1, 0) * gsl_matrix_get(L_2, 0, 1) + gsl_matrix_get(L_1, 0, 1) * gsl_matrix_get(L_2, 1, 1));
		
		// L_1 <- L_3
		gsl_matrix_memcpy (L_1, L_3);
		
		
		
		
		
	}
	
	// L <- L_3 - L CALCOLATA
	gsl_matrix_memcpy (L, L_3);
	
	// calcolo di A_n ampiezza del campo rifratto nell'ultimo mezzo
	
	double eta_1, eta_n;
	eta_1 = eta_n = in.n_st*in.k_onda;
	double c11   = gsl_matrix_get(L, 0, 0);
	double c12   = gsl_matrix_get(L, 0, 1);
	double c21   = gsl_matrix_get(L, 1, 0);
	double c22   = gsl_matrix_get(L, 1, 1);
	
	//gsl_complex A_n = gsl_complex_rect(0, 0); 
	
	gsl_complex A_n = gsl_complex_div(gsl_complex_rect(2* eta_1 * in.A_1, 0),  
										gsl_complex_rect(eta_1 * c11+ eta_n * c22, eta_1 * eta_n * c12 - c21));
	
	// output - cacolo intensita nelle membrane (a ritroso)
	gsl_complex v_2[2];
	gsl_complex v_1[2];
	
	// primo passo, calcolo v_n, corrisponde allo stroma
	v_2[0] = A_n;
	v_2[1] = gsl_complex_mul_real(A_n, eta_n);
	
	// calcolo v_{n-1}, corrisponde alla membrana
	v_1[0] = v_2[0];
	v_1[1] = v_2[1];
	
	gsl_complex P, Q;
	P = v_1[0];
	Q = v_1[1];
	
	gsl_complex A, B; // ampiezze del campo ottico
	eta = in.n_me * in.k_onda;
	A = gsl_complex_div_real(gsl_complex_add(P, gsl_complex_div_real(Q, -eta)), 2);
	B = gsl_complex_div_real(gsl_complex_add(P, gsl_complex_div_real(Q, eta)), 2);
	
	// A e B determinano il campo
	// calcolo dell'intensita media
	g->tilacoidi_figli[(g->N-1)].me2->I_luce = intensita(A, B);
	
	
	v_2[0] = v_1[0];
	v_2[1] = v_1[1];
	
	// calcolo ultimo lumen
	eta = in.n_lu * in.k_onda;
	//h = g->tilacoidi_figli[(g->N-1)]->lu->dim;
	h = g->tilacoidi_figli[(g->N-1)].lu->dim;
	fi = in.n_lu * in.k_onda * h;
	C = cos(fi);
	S = sin(fi);
	
	gsl_matrix_set(L, 0, 0, C);
	gsl_matrix_set(L, 0, 1, -S/eta);
	gsl_matrix_set(L, 1, 0, eta*S);
	gsl_matrix_set(L, 1, 1, C);
	
	// prodotto v_1 <- L x v_2
	
	v_1[0] = gsl_complex_add(gsl_complex_mul_real(v_2[0], gsl_matrix_get(L, 0, 0)), gsl_complex_mul_real(v_2[1], gsl_matrix_get(L_1, 0, 1)));
	v_1[1] = gsl_complex_add(gsl_complex_mul_real(v_2[0], gsl_matrix_get(L, 1, 0)), gsl_complex_mul_real(v_2[1], gsl_matrix_get(L_1, 1, 1)));
	
	v_2[0] = v_1[0];
	v_2[1] = v_1[1];
	
	// penultima membrana
	eta = in.n_me * in.k_onda;
	h = in.D_me;
	fi = in.n_me * in.k_onda * h;
	C = cos(fi);
	S = sin(fi);
	
	gsl_matrix_set(L_2, 0, 0, C);
	gsl_matrix_set(L_2, 0, 1, -S/eta);
	gsl_matrix_set(L_2, 1, 0, eta*S);
	gsl_matrix_set(L_2, 1, 1, C);
	
	// prodotto v_1 <- L x v_2
	v_1[0] = gsl_complex_add(gsl_complex_mul_real(v_2[0], gsl_matrix_get(L, 0, 0)), gsl_complex_mul_real(v_2[1], gsl_matrix_get(L_1, 0, 1)));
	v_1[1] = gsl_complex_add(gsl_complex_mul_real(v_2[0], gsl_matrix_get(L, 1, 0)), gsl_complex_mul_real(v_2[1], gsl_matrix_get(L_1, 1, 1)));
	
	//gsl_complex P, Q;
	P = v_1[0];
	Q = v_1[1];
	
	//gsl_complex A, B; // ampiezze del campo ottico
	eta = in.n_me * in.k_onda;
	A = gsl_complex_div_real(gsl_complex_add(P, gsl_complex_div_real(Q, -eta)), 2);
	B = gsl_complex_div_real(gsl_complex_add(P, gsl_complex_div_real(Q, eta)), 2);
	
	
	// A e B determinano il campo
	// calcolo dell'intensita media
	//g->tilacoidi_figli[(g->N-1)]->me1->I_luce = intensita(A, B);
	g->tilacoidi_figli[(g->N-1)].me1->I_luce = intensita(A, B);
	
	
	v_2[0] = v_1[0];
	v_2[1] = v_1[1];
	
	// su tutti i tilacoidi
	for(i = g->N-2; i >= 0; i--) {
		tilacoide_figlio = g->tilacoidi_figli+i;//[i];
		
		// stroma
		eta = in.n_st * in.k_onda;
		h = tilacoide_figlio->st_succ->dim;
		fi = in.n_me * in.k_onda * h;
		C = cos(fi);
		S = sin(fi);
		
		gsl_matrix_set(L, 0, 0, C);
		gsl_matrix_set(L, 0, 1, -S/eta);
		gsl_matrix_set(L, 1, 0, eta*S);
		gsl_matrix_set(L, 1, 1, C);
		
		// prodotto v_1 <- L x v_2
		
		v_1[0] = gsl_complex_add(gsl_complex_mul_real(v_2[0], gsl_matrix_get(L, 0, 0)), gsl_complex_mul_real(v_2[1], gsl_matrix_get(L_1, 0, 1)));
		v_1[1] = gsl_complex_add(gsl_complex_mul_real(v_2[0], gsl_matrix_get(L, 1, 0)), gsl_complex_mul_real(v_2[1], gsl_matrix_get(L_1, 1, 1)));
		
		v_2[0] = v_1[0];
		v_2[1] = v_1[1];
		
		//gsl_complex P, Q;
		P = v_1[0];
		Q = v_1[1];
		
		//gsl_complex A, B; // ampiezze del campo ottico
		eta = in.n_me * in.k_onda;
		A = gsl_complex_div_real(gsl_complex_add(P, gsl_complex_div_real(Q, -eta)), 2);
		B = gsl_complex_div_real(gsl_complex_add(P, gsl_complex_div_real(Q, eta)), 2);
	
		
		// A e B determinano il campo
		// calcolo dell'intensita media
		tilacoide_figlio->me2->I_luce = intensita(A, B);
		
		v_2[0] = v_1[0];
		v_2[1] = v_1[1];
		
		// calcolo ultimo lumen
		eta = in.n_lu * in.k_onda;
		h = tilacoide_figlio->lu->dim;
		fi = in.n_lu * in.k_onda * h;
		C = cos(fi);
		S = sin(fi);
		
		gsl_matrix_set(L, 0, 0, C);
		gsl_matrix_set(L, 0, 1, -S/eta);
		gsl_matrix_set(L, 1, 0, eta*S);
		gsl_matrix_set(L, 1, 1, C);
		
		// prodotto v_1 <- L x v_2
		
		v_1[0] = gsl_complex_add(gsl_complex_mul_real(v_2[0], gsl_matrix_get(L, 0, 0)), gsl_complex_mul_real(v_2[1], gsl_matrix_get(L_1, 0, 1)));
		v_1[1] = gsl_complex_add(gsl_complex_mul_real(v_2[0], gsl_matrix_get(L, 1, 0)), gsl_complex_mul_real(v_2[1], gsl_matrix_get(L_1, 1, 1)));
		
		v_2[0] = v_1[0];
		v_2[1] = v_1[1];
		
		// penultima membrana
		eta = in.n_me * in.k_onda;
		h = in.D_me;
		fi = in.n_me * in.k_onda * h;
		C = cos(fi);
		S = sin(fi);
		
		gsl_matrix_set(L_2, 0, 0, C);
		gsl_matrix_set(L_2, 0, 1, -S/eta);
		gsl_matrix_set(L_2, 1, 0, eta*S);
		gsl_matrix_set(L_2, 1, 1, C);
		
		// prodotto v_1 <- L x v_2
		v_1[0] = gsl_complex_add(gsl_complex_mul_real(v_2[0], gsl_matrix_get(L, 0, 0)), gsl_complex_mul_real(v_2[1], gsl_matrix_get(L_1, 0, 1)));
		v_1[1] = gsl_complex_add(gsl_complex_mul_real(v_2[0], gsl_matrix_get(L, 1, 0)), gsl_complex_mul_real(v_2[1], gsl_matrix_get(L_1, 1, 1)));
		
		//gsl_complex P, Q;
		P = v_1[0];
		Q = v_1[1];
		
		//gsl_complex A, B; // ampiezze del campo ottico
		eta = in.n_me * in.k_onda;
		A = gsl_complex_div_real(gsl_complex_add(P, gsl_complex_div_real(Q, -eta)), 2);
		B = gsl_complex_div_real(gsl_complex_add(P, gsl_complex_div_real(Q, eta)), 2);
	
	
		// A e B determinano il campo
		// calcolo dell'intensita media
		tilacoide_figlio->me1->I_luce = intensita(A, B);
		
		
		v_2[0] = v_1[0];
		v_2[1] = v_1[1];
	}
	
	
	// libera memoria matrici
	gsl_matrix_free(L_1);
	gsl_matrix_free(L_2);
	gsl_matrix_free(L_3);
	gsl_matrix_free(L);
}

//calcolo le correnti prodotte nel grano da ciascuna membrana, e le conseguenti variazioni di ioni
void variazioni_ioni(grano* g){
	// in ogni tilacoide, 
	int i;
	tilacoide* tilacoide_figlio;
	double I_luce;
	double i_h, i_e;
	i_h = i_e = 0;
	double delta_H;
	for(i = 0; i < g->N; i++) {
		tilacoide_figlio = g->tilacoidi_figli+i;//[i];
		// membrana 1
		I_luce = tilacoide_figlio->me1->I_luce;
		
		i_e += in.alfa_PSII * I_luce;
		i_h += (in.gamma_b6f + in.gamma_OEC) * i_e;
		
		// membrana 2
		I_luce = tilacoide_figlio->me2->I_luce;
		
		i_e += in.alfa_PSII * I_luce;
		i_h += (in.gamma_b6f + in.gamma_OEC) * i_e;
		
		// memorizzo nel tilacoide
		tilacoide_figlio->i_e = i_e;
		tilacoide_figlio->i_h = i_h;
		
		// calcolo variazione H
		delta_H = in.eta_H * i_h * in.R_lu;
		tilacoide_figlio->lu->H = in.H + delta_H; // base + variazione
		
		// calcolo variazioni ioni
		// lumen
		tilacoide_figlio->lu->c_mg = in.c_mg - in.ro_mg * delta_H; 		// concentrazione magnesio
		tilacoide_figlio->lu->c_k = in. c_k - in.ro_k * delta_H;		// concentrazione potassio
		tilacoide_figlio->lu->c_cl = in.c_cl + in.ro_cl * delta_H ;		// concentrazione cloro
		
		// tralasciamo il primo e l'ultimo stroma
		if(i > 0){
			// stroma (variazioni di segno opposto), dipende dalla media dei delta_H
			//double delta_H_medio = (delta_H + (g->tilacoidi_figli[i-1])->lu->H) / 2;
			double delta_H_medio = (delta_H + (g->tilacoidi_figli[i-1]).lu->H - in.H) / 2;
			tilacoide_figlio->st_prec->c_mg = in.c_mg + in.ro_mg * delta_H_medio; 		// concentrazione magnesio
			tilacoide_figlio->st_prec->c_k = in. c_k + in.ro_k * delta_H_medio;		// concentrazione potassio
			tilacoide_figlio->st_prec->c_cl = in.c_cl - in.ro_cl * delta_H_medio ;		// concentrazione cloro
		}
	}
}

void variazioni_dimensioni(grano* g){		
	// calcolo equilibrio... risultato: soluzione equazione equilibrio pressioni, oppure distanza minima
	double cost_debye_quadra;
	int i, status;
	tilacoide* tilacoide_figlio;
	double D_st_eq;	// equilibrio
	
	// COSTANTE DI HAMAKER
	
	// calcolo da eseguire solo una volta
	if(!cost_Hamaker_calcolata) {
		cost_Hamaker = calcola_Hamaker();
		
		#if defined(DEBUG) || defined(DEBUG_EQ)
		printf("cost_hamaker=%e\n", cost_Hamaker);
		#endif
			
		cost_Hamaker_calcolata = 1;
	}
	
	
	
	//per ogni tilacoide
	for(i = 0; i < g->N; i++) {
		tilacoide_figlio = g->tilacoidi_figli+i;//[i];
	
		// calcolo stroma successivo, tranne ultimo tilacoide
		if(i < g->N-1) {
			stroma* st = tilacoide_figlio->st_succ;
			
			// COSTANTE DI DEBYE
			
			double termine_sommatoria = ((in.H) + pow((in.z_mg),2)*(st->c_mg) + pow((in.z_cl),2) * (st->c_cl) + pow((in.z_k),2)*(st->c_k))*GSL_CONST_NUM_AVOGADRO;
			
			#ifdef DEBUG_DEBYE
			printf("termine_sommatoria=%e\n",termine_sommatoria);
			sqrt(termine_sommatoria);
			#endif
			
			double termine_denominatore = in.epsilon_st *  GSL_CONST_MKSA_BOLTZMANN * (in.T);
			
			#ifdef DEBUG_DEBYE
			printf("termine_denominatore=%e\n",termine_denominatore);
			sqrt(termine_denominatore);
			#endif
			
			double elettrone_quadro = GSL_CONST_MKSA_ELECTRON_CHARGE*GSL_CONST_MKSA_ELECTRON_CHARGE;
			
			#ifdef DEBUG_DEBYE
			printf("elettrone_quadro=%e\n",elettrone_quadro);
			sqrt(termine_denominatore);
			#endif
			
			cost_debye_quadra = elettrone_quadro * termine_sommatoria / termine_denominatore;
			
			#ifdef DEBUG_DEBYE
			printf("debye_quadra=%e\n", cost_debye_quadra);
			sqrt(termine_denominatore);
			#endif
			
			/*
			cost_debye_quadra = (1.0/(in.epsilon_st * GSL_CONST_MKSA_BOLTZMANN * (in.T)))*pow(GSL_CONST_MKSA_ELECTRON_CHARGE,2)* ((in.H) + pow((in.z_mg),2) * (st->c_mg) + pow((in.z_cl),2) * (st->c_cl) + pow((in.z_k),2)*(st->c_k))*GSL_CONST_NUM_AVOGADRO;
			*/
			
			//cost_debye_quadra = 9E18; // <----<< costante! valore medio
			
			
			
			#if defined(DEBUG) || defined(DEBUG_EQ)
			printf("debye_quadra=%e\n", cost_debye_quadra);
			#endif
			
			// SOLUZIONE EQUILIBRIO
			
			status = calcola_equilibrio(cost_Hamaker, cost_debye_quadra, &D_st_eq);
			assert(status == GSL_SUCCESS);
			
			#if defined(DEBUG) || defined(DEBUG_EQ)
			printf("D_st_eq=%e\n", D_st_eq);
			#endif
			
			
			// memorizza nuova dimensione stroma successivo al tilacoide
			st->dim = D_st_eq;
		}
		
		// calcolo lumen
		if(i == g->N-1) {
			// ultimo tilacoide lumen prop allo stroma precedente
			tilacoide_figlio->lu->dim = (in.d_lu) * (tilacoide_figlio->st_prec->dim);
			
		}
		else if(i == 0) {
			// primo tilacoide, prop all successivo
			tilacoide_figlio->lu->dim = (in.d_lu) * (tilacoide_figlio->st_succ->dim);
		}
		else {
			// caso generico, prop alla media
			tilacoide_figlio->lu->dim = (in.d_lu) * (tilacoide_figlio->st_prec->dim + tilacoide_figlio->st_succ->dim) / 2;
		}
	}
}

int calcola_equilibrio(double cost_Hamaker, double cost_Debye_quadra, double* risultato) {
	// da esempio: http://www.gnu.org/software/gsl/manual/html_node/Root-Finding-Examples.html
	int status;
    int iter = 0, max_iter = 100;
	const gsl_root_fsolver_type *T;
	gsl_root_fsolver *s;
	double r = 0;
	double x_lo = 1E-10, x_hi = 8E-9;
	gsl_function F;
	params parametri = {cost_Hamaker, cost_Debye_quadra, in.sigma, in.D_me, in.teta};
	
	F.function = &funzione_equilibrio;
	F.params = &parametri;
	
	//T = gsl_root_fsolver_brent;
	T = gsl_root_fsolver_bisection;
	s = gsl_root_fsolver_alloc (T);
	gsl_root_fsolver_set (s, &F, x_lo, x_hi);
	
	#ifdef DEBUG_ROOT
	printf ("using %s method\n", 
		   gsl_root_fsolver_name (s));
	
	printf ("%5s [%12s, %12s] %12s %12s\n",
		   "iter", "lower", "upper", "root", "err(est)");
	#endif
	
	do
	 {
	   iter++;
	   status = gsl_root_fsolver_iterate (s);
	   r = gsl_root_fsolver_root (s);
	   x_lo = gsl_root_fsolver_x_lower (s);
	   x_hi = gsl_root_fsolver_x_upper (s);
	   status = gsl_root_test_interval (x_lo, x_hi,
										0, 0.001);
										
		// controllo su dimensione minima
		if(x_hi < in.D_min){
			*risultato = in.D_min;
			#ifdef DEBUG_ROOT
			printf ("limite superiore= %e < %e=minimo consentito\n",x_hi,in.D_min);
			#endif
			gsl_root_fsolver_free (s);
			return GSL_SUCCESS;
		}
		
		#ifdef DEBUG_ROOT
		if (status == GSL_SUCCESS)
		 printf ("Converged:\n");

		
		printf ("%5d [%e, %e] %e %e\n",
                   iter, x_lo, x_hi,
                   r, 
                   x_hi - x_lo);
		#endif
	 }
	while (status == GSL_CONTINUE && iter < max_iter);
	
	gsl_root_fsolver_free (s);
	
	*risultato = r;
	return status;
}

void salva_risultati_relazione_luce_dimensioni(double A_1, rete* r){
	// salva le coppie input - output
	
	double dim = dim_grano(r->grani_figli+0);
	
	io_util_coppia* coppia = malloc(sizeof(io_util_coppia));
	assert(coppia != NULL);
	coppia->x = A_1;
	coppia->y = dim;
	
	io_util_output(coppia, 1, "luce-dimensione");
	
	free(coppia);
}

double dim_grano(grano* gra){
		double dim; // ordinata
		int i;
		
		// le membrane non variano
		dim = gra->N*2*in.D_me;
		
		for(i = 0; i < gra->N; i++) {
			// sequenza: lumen - stroma
			
			dim+=gra->tilacoidi_figli[i].lu->dim;
			
			// stroma, salto l'ultimo
			if(i <= gra->N -1){
				dim+=gra->tilacoidi_figli[i].st_succ->dim;
				#ifdef DEBUG
				printf("dim_stroma_%d=%e\n",i, gra->tilacoidi_figli[i].st_succ->dim);
				#endif
			}
		}
	return dim;
}

// costruisce la rete di tilacoidi, la inizializza con i dati della struct input in
rete* costruisci_rete() {
	rete* r;
	int i,j, k;
	grano* grano_figlio;
	tilacoide* tilacoide_figlio;
	stroma* stroma_figlio;
	lumen* lumen_figlio;
	
	// crea rete
	
	r = malloc(sizeof(rete));
	assert(r != NULL);
	
	// crea grani
	r->N_gr = in.N_gr;
	r->grani_figli = malloc(sizeof(grano) * in.N_gr);
	assert(r->grani_figli != NULL);
	
	// inizializza grani
	for(i = 0; i < r->N_gr; i++) {
		grano_figlio = r->grani_figli + i;
		grano_figlio->N = in.N_til;	
		grano_figlio->tilacoidi_figli = malloc(sizeof(tilacoide) * grano_figlio->N);
		assert(grano_figlio->tilacoidi_figli != NULL);
		
		// crea tilacoidi
		for(j = 0; j < grano_figlio->N; j++) {
			tilacoide_figlio = grano_figlio->tilacoidi_figli +j;
			
			if(j == 0){
				// caso particolare
				tilacoide_figlio->st_prec = malloc(sizeof(stroma));
				assert(tilacoide_figlio->st_prec != NULL);
				
				
				//inizializzazione
				// stroma precedente
				stroma_figlio = tilacoide_figlio->st_prec;
				stroma_figlio->dim = in.D_st; 				// dimensione
				stroma_figlio->c_mg = in.c_mg; 				// concentrazione magnesio
				stroma_figlio->c_k = in.c_k;					// concentrazione potassio
				stroma_figlio->c_cl = in.c_cl;				// concentrazione cloro
				
				
			}
			else {
				tilacoide_figlio->st_prec = ((grano_figlio->tilacoidi_figli)+(j-1))->st_succ; // il successivo del precedente diventa il precedente del corrente 
			}
			tilacoide_figlio->me1 = malloc(sizeof(membrana));
			assert(tilacoide_figlio->me1 != NULL);
			tilacoide_figlio->lu = malloc(sizeof(lumen));
			assert(tilacoide_figlio->lu != NULL);
			tilacoide_figlio->me2 = malloc(sizeof(membrana));
			assert(tilacoide_figlio->me2 != NULL);
			tilacoide_figlio->st_succ = malloc(sizeof(stroma)); // condiviso con il successivo 
			assert(tilacoide_figlio->st_succ != NULL);
			
			//inizializzazione
			
			// lumen
			lumen_figlio = tilacoide_figlio->lu;
			lumen_figlio->dim = in.D_lu; 				// dimensione
			lumen_figlio->c_mg = in.c_mg; 				// concentrazione magnesio
			lumen_figlio->c_k = in.c_k;					// concentrazione potassio
			lumen_figlio->c_cl = in.c_cl;				// concentrazione cloro
			lumen_figlio->H = in.H;						// concentrazione idrogen
			
			// stroma successivo
			stroma_figlio = tilacoide_figlio->st_succ;
			stroma_figlio->dim = in.D_st; 				// dimensione
			stroma_figlio->c_mg = in.c_mg; 				// concentrazione magnesio
			stroma_figlio->c_k = in.c_k;					// concentrazione potassio
			stroma_figlio->c_cl = in.c_cl;				// concentrazione cloro
			
		}
	}

	
	// crea lamelle
	// collegano tutti i grani con tutti i grani, compresi se stessi
	r->N_lm = r->N_gr * r->N_gr; // le lamelle sono il quadrato dei grani
	r->lamelle_figlie = malloc(sizeof(lamella) * r->N_lm);
	assert(r->lamelle_figlie != NULL);
	k = 0;
	for(i = 0; i < r->N_gr; i++){
		for(j = 0; j < r->N_gr; j++) {
			r->lamelle_figlie[k].grano1 = r->grani_figli +i;
			r->lamelle_figlie[k].grano2 = r->grani_figli +j;
		}
	}
	return r;
}

void distruggi_rete(rete* r) {
	int i, j;
	
	tilacoide* til;
	grano* gra;
	
	// per ogni grano
	for(i = 0; i < r->N_gr; i++) {
		// grano corrente
		gra = r->grani_figli + i; 
		// per ogni tilacoide
		for(j = 0; j < gra->N; j++) {
				// tilacoide corrente
				til = gra->tilacoidi_figli + j;
				if(j == 0) {
					free(til->st_prec);
				}
				free(til->me1);
				free(til->lu);
				free(til->me2);
				free(til->st_succ);
			
		}
		// libero tutti i tilacoidi del grano
		free(gra->tilacoidi_figli);
	}
	// libero tutti i grani della rete
	free(r->grani_figli);
	// libero tutte le lamelle della rete
	free(r->lamelle_figlie);
	// libera la rete
	free(r);
}

double intensita(gsl_complex A, gsl_complex B) {
	// PROVVISORIO
	//~ double intensita = (gsl_complex_abs2(A) + gsl_complex_abs2(B) )/(2*in.Z_st);
	double intensita = gsl_complex_abs(gsl_complex_mul(gsl_complex_add(A,B) , gsl_complex_conjugate(gsl_complex_add(A,B))));
	printf("I=%e\n",intensita);
	#ifdef DEBUG
	printf("|A|=%e\n",gsl_complex_abs(A));
	printf("|B|=%e\n",gsl_complex_abs(B));
	printf("intensita=%e\n",intensita);
	#endif
	return intensita;
}

double calcola_Hamaker() {
	double cost_Hamaker, risultato_sommatoria;
	
	Hamaker_sommatoria(&risultato_sommatoria);
	
	cost_Hamaker = 3* GSL_CONST_MKSA_BOLTZMANN *(in.T) / 2 * risultato_sommatoria;
	
	return cost_Hamaker;
}

void Hamaker_sommatoria(double* risultato) {
	// da: http://www.gnu.org/software/gsl/manual/html_node/Example-of-accelerating-a-series.html
	double epsilon_st, epsilon_me, epsilon_p, epsilon_h;
	int xi, xi_n;
	
	xi = (int) (2* M_PI * GSL_CONST_MKSA_BOLTZMANN * (in.T) / GSL_CONST_MKSA_PLANCKS_CONSTANT_HBAR);
	
	int numero_iterazioni = 100;
	double t[numero_iterazioni];
	double sum_accel, err;
	//double sum = 0;
	int n;
	
	gsl_sum_levin_u_workspace * w 
	 = gsl_sum_levin_u_alloc (numero_iterazioni);
	
	for (n = 0; n < numero_iterazioni; n++)
	 {
		 xi_n = xi * n;
		 // assumo le parentesi quadre della formula come parti intere
		
		epsilon_st = (74.8/(int)(1+xi_n/6.5E-5)) + (1.464/(int)(1+pow((xi_n/0.021), 2)))+ (0.737/(int)(1+pow((xi_n/0.067),2)))+ (0.153/(int) pow(1+(xi_n/0.092),2))+(0.014/(int)(1+pow((xi_n/0.2),2)))+(0.075/(int)(1+pow((xi_n/0.42),2)))+(0.078/(int)(1+pow((xi_n/12.7),2)))+ 1;
		
		
		epsilon_p = (1+(in.C1)/(int)(1+(xi_n/10)*(xi_n/10)));
		epsilon_h = (1+1/(int)(1+(xi_n/10.4)*(xi_n/10.4)));
		
		epsilon_me = (in.ni_p) * epsilon_p + (1-(in.ni_p)) * epsilon_h;
		
		t[n] = pow((epsilon_st-epsilon_me)/(epsilon_st+epsilon_me), 2);
	  }
	
	gsl_sum_levin_u_accel (t, numero_iterazioni, w, &sum_accel, &err);
	
	gsl_sum_levin_u_free (w);
	
	*risultato = sum_accel;
}

void produzione_energetica(rete* r){
		//DA FARE
}
