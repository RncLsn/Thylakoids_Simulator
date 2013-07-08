#ifndef __IO_UTIL__
#define __IO_UTIL__


typedef struct io_util_dati_input{

// OTTICA

// intensita luce in ingresso
double A_1; 						// ampiezza onda di luce in ingresso al sistema
double k_onda; 		// vettore d'onda	 = 2*M_PI / 700E-9	

// indici di rifrazione
double n_st; 				// stroma
double n_me; 				// membrana
double n_lu; 				// lumen

// impedenze caratteristiche
double Z_st; 				// impedenza dello stroma
double Z_me; 				// impendenza della membrana
double Z_lu; 				// impedenza del lumen

// STRUTTURA

double N_gr; 			// numero grani
double N_til; 			// numero tilacoidi per grano

double D_me; 				// dimensione membrana
double D_min; 				// minima distanza tra le membrane
double D_st;				// dimensione iniziale lumen
double D_lu;				// dimensione iniziale stroma
double d_lu; 				// coefficiente di proporzionalita per la D_lu

// MECCANICA - ELETTROSTATICA

double T; 					// temperatura
//double K_B; 				// costante di Boltzmann
double sigma; 				// carica di superficie
double C1; 				// costante per il calcolo dell'attrazione di van der Waals
double ni_p; 				// fattore di combinazione lineare delle costanti dielettriche della membrana
double epsilon_st;			// costante dielettrica stroma
//double elettrone;			// carica elettrone
//h_tagliato 				// costante di planck

// ELETTROCHIMICA

double z_mg; 		// valenza magnesio
double z_k; 		// valenza potassio
double z_cl; 		// valenza cloro

double c_mg; 				// concentrazione iniziale magnesio
double c_k; 				// concentrazione iniziale potassio
double c_cl; 				// concentrazione iniziale cloro

// fattori di proporzionalita alla variazione di H
double ro_mg;
double ro_k;
double ro_cl;
double H; 					// concentrazione  iniziale idrogeno

// costanti per i circuiti elettrici
double R_ATP;
double R_NADPH;
double R_lu; 				// resistenza del lumen
double alfa_PSII; 			// efficienza fotoelettrica PSII
double beta_PSI;			// efficienze fotovoltaica PSI
double gamma_b6f; 			// guadagno in corrente ionica nel b6f
double gamma_OEC; 			// guadagno in corrente ionica nel OEC
double eta_H; 				// rapporto incremento_H - tensione ai capie della R_lu

double teta;				// coefficiente forza repulsiva
int N_iter;					// numero iterazioni da eseguire

} io_util_dati_input;

typedef struct io_util_coppia{
	double x; 
	double y;
} io_util_coppia;


io_util_dati_input* io_util_input(char* percorso_file);

void io_util_output(io_util_coppia* array_coppie, int n, char* nome_relazione);

#endif
