#include "io_util.h"

#include <unistd.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <fcntl.h>
#include <sys/types.h> 
#include <stdlib.h>


#define BUF_DIM 256

#define CARICA(param, param_str) do{   										\
		printf("carica %s\n",param_str);								\
		rewind(fd);  												\
		trovato = 0;												\
		do {														\
			nread = fscanf(fd, "%s", buf);							\
			if(nread != 0 && !strcmp(param_str, buf)) {				\
				assert( fscanf(fd, "%s", buf) != 0 );				\
				in->param = atof(buf);						\
				trovato = 1;								\
				break;												\
			}														\
		}while(nread != 0);											\
			if(!trovato) {											\
				fprintf(stderr, "errore: parametro %s mancante!\n", param_str); 	\
				exit(EXIT_SUCCESS); 								\
			}														\
	}while(0);



io_util_dati_input* io_util_input(char* percorso_file) {
	io_util_dati_input* in;
	
	char buf [BUF_DIM];
	int nread, trovato;
	
	printf("\n***************\n\nlettura input\n\n***************\n\n");
	
	in = malloc(sizeof(io_util_dati_input));
	assert(in != NULL);
	
	FILE* fd = fopen(percorso_file, "r");
	
	CARICA(A_1,"A_1") 						// ampiezza onda di luce in ingresso al sistema
	CARICA(k_onda,"k_onda") 		// vettore d'onda	 = 2*M_PI / 700E-9	
	
	// indici di rifrazione
	CARICA(n_st,"n_st") 				// stroma
	CARICA(n_me,"n_me") 				// membrana
	CARICA(n_lu,"n_lu") 				// lumen
	
	// impedenze caratteristiche
	CARICA(Z_st,"Z_st") 				// impedenza dello stroma
	CARICA(Z_me,"Z_me") 				// impendenza della membrana
	CARICA(Z_lu,"Z_lu") 				// impedenza del lumen
	
	// STRUTTURA
	
	CARICA(N_gr,"N_gr") 			// numero grani
	CARICA(N_til,"N_til") 			// numero tilacoidi per grano
	
	CARICA(D_me,"D_me") 				// dimensione membrana
	CARICA(D_min,"D_min") 				// minima distanza tra le membrane
	CARICA(D_st,"D_st")				// dimensione iniziale lumen
	CARICA(D_lu,"D_lu")				// dimensione iniziale stroma
	CARICA(d_lu,"d_lu") 				// coefficiente di proporzionalita per la D_lu
	
	// MECCANICA - ELETTROSTATICA
	
	CARICA(T,"T") 					// temperatura
	CARICA(sigma,"sigma") 				// carica di superficie
	CARICA(C1,"C1") 				// costante per il calcolo dell'attrazione di van der Waals
	CARICA(ni_p,"ni_p") 				// fattore di combinazione lineare delle costanti dielettriche della membrana
	CARICA(epsilon_st,"epsilon_st")			// costante dielettrica stroma
	//CARICA(elettrone,"")			// carica elettrone
	//h_tagliato 				// costante di planck
	
	// ELETTROCHIMICA
	
	CARICA(z_mg,"z_mg") 		// valenza magnesio
	CARICA(z_k,"z_k") 		// valenza potassio
	CARICA(z_cl,"z_cl") 		// valenza cloro
	
	CARICA(c_mg,"c_mg") 				// concentrazione iniziale magnesio
	CARICA(c_k,"c_k") 				// concentrazione iniziale potassio
	CARICA(c_cl,"c_cl") 				// concentrazione iniziale cloro
	
	// fattori di proporzionalita alla variazione di H
	CARICA(ro_mg,"ro_mg")
	CARICA(ro_k,"ro_k")
	CARICA(ro_cl,"ro_cl")
	CARICA(H,"H") 					// concentrazione  iniziale idrogeno
	
	// costanti per i circuiti elettrici
	CARICA(R_ATP,"R_ATP")
	CARICA(R_NADPH,"R_NADPH")
	CARICA(R_lu,"R_lu") 				// resistenza del lumen
	CARICA(alfa_PSII,"alfa_PSII") 			// efficienza fotoelettrica PSII
	CARICA(beta_PSI,"beta_PSI")			// efficienze fotovoltaica PSI
	CARICA(gamma_b6f,"gamma_b6f") 			// guadagno in corrente ionica nel b6f
	CARICA(gamma_OEC,"gamma_OEC") 			// guadagno in corrente ionica nel OEC
	CARICA(eta_H,"eta_H") 				// rapporto incremento_H - tensione ai capie della R_lu
	
	CARICA(teta,"teta") 
	
	CARICA(N_iter,"N_iter") 
	
	fclose(fd);
	
	return in;
}


void io_util_output(io_util_coppia* array_coppie, int n, char* nome_relazione) {
	// crea file .dat per gnuplot
	
	int fd, i;
	char buf [BUF_DIM];
	
	sprintf(buf, "./results/%s.dat", nome_relazione);
	
	fd = open(buf, O_CREAT|O_APPEND|O_WRONLY, S_IRWXU);
	assert(fd != 0);
	
	sprintf(buf, "# id processo: %d\n",getpid()); 
	
	for(i = 0; i < n; i++) {
	
		sprintf(buf, "%e\t%e\n",(array_coppie[i]).x, (array_coppie[i]).y); 
		#ifdef DEBUG_UTIL
		int wb = write(fd, buf, strlen(buf));
		printf("scritti %d byte\n",wb);
		#else
		write(fd, buf, strlen(buf));
		#endif
	
	}
	
	close(fd);
}
