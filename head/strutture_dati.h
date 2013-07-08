#ifndef __STRUTTRE_DATI_
#define __STRUTTRE_DATI_

typedef struct lumen{
	double dim; 	// dimensione
	double c_mg; 	// concentrazione magnesio
	double c_k;	// concentrazione potassio
	double c_cl;	// concentrazione cloro
	double H;		// concentrazione idrogeno
} lumen;

typedef struct membrana{
	double I_luce; // intensita luce
} membrana;

typedef struct stroma{
	double dim; 	// dimensione
	double c_mg; 	// concentrazione magnesio
	double c_k;	// concentrazione potassio
	double c_cl;	// concentrazione cloro
} stroma;

typedef struct tilacoide{
	stroma* st_prec; // stroma precedente
	membrana* me1; // prima membrana
	lumen* lu;
	membrana* me2; // seconda membrana
	stroma* st_succ;
	double i_e;
	double i_h;
} tilacoide;

typedef struct grano{
	unsigned N; // numero tilacoidi
	tilacoide* tilacoidi_figli; // array di tilacoidi
} grano;

typedef struct lamella{
	grano* grano1;
	grano* grano2;
	double ener_ATP;
	double ener_NADPH;
} lamella;

typedef struct rete{
	unsigned N_gr; // numero grani
	grano* grani_figli; // array di struct grano
	unsigned N_lm; // numero lamelle
	lamella* lamelle_figlie; // array di struct lamella
} rete;

#endif
