//
//  functions.h
//  
//
//  Created by Ari Silburt on 2015-06-12.
//
//

#ifndef ____functions__
#define ____functions__

#include <stdio.h>
#include "../../src/rebound.h"

//FUNCTIONS******************************
void legend(char* planetdir, char* legenddir, char* charizard, struct reb_simulation* r, double tmax, double m_planetesimal, double total_planetesimal_mass, int N_planetesimals, double inner, double outer, double powerlaw, double mp, double a, double e, double Ms, double drh, int HYBRID_ON);

double calc_dt(struct reb_simulation* r, double mp, double Ms, double a, double dRHill);

void calc_Hill2(struct reb_simulation* r);

double calc_Etot(struct reb_simulation* a, double* K1, double* U1);

//void calc_ELtot(double* Etot, double* Ktot, double* Utot, double* Ltot, double planetesimal_mass, struct reb_simulation* r);

void calc_ae(double* a, double* e, double* d, struct reb_simulation* r, int i, double t);

void planetesimal_forces_routine();

void planetesimal_forces_global(struct reb_simulation *a);

void planetesimal_forces_mini(struct reb_simulation *a);

void check_for_encounter(struct reb_simulation* const r, int* N_encounters, double* minimum_ratio);

void ini_mini(struct reb_simulation* const r, struct reb_simulation* s, int turn_planetesimal_forces_on);

void update_global(struct reb_simulation* const s, struct reb_simulation* const r, int N_encounters_previous, int N_encounters);

void add_or_subtract_particles(struct reb_simulation* r, struct reb_simulation* s, int N_encounters,int N_encounters_previous, int dN);

void update_previous_global_positions(struct reb_simulation* r, int N_encounters);

void update_encounter_indices(int* N_encounters, int* N_encounters_previous);

void clock_finish(clock_t timer, int N_encounters, char* legenddir);

void global_free();


//EXTERNAL VARIABLES******************************
extern double planetesimal_mass;
extern int* encounter_index;
extern int* previous_encounter_index;
extern double* Hill2;
extern double* x_prev; extern double* y_prev; extern double* z_prev;
extern double t_prev;
extern int N_encounters_tot;
extern struct reb_simulation* r;
extern struct reb_simulation* s;


#endif /* defined(____functions__) */
