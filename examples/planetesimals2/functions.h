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
void legend(char* planetdir, char* legenddir, struct reb_simulation* r, double tmax, double m_planetesimal, double total_planetesimal_mass, double inner, double outer, double powerlaw, double mp, double a, double e, double Ms, double drh);

double calc_dt(struct reb_simulation* r, double mp, double Ms, double a, double dRHill);

void calc_ELtot(double* Etot, double* Ltot, double planetesimal_mass, struct reb_simulation* r);

void calc_ae(double* a, double* e, double* d, struct reb_simulation* r, int i);

void planetesimal_forces(struct reb_simulation *a);

int check_for_encounter(struct reb_simulation* const r);

void ini_mini(struct reb_simulation* const r, struct reb_simulation* s);

void update_mini(struct reb_simulation* const r, struct reb_simulation* const s, int encounter_index);

void update_global(struct reb_simulation* const s, struct reb_simulation* const r, int encounter_index);

void clock_finish(clock_t timer, int N_encounters, char* legenddir);

//EXTERNAL VARIABLES******************************
extern double planetesimal_mass;



#endif /* defined(____functions__) */
