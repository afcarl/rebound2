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

void calc_ELtot(double* Etot, double* Ltot, double planetesimal_mass, struct reb_simulation* r);

void calc_ae(double* a, double* e, struct reb_simulation* r);

void time_init(void);

void time_end(char* dir);

#endif /* defined(____functions__) */
