/**
 * Spreading ring
 *
 * A narrow ring of collisional particles is spreading.
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "rebound.h"
#include "tools.h"
#include "../examples/planetesimals2/functions.h"

void heartbeat(struct reb_simulation* r);

double tmax = 10000, planetesimal_mass;
int n_output = 30000;
char plntdir[200] = "output/planet_", lgnddir[200] = "output/planet_";

int main(int argc, char* argv[]){
	struct reb_simulation* r = reb_create_simulation();
	// Setup constants
	r->integrator	= REB_INTEGRATOR_WHFAST;
	r->collision	= REB_COLLISION_NONE;
	r->boundary	= REB_BOUNDARY_OPEN;
	r->G 		= 1;		
	r->N_active	= 2;
	r->softening 	= 0.01;		
	r->dt 		= 1e-3;
	r->heartbeat	= heartbeat;
    
    //other constants
    double N_planetesimals = 20;
    double M_planetesimals = 3e-6; //Total Mass of all planetesimals (default = Earth mass, 3e-6)
	double boxsize = 4.8;
	reb_configure_box(r, boxsize, 1, 1, 1);

	// Setup particles
	int _N = 1000;
	// Initial conditions
	struct reb_particle star = {0};
	star.m 		= 1;
	star.r		= 0.01;
	reb_add(r, star);

    //planet 1
    double a1=1, m1=1e-5, e1=0.01;
    struct reb_particle p1 = reb_tools_init_orbit2d(r->G,star.m,m1,a1,e1,0,0);
    reb_add(r, p1);
    
    //planetesimals
    double inner = 12, outer = 4, powerlaw = 0.5;  //higher the inner number, closer to the star
    int seed = 13;          //seed was 11
    srand(seed);
    planetesimal_mass = M_planetesimals / N_planetesimals;  //mass of each planetesimal
	while(r->N<N_planetesimals + r->N_active){
		struct reb_particle pt = {0};
		double a	= reb_random_powerlaw(boxsize/outer,boxsize/inner,powerlaw);
		double phi 	= reb_random_uniform(0,2.*M_PI);
		pt.x 		= a*cos(phi);
		pt.y 		= a*sin(phi);
		pt.z 		= a*reb_random_normal(0.0001);
		double vkep 	= sqrt(r->G*star.m/a);
		pt.vx 		= -vkep * sin(phi);
		pt.vy 		= vkep * cos(phi);
		pt.m 		= 0.0001;
		pt.r 		= .3/sqrt((double)_N);
		reb_add(r, pt);
	}
    
    // The WH integrator assumes a heliocentric coordinatesystem.
    if(r->integrator != REB_INTEGRATOR_WH) reb_move_to_com(r);
    
    //Integrate!!!
	reb_integrate(r, tmax); //INFINITY is a choice
}

void heartbeat(struct reb_simulation* r){
	if (reb_output_check(r, tmax/n_output)){
        double a_p=0, e_p=0, Etot=0, Ltot=0;
        calc_ae(&a_p, &e_p, r);
        calc_ELtot(&Etot, &Ltot, planetesimal_mass, r);
        
        FILE *append;
        //append=fopen(plntdir, "a");
        append=fopen("test.txt", "a");
        fprintf(append,"%f,%.8f,%.8f,%.15f,%.15f\n",r->t,a_p,e_p,Etot,Ltot);
        fclose(append);
        
        if(r->integrator != REB_INTEGRATOR_WH) reb_move_to_com(r);
        
		reb_output_timing(r, 0);
	}
}
