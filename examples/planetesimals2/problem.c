/**
 * A.S. This is my planetesimal disk with special integrator for close encounters.
 *      Particle id's: 0 = star, 1 = massive body, 2 = planetesimal, 3 = CLOSE ENCOUNTER
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "rebound.h"
#include "tools.h"
#include "../examples/planetesimals2/functions.h"

void heartbeat(struct reb_simulation* r);
double tmax, planetesimal_mass, CE_exit_time = 0;
int n_output, CE_index = 0, encounter_index = 0, previous_encounter_index = 0;
char plntdir[200] = "output/planet_", lgnddir[200] = "output/planet_";
struct reb_simulation* s;

int main(int argc, char* argv[]){
	struct reb_simulation* r = reb_create_simulation();
	// Setup constants
    tmax = 50;
	r->integrator	= REB_INTEGRATOR_IAS15;    //REB_INTEGRATOR_IAS15 = 0, WHFAST = 1, HYBRID = 5
	r->collision	= REB_COLLISION_NONE;
	r->boundary     = REB_BOUNDARY_OPEN;
	r->heartbeat	= heartbeat;
    r->additional_forces = planetesimal_forces;
    //r->usleep   = 20000; //larger the number, slower OpenGL simulation
    
    // System constants
    double N_Rhill = atof(argv[1]);     //# hill radii for boundary between switch. Try 3?
    double dRHill = atof(argv[2]);      //Number of hill radii buffer. Default = 2?
    double N_planetesimals = 1;
    double M_planetesimals = 5e-5; //Total Mass of all planetesimals (default = Earth mass, 3e-6)
	
    // Other constants
    n_output = 50000;
    double boxsize = 5;
    double kicksperorbit = 50.;
	reb_configure_box(r, boxsize, 1, 1, 1);

	// Initial conditions
	struct reb_particle star = {0};
	star.m 		= 1;
	star.r		= 0.01;
    star.id     = 0;
	reb_add(r, star);
    r->N_active = 1;

    //planet 1
    double a1=0.7, m1=5e-5, e1=0.01;
    struct reb_particle p1 = reb_tools_init_orbit2d(r->G,star.m,m1,a1,e1,0,0);
    p1.id   = 1;
    reb_add(r, p1);
    r->N_active++;
    
    //planet 2
    double a2=1, m2=5e-5, e2=0.01;
    struct reb_particle p2 = reb_tools_init_orbit2d(r->G,star.m,m2,a2,e2,0,0.1);
    p1.id   = 1;
    reb_add(r, p2);
    //r->N_active++;
    
    //calc dt
    if(r->integrator == REB_INTEGRATOR_IAS15){
        r->dt = sqrt(4.0*M_PI*pow(a1,3)/(r->G*star.m))/kicksperorbit;
        printf("dt = %f \n",r->dt);
        N_Rhill = -1;
        dRHill = -1;
    } else r->dt = calc_dt(r, m1, star.m, a1, N_Rhill, dRHill);
    r->dt = 0.015;

    //planetesimals
    //double outer = 3, inner = 14, powerlaw = 0.5;  //higher the inner number, closer to the star
    double outer = 7.05, inner = 7.25, powerlaw = 0.5;
    int seed = 12;          //seed was 11
    srand(seed);
    planetesimal_mass = M_planetesimals / N_planetesimals;  //mass of each planetesimal
/*
	while(r->N<N_planetesimals + r->N_active){
		struct reb_particle pt = {0};
		double a	= reb_random_powerlaw(boxsize/outer,boxsize/inner,powerlaw);
		double phi 	= reb_random_uniform(0,2.*M_PI);
		pt.x 		= a*cos(phi);
		pt.y 		= a*sin(phi);
		pt.z 		= a*reb_random_normal(0.0001);
		double vkep = sqrt(r->G*star.m/a);
		pt.vx 		= -vkep * sin(phi);
		pt.vy 		= vkep * cos(phi);
		pt.m 		= 0;
		pt.r 		= .3/sqrt((double)N_planetesimals);
        pt.id       = 2;
		reb_add(r, pt);
	}
    */
    //naming stuff
    legend(plntdir, lgnddir, r, tmax, r->N_active, r->N, planetesimal_mass, M_planetesimals, inner, outer, powerlaw, m1, a1, e1, star.m, N_Rhill, dRHill);
    
    //initialize mini simulation (IAS15)
    s = reb_create_simulation();
    ini_mini(r,s);
    
    // The WH integrator assumes a heliocentric coordinatesystem.
    if(r->integrator != REB_INTEGRATOR_WH) reb_move_to_com(r);
	reb_integrate(r, tmax); //Integrate! INFINITY is a choice
}

void heartbeat(struct reb_simulation* r){
    if(r->integrator == REB_INTEGRATOR_WHFAST){
        encounter_index = check_for_encounter(r);
        if(encounter_index != 0){
            if(previous_encounter_index == 0){ //initialize mini simulation
                fprintf(stderr,"\n\033[1mParticle %d entering\033[0m Hill Sphere at t=%f. \n",encounter_index, r->t);
                update_mini(r,s,encounter_index);
                previous_encounter_index = encounter_index;
                
            } else if(previous_encounter_index == encounter_index){//ensure same particle
                reb_integrate(s, s->t + s->dt);
                update_global(s,r,encounter_index);
                //printf("time comparison outer: r->t=%f, s->t=%f \n",r->t, s->t);
            }
            else{
                fprintf(stderr,"\n\033[1mWarning!\033[0m %d Just entered the Hill Sphere, but %d still inside too. Multiple Particles in Hill Sphere at t=%f. \n",encounter_index, previous_encounter_index, r->t);
            }
        } else if(previous_encounter_index != 0){//particle left hill sphere, update final time
             fprintf(stderr,"\n\033[1mParticle %d leaving\033[0m Hill Sphere at t=%f.\n",previous_encounter_index, r->t);
            reb_integrate(s, s->t + s->dt);
            update_global(s,r,previous_encounter_index);
            previous_encounter_index = 0;
        }
        
    }
    
    //output stuff
    if (reb_output_check(r, tmax/n_output)){
        double a_p=0, e_p=0, Etot=0, Ltot=0;
        calc_ae(&a_p, &e_p, r);
        calc_ELtot(&Etot, &Ltot, planetesimal_mass, r);
        
        FILE *append;
        append=fopen(plntdir, "a");
        fprintf(append,"%f,%.8f,%.8f,%.15f,%.15f\n",r->t,a_p,e_p,Etot,Ltot);
        fclose(append);
        if(r->integrator != REB_INTEGRATOR_WH) reb_move_to_com(r);
        
		reb_output_timing(r, 0);
	}
}
