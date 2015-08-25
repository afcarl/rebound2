/**
 * A.S. This is my planetesimal disk with special integrator for close encounters.
 *      Particle id's: 0 = star, 1 = massive body, 2 = planetesimal, 3 = CLOSE ENCOUNTER
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include "rebound.h"
#include "tools.h"
#include "../examples/planetesimals2/functions.h"

void heartbeat(struct reb_simulation* r);
double tmax, planetesimal_mass, CE_exit_time = 0, E_ini, L_ini;
int n_output, CE_index = 0, warning = 0, N_encounters = 0, N_encounters_previous;
int* encounter_index = NULL;
int* previous_encounter_index = NULL;
char plntdir[200] = "output/planet_", lgnddir[200] = "output/planet_";
struct reb_simulation* s;

int main(int argc, char* argv[]){
	struct reb_simulation* r = reb_create_simulation();
	// Setup constants
    tmax = 500;
	r->integrator	= atoi(argv[3]);    //REB_INTEGRATOR_IAS15 = 0, WHFAST = 1, HYBRID = 5
	r->collision	= REB_COLLISION_NONE;
	r->boundary     = REB_BOUNDARY_OPEN;
	r->heartbeat	= heartbeat;
    r->additional_forces = planetesimal_forces;
    r->ri_hybrid.switch_ratio = atof(argv[1]);     //# hill radii for boundary between switch. Try 3?
    //r->usleep   = 20000; //larger the number, slower OpenGL simulation
    
    // System constants
    double dRHill = atof(argv[2]);      //Number of hill radii buffer. Sets the timestep. Smaller = stricter
    double N_planetesimals = 50;
    double M_planetesimals = 3e-6; //Total Mass of all planetesimals (default = Earth mass, 3e-6)
	
    // Other constants
    n_output = 50000;
    double boxsize = 5;
    double kicksperorbit = 50.;
	reb_configure_box(r, boxsize, 1, 1, 1);

	// Initial conditions
	struct reb_particle star = {0};
	star.m 		= 1;
	star.r		= 0.01;
    star.id     = 0;        // 0 = star
	reb_add(r, star);
    r->N_active = 1;

    //planet 1
    double a1=0.7, m1=5e-5, e1=0.01;
    struct reb_particle p1 = reb_tools_init_orbit2d(r->G,star.m,m1,a1,e1,0,0);
    p1.id = 1;              //1 = planet
    reb_add(r, p1);
    r->N_active++;
    
    //planet 2
    double a2=1, m2=5e-5, e2=0.01;
    struct reb_particle p2 = reb_tools_init_orbit2d(r->G,star.m,m2,a2,e2,0,0);
    p2.id = 1;              //1 = planet
    reb_add(r, p2);
    r->N_active++;
    
    //calc dt
    if(r->integrator == REB_INTEGRATOR_IAS15){
        r->dt = sqrt(4.0*M_PI*pow(a1,3)/(r->G*star.m))/kicksperorbit;
        printf("dt = %f \n",r->dt);
        dRHill = -1;
    } else r->dt = calc_dt(r, m1, star.m, a1, dRHill);

    //planetesimals
    double outer = 3, inner = 15, powerlaw = 0.5;  //higher the inner number, closer to the star
    int seed = atoi(argv[4]);          //seed was 11
    srand(seed);
    planetesimal_mass = M_planetesimals / N_planetesimals;  //mass of each planetesimal
    while(r->N<N_planetesimals + r->N_active){
		struct reb_particle pt = {0};
		double a	= reb_random_powerlaw(boxsize/outer,boxsize/inner,powerlaw);
		double phi 	= reb_random_uniform(0,2.*M_PI);
        //double a = 0.664171;
        //double phi=0.948809;
		pt.x 		= a*cos(phi);
		pt.y 		= a*sin(phi);
		pt.z 		= a*reb_random_normal(0.0001);
		double vkep = sqrt(r->G*star.m/a);
		pt.vx 		= -vkep * sin(phi);
		pt.vy 		= vkep * cos(phi);
		pt.m 		= 0;
		pt.r 		= .3/sqrt((double)N_planetesimals);
        pt.id       = 2;    //2 = planetesimal
		reb_add(r, pt);
	}
    
    //Initializing stuff
    legend(plntdir, lgnddir, r, tmax, planetesimal_mass, M_planetesimals, inner, outer, powerlaw, m1, a1, e1, star.m, dRHill);
    s = reb_create_simulation();    //initialize mini simulation (IAS15)
    ini_mini(r,s);
    if(r->integrator != REB_INTEGRATOR_WH) reb_move_to_com(r);
    calc_ELtot(&E_ini, &L_ini, planetesimal_mass, r);
    clock_t timer = clock();
    encounter_index = malloc(sizeof(int));
    previous_encounter_index = malloc(sizeof(int));
    
    //Integrate!
	reb_integrate(r, tmax);
    
    //finish
    clock_finish(timer,N_encounters,lgnddir);
    free(encounter_index);
    free(previous_encounter_index);
}

void heartbeat(struct reb_simulation* r){
    if(r->integrator == REB_INTEGRATOR_WHFAST){
        check_for_encounter(r, &encounter_index, &N_encounters);
        int dN = N_encounters - N_encounters_previous;
        if(abs(dN) > 1){
            fprintf(stderr,"\n\033[1mAlert!\033[0m %d Particles entering/leaving Hill sphere simultaneously. Exiting.\n",dN);
            exit(0);
        }
        switch(dN){     //assume that particles can only enter/leave Hill sphere 1 at a time.
            case 0:{    //no new particles entering Hill Sphere
                if(N_encounters != 0){          //particle(s) inside Hill sphere
                    reb_integrate(s, r->t);
                    update_global(s,r,encounter_index,N_encounters);
                }       //else do nothing I think....
            }break;
            case 1:{    //new particle entering, update mini
                if(N_encounters > 1)reb_integrate(s, r->t); //Just update for 1st particle in Hill
                update_and_add_mini(r,s,encounter_index,N_encounters);
            } break;
            case -1:{    //dN < 0, old particle leaving
                reb_integrate(s, r->t);
                int removal_id = 3;
                compare_encounter_indices(s,encounter_index,previous_encounter_index,N_encounters,removal_id);
                update_and_subtract_mini(r,s,previous_encounter_index,N_encounters,removal_id);
            } break;
        }
        update_encounter_indices(&encounter_index, &previous_encounter_index, &N_encounters, &N_encounters_previous);

        /*
         struct reb_particle* global = r->particles;
         struct reb_particle* mini = s->particles;
         for(int i=0;i<r->N_active+N_encounters;i++){
         printf("mini %d, x,y,vx,vy=%f,%f,%f,%f\n",i,mini[i].x,mini[i].y,mini[i].vx,mini[i].vy);
         printf("mini %d, x,y,vx,vy=%f,%f,%f,%f\n",i,global[i].x,global[i].y,global[i].vx,global[i].vy);
         }
         int i = encounter_index[0];
         printf("mini %d, x,y,vx,vy=%f,%f,%f,%f\n",i,global[i].x,global[i].y,global[i].vx,global[i].vy);
         exit(0);*/
        /*
            if(previous_encounter_index == 0){ //initialize mini simulation
                //fprintf(stderr,"\n\033[1mParticle %d entering\033[0m Hill Sphere at t=%f. \n",encounter_index, r->t);
                update_mini(r,s,encounter_index);
                previous_encounter_index = encounter_index;
                
            } else if(previous_encounter_index == encounter_index){//ensure same particle
                reb_integrate(s, r->t);
                update_global(s,r,encounter_index);
                //printf("time comparison: r->t=%f, s->t=%f, s->dt=%f \n",r->t, s->t, s->dt);
            }
            else{
                if(warning == 0){
                    fprintf(stderr,"\n\033[1mWarning!\033[0m %d Just entered the Hill Sphere, but %d still inside too. Multiple Particles in Hill Sphere at t=%f. \n",encounter_index, previous_encounter_index, r->t);
                    warning++;
                }
            }
        } else if(previous_encounter_index != 0){//particle left hill sphere, update final time
            fprintf(stderr,"\n\033[1mParticle %d leaving\033[0m Hill Sphere at t=%f.\n",previous_encounter_index, r->t);
            reb_integrate(s, r->t);
            update_global(s,r,previous_encounter_index);
            previous_encounter_index = 0;
            N_encounters++;
        }*/
    }
    
    //output stuff
    if (reb_output_check(r, tmax/n_output)){
        double E_curr = 0, L_curr = 0, a_p = 0, d_p = 0, e_p = 0, t = r->t;
        calc_ELtot(&E_curr, &L_curr, planetesimal_mass, r); //calcs Etot all in one go.
        for(int i=1;i<r->N_active;i++){
            calc_ae(&a_p, &e_p, &d_p, r, i);
            
            FILE *append;
            append=fopen(plntdir, "a");
            fprintf(append,"%f,%.8f,%.8f,%e,%.16f,%.8f\n",t,a_p,e_p,fabs((E_ini - E_curr)/E_ini),fabs((L_ini - L_curr)/L_ini),d_p);
            fclose(append);
            E_curr = E_ini; L_curr = L_ini;
        }
		reb_output_timing(r, 0);
	}
}
