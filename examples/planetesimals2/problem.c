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
#include <string.h>
#include "rebound.h"
#include "integrator_whfast.h"
#include "tools.h"
#include "../examples/planetesimals2/functions.h"

void heartbeat(struct reb_simulation* r);
char plntdir[200] = "output/planet_", lgnddir[200] = "output/planet_", xyz_check[200]="output/planet_", CEprint[200]="output/planet_";

double tmax, planetesimal_mass, E0, K0, U0, n_output, dt_ini;
int N_encounters = 0, N_encounters_previous, N_encounters_tot = 0, HYBRID_ON, err_print_msg = 0;
int* encounter_index; int* previous_encounter_index; double* Hill2; double* x_prev; double* y_prev; double* z_prev; double t_prev;
struct reb_simulation* s; struct reb_simulation* r;

int main(int argc, char* argv[]){
    //switches
    int turn_planetesimal_forces_on = 1;
    HYBRID_ON = 1;
    
    //System constants
    tmax = atoi(argv[1]);
    int N_planetesimals = atoi(argv[2]);
    double M_planetesimals = 3e-6;                        //Total Mass of all planetesimals (default = Earth mass, 3e-6)
    planetesimal_mass = M_planetesimals/N_planetesimals;  //mass of each planetesimal
    double ias_epsilon = 1e-7;                            //sets precision of ias15
    double HSR2 = atof(argv[3]);                          //Transition boundary between WHFAST and IAS15. Units of Hill^2
    double dRHill = atof(argv[4]);                        //Sets the timestep - max # Hill radii/timestep. Smaller=stricter
    int seed = atoi(argv[5]);
    
	//Simulation Setup
    r = reb_create_simulation();
	r->integrator	= 1;                                  //REB_INTEGRATOR_IAS15 = 0, WHFAST = 1, WH=3, HYBRID = 5
	r->collision	= REB_COLLISION_NONE;
	r->heartbeat	= heartbeat;
    r->ri_hybrid.switch_ratio = HSR2;
    if(turn_planetesimal_forces_on==1)r->additional_forces = planetesimal_forces_global;
    //r->ri_whfast.corrector 	= 11;
    //r->usleep   = 5000; //larger the number, slower OpenGL simulation
	
    // Other setup stuff
    srand(seed);
    n_output = 100000;
    //r->boundary     = REB_BOUNDARY_OPEN;
    //double boxsize = 20;
	//reb_configure_box(r, boxsize, 3, 3, 1);

	// Initial conditions
	struct reb_particle star = {0};
	star.m 		= 1;
	star.r		= 0.01;
    star.id     = 0;        // 0 = star
	reb_add(r, star);

    //planet 1
    double a1=0.7, m1=5e-5, e1=0.01, inc1 = reb_random_normal(0.00001);
    struct reb_particle p1 = {0};
    p1 = reb_tools_orbit_to_particle(r->G, star, m1, a1, e1, inc1, 0, 0, 0);
    p1.r = 0.1;
    p1.id = 1;              //1 = planet
    reb_add(r, p1);
    
    //planet 2
    double a2=1, m2=5e-5, e2=0.01, inc2=reb_random_normal(0.00001);
    struct reb_particle p2 = {0};
    p2 = reb_tools_orbit_to_particle(r->G, star, m2, a2, e2, inc2, 0, 0, 0);
    p2.r = 0.1;
    p2.id = 1;              //1 = planet
    reb_add(r, p2);
    
    //N_active and move to COM
    r->N_active = r->N;
    if(r->integrator != REB_INTEGRATOR_WH) reb_move_to_com(r);
    
    //calc dt
    dt_ini = calc_dt(r, m1, star.m, a1, dRHill);
    if(r->integrator == REB_INTEGRATOR_IAS15){
        //double kicksperorbit = 50.;
        //dt_ini = sqrt(4.0*M_PI*pow(a1,3)/(r->G*star.m))/kicksperorbit;
        dt_ini /= 3.;
        printf("dt = %f \n",dt_ini);
        dRHill = -1;
    }
    r->dt = dt_ini;
    
    //planetesimals
    double planetesimal_buffer = 0.1;   //Chatterjee & Ford use 0.01
    double inner = a1 - planetesimal_buffer, outer = a2 + planetesimal_buffer, powerlaw = 0.5;
    while(r->N<N_planetesimals + r->N_active){
		struct reb_particle pt = {0};
		double a	= reb_random_powerlaw(inner,outer,powerlaw);
        double phi 	= reb_random_uniform(0,2.*M_PI);
        double inc = reb_random_normal(0.0001);
        double Omega = reb_random_uniform(0,2.*M_PI);
        double apsis = reb_random_uniform(0,2.*M_PI);
        pt = reb_tools_orbit_to_particle(r->G, star, planetesimal_mass, a, 0, inc, Omega, apsis,phi);
		pt.r 		= 0.005;
        pt.id       = r->N;
		reb_add(r, pt);
	}
    
    //Ini malloc arrays
    x_prev = calloc(sizeof(double),r->N);           //Previous global positions for interpolating
    y_prev = calloc(sizeof(double),r->N);
    z_prev = calloc(sizeof(double),r->N);
    encounter_index = malloc(sizeof(int));          //encounter index
    previous_encounter_index = malloc(sizeof(int));
    Hill2 = calloc(sizeof(double),r->N);             //Hill radius squared for fast calc.
    calc_Hill2(r);
    
    //Initializing stuff
    legend(plntdir, lgnddir, xyz_check, CEprint, r, tmax, planetesimal_mass, M_planetesimals, N_planetesimals,inner, outer, powerlaw, m1, a1, e1, star.m, dRHill,ias_epsilon,seed,HYBRID_ON);
    s = reb_create_simulation();    //initialize mini simulation (IAS15)
    ini_mini(r,s,ias_epsilon,turn_planetesimal_forces_on);
    E0 = calc_Etot(r, &K0, &U0);
    
    //Integrate!
    clock_t t_ini = clock_start();
    reb_integrate(r, tmax);
    
    //finish
    clock_finish(t_ini,N_encounters_tot,lgnddir);
    global_free();
}

void heartbeat(struct reb_simulation* r){
    double K = 0, U = 0, min_r = 0, max_val = 0;
    double E1 = calc_Etot(r,&K,&U);
    
    if(HYBRID_ON == 1){
        check_for_encounter(r, s, &N_encounters, &min_r, &max_val, dt_ini);
        int dN = N_encounters - N_encounters_previous;
        if(N_encounters_previous == 0){
            if(N_encounters > 0){//1st update in a while, update mini massive bodies, add particles, no int
                s->t = r->t;
                int N_active = s->N_active;
                struct reb_particle* global = r->particles;
                struct reb_particle* mini = s->particles;
                for(int i=0; i<N_active; i++) mini[i] = global[i];
                add_or_subtract_particles(r,s,N_encounters,N_encounters_previous,dN,CEprint);
                update_previous_global_positions(r, N_encounters);
            } //otherwise do nothing.
        } else { //integrate existing mini, update global, add/remove new/old particles.
            reb_integrate(s, r->t);
            update_global(s,r,N_encounters_previous, N_encounters);
            add_or_subtract_particles(r,s,N_encounters,N_encounters_previous,dN,CEprint);
            update_previous_global_positions(r, N_encounters);
        }
        update_encounter_indices(&N_encounters, &N_encounters_previous);
    }
    
    //OUTPUT stuff*******
    if(reb_output_check(r,tmax/n_output)){
        FILE *append;
        append = fopen(plntdir, "a");
        fprintf(append, "%.16f,%.16f, %d, %.12f,%.12f,%.16f,%.16f,%.16f\n",r->t,s->t,N_encounters_previous,min_r,max_val,fabs((E1 - E0)/E0),fabs((K - K0)/K0),fabs((U - U0)/U0));
        fclose(append);
        
        //output error stuff
        if(fabs((E1 - E0)/E0) > 1e-5 && err_print_msg == 0){
            err_print_msg++;
            fprintf(stderr,"\n\033[1mERROR EXCEEDED for %s.\033[0m.\n",plntdir);
        }
        reb_output_timing(r, 0);    //output only when outputting values. Saves some time
    }
    /*
    if(reb_output_check(r,tmax/100)){
        FILE *xyz_output;
        xyz_output = fopen(xyz_check, "a");
        struct reb_particle* global = r->particles;
        fprintf(xyz_output, "%.16f\n",r->t);
        for(int i=0;i<r->N;i++){
            fprintf(xyz_output, "%d,%.16f,%.16f,%.16f,%.16f,%.16f,%.16f,%.16f,%.16f,%.16f\n",i,global[i].x,global[i].y,global[i].z,global[i].vx,global[i].vy,global[i].vz,global[i].ax,global[i].ay,global[i].az);
        }
        fclose(xyz_output);
    }*/
}


//OUTPUT stuff*******
/*double E_curr = 0, K_curr = 0, U_curr = 0, L_curr = 0, a_p = 0, d_p = 0, e_p = 0, t = r->t;
 calc_ELtot(&E_curr, &K_curr, &U_curr, &L_curr, planetesimal_mass, r); //calcs Etot all in one go.
 for(int i=1;i<r->N_active;i++){
 calc_ae(&a_p, &e_p, &d_p, r, i, t);
 
 FILE *append;
 append=fopen(plntdir, "a");
 fprintf(append,"%f,%.8f,%.8f,%.16f,%.16f,%.16f,%.16f,%.16f\n",t,a_p,e_p,fabs((E_ini - E_curr)/E_ini),fabs((K_ini - K_curr)/K_ini), fabs((U_ini - U_curr)/U_ini),fabs((L_ini - L_curr)/L_ini),d_p);
 fclose(append);
 E_curr = E_ini; L_curr = L_ini;
 }*/
