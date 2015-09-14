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
#include "integrator_whfast.h"
#include "tools.h"
#include "../examples/planetesimals2/functions.h"

void heartbeat(struct reb_simulation* r);
double tmax, planetesimal_mass, CE_exit_time = 0, E_ini, K_ini, U_ini, L_ini, n_output;
int N_encounters = 0, N_encounters_previous, N_encounters_tot = 0, HYBRID_ON, output_counter = 0, temp = 0;
int* encounter_index = NULL; int* previous_encounter_index = NULL; double* Hill = NULL;
char plntdir[200] = "output/planet_", lgnddir[200] = "output/planet_";
struct reb_simulation* s;

int main(int argc, char* argv[]){
    // System constants
    tmax = 500;
    HYBRID_ON = 1;
    double dRHill = atof(argv[2]);      //Number of hill radii buffer. Sets the timestep. Smaller = stricter
    double N_planetesimals = 20;
    double M_planetesimals = 3e-6; //Total Mass of all planetesimals (default = Earth mass, 3e-6)
    
    struct reb_simulation* r = reb_create_simulation();
	// Setup constants
	r->integrator	= atoi(argv[3]);    //REB_INTEGRATOR_IAS15 = 0, WHFAST = 1, WH=3, HYBRID = 5
	r->collision	= REB_COLLISION_NONE;
	r->boundary     = REB_BOUNDARY_OPEN;
	r->heartbeat	= heartbeat;
    r->additional_forces = planetesimal_forces;
    r->ri_hybrid.switch_ratio = atof(argv[1]);     //# hill radii for boundary between switch. Try 3?
    //r->usleep   = 5000; //larger the number, slower OpenGL simulation
	
    // Other constants
    n_output = 10000;
    double boxsize = 5;
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
    struct reb_particle p1 = reb_tools_init_orbit3d(r->G, star.m, m1, a1, e1, reb_random_normal(0.0001), 0, 0, 0);
    p1.id = 1;              //1 = planet
    reb_add(r, p1);
    r->N_active++;
    
    //planet 2
    double a2=1, m2=5e-5, e2=0.01;
    struct reb_particle p2 = reb_tools_init_orbit3d(r->G, star.m, m2, a2, e2, reb_random_normal(0.0001), 0, 0, 0);
    p2.id = 1;              //1 = planet
    reb_add(r, p2);
    r->N_active++;
    
    //calc dt
    if(r->integrator == REB_INTEGRATOR_IAS15){
        double kicksperorbit = 50.;
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
        double inc = reb_random_normal(0.0001);
        //double Omega = reb_random_uniform(0,2.*M_PI);
        //double apsis = reb_random_uniform(0,2.*M_PI);
        //double a = 0.664171;
        //double phi=0.5;
        double Omega = 0;
        double apsis = 0;
        //double a = 0.695;
        //double phi = 0.03;
        //double a = 0.3;
        //double phi=0.5;
        pt.m 		= planetesimal_mass;
        pt = reb_tools_init_orbit3d(r->G, star.m, planetesimal_mass, a, 0, inc, Omega, apsis,phi);
		pt.r 		= 0.3;
        pt.id       = r->N;
		reb_add(r, pt);
	}
    
    //move to COM - before planetesimals?
    if(r->integrator != REB_INTEGRATOR_WH) reb_move_to_com(r);
    
    //Hill Sphere (for speed in check_for_encounter)
    Hill = calloc(sizeof(double),r->N);
    struct reb_particle* restrict const particles = r->particles;
    struct reb_particle p0 = particles[0];
    for(int i=1;i<r->N;i++){
        struct reb_particle body = particles[i];
        double mp;
        if(i>=r->N_active) mp = planetesimal_mass; else mp = body.m;
        Hill[i] = pow((mp/(p0.m*3.)), 2./3.);
    }
    
    //Initializing stuff
    legend(plntdir, lgnddir, r, tmax, planetesimal_mass, M_planetesimals, inner, outer, powerlaw, m1, a1, e1, star.m, dRHill,HYBRID_ON);
    s = reb_create_simulation();    //initialize mini simulation (IAS15)
    ini_mini(r,s);
    calc_ELtot(&E_ini, &K_ini, &U_ini, &L_ini, planetesimal_mass, r);
    clock_t timer = clock();
    encounter_index = malloc(sizeof(int));
    previous_encounter_index = malloc(sizeof(int));
    
    //Integrate!
    reb_integrate(r, tmax);
    
    //finish
    clock_finish(timer,N_encounters_tot,lgnddir);
    free(encounter_index);
    free(previous_encounter_index);
}

void heartbeat(struct reb_simulation* r){
    if(HYBRID_ON == 1){
        check_for_encounter(r, &N_encounters);
        int dN = N_encounters - N_encounters_previous;
        
        if(N_encounters_previous == 0){
            if(N_encounters > 0){
                //first update in a while, only update massive bodies in mini and add any particles
                s->t = r->t;
                int N_active = s->N_active;
                struct reb_particle* global = r->particles;
                struct reb_particle* mini = s->particles;
                for(int i=0; i<N_active; i++) mini[i] = global[i];
                add_or_subtract_particles(r,s,N_encounters,N_encounters_previous,dN);
            } //otherwise do nothing.
        } else {
            //integrate existing mini, update global, add/remove new/old particles.
            reb_integrate(s, r->t);
            update_global(s,r,N_encounters_previous, N_encounters);
            add_or_subtract_particles(r,s,N_encounters,N_encounters_previous,dN);
        }
        update_encounter_indices(&N_encounters, &N_encounters_previous);
    }
    
    //output stuff
    if(r->t > output_counter*tmax/n_output){
        output_counter++;
        double E_curr = 0, K_curr = 0, U_curr = 0, L_curr = 0, a_p = 0, d_p = 0, e_p = 0, t = r->t;
        calc_ELtot(&E_curr, &K_curr, &U_curr, &L_curr, planetesimal_mass, r); //calcs Etot all in one go.
        for(int i=1;i<r->N_active;i++){
            calc_ae(&a_p, &e_p, &d_p, r, i, t);
            
            FILE *append;
            append=fopen(plntdir, "a");
            fprintf(append,"%f,%.8f,%.8f,%.16f,%.16f,%.16f,%.16f,%.16f\n",t,a_p,e_p,fabs((E_ini - E_curr)/E_ini),fabs((K_ini - K_curr)/K_ini), fabs((U_ini - U_curr)/U_ini),fabs((L_ini - L_curr)/L_ini),d_p);
            fclose(append);
            E_curr = E_ini; L_curr = L_ini;
        }
        reb_output_timing(r, 0);
        
        /*
        FILE* output;
        output=fopen("debug/HYB_sept14_par11_massless.txt","a");
        struct reb_particle* global = r->particles;
        for(int i=0;i<r->N_active;i++) fprintf(output,"%f,%.15f,%.15f,%.15f,%.15f,%.15f,%.15f\n",r->t,global[i].x,global[i].y,global[i].z,global[i].vx,global[i].vy,global[i].vz);
        fprintf(output,"%f,%.15f,%.15f,%.15f,%.15f,%.15f,%.15f\n",r->t,global[11].x,global[11].y,global[11].z,global[11].vx,global[11].vy,global[11].vz);
        fclose(output);
        */
         
        /*
        E_curr = 0, K_curr = 0, U_curr = 0, L_curr = 0, a_p = 0, d_p = 0, e_p = 0;
        calc_ELtot(&E_curr, &K_curr, &U_curr, &L_curr, planetesimal_mass, s); //calcs Etot all in one go.
        for(int i=1;i<s->N_active;i++){
            calc_ae(&a_p, &e_p, &d_p, s, i);
            
            FILE *append;
            append=fopen("output/mini_output.txt", "a");
            fprintf(append,"%f,%.8f,%.8f,%.16f,%.16f,%.16f,%.16f,%.8f\n",t,a_p,e_p,fabs((E_ini - E_curr)/E_ini),fabs((K_ini - K_curr)/K_ini), fabs((U_ini - U_curr)/U_ini),fabs((L_ini - L_curr)/L_ini),d_p);
            fclose(append);
            E_curr = E_ini; L_curr = L_ini;
        }*/
    }
}
