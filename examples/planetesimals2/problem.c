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

double dxold1=0,dyold1=0,dzold1=0,dxold2=0,dyold2=0,dzold2=0;
double tmax, planetesimal_mass, E0, K0, U0, n_output, dt_ini, t_output, t_log_output;
int N_encounters = 0, N_encounters_previous, N_encounters_tot = 0, HYBRID_ON, err_print_msg = 0, n_o=0;
int* encounter_index; int* previous_encounter_index; double* Hill2; double* x_prev; double* y_prev; double* z_prev; double t_prev;
struct reb_simulation* s; struct reb_simulation* r;

int main(int argc, char* argv[]){
    //switches
    int turn_planetesimal_forces_on = 1;
    int p1_satellite_on = 0;
    HYBRID_ON = 1;
    
    //System constants
    tmax = atof(argv[1]);
    int N_planetesimals = atoi(argv[2]);
    double M_planetesimals = 3e-6;                        //Tot. Mass of all planetesimals (Earth mass, 3e-6)
    planetesimal_mass = M_planetesimals/N_planetesimals;  //mass of each planetesimal
    double ias_epsilon = 1e-9;                            //sets precision of ias15
    double HSR2 = atof(argv[3]);                          //Transition boundary bet. WHFAST & IAS15. Units of Hill^2
    double dRHill = atof(argv[4]);                        //Sets the timestep - max # Hill radii/timestep.
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

    //Boundary stuff
    //r->boundary     = REB_BOUNDARY_OPEN;
    //double boxsize = 20;
    //reb_configure_box(r, boxsize, 3, 3, 1);
    
    srand(seed);        //This needs to be here, just before rand() is called
    
	// Initial conditions
	struct reb_particle star = {0};
	star.m 		= 1;
	star.r		= 0.005;        //I think radius of particle is in AU!
    star.id     = 0;            // 0 = star
	reb_add(r, star);

    //planet 1
    double a1=0.7, m1=5e-4, e1=0, inc1 = reb_random_normal(0.00001);
    struct reb_particle p1 = {0};
    p1 = reb_tools_orbit_to_particle(r->G, star, m1, a1, e1, inc1, 0, 0, 0);
    p1.r = 1.6e-4;              //I think radius of particle is in AU!
    p1.id = r->N;
    reb_add(r, p1);
    
    //planet 2
    double a2=1, m2=5e-5, e2=0.01, inc2=reb_random_normal(0.00001);
    struct reb_particle p2 = {0};
    p2 = reb_tools_orbit_to_particle(r->G, star, m2, a2, e2, inc2, 0, 0, 0);
    p2.r = 0.1;
    p2.id = r->N;
    reb_add(r, p2);
    
    //N_active and move to COM
    r->N_active = r->N;
    if(r->integrator != REB_INTEGRATOR_WH) reb_move_to_com(r);
    
    //calc dt
    dt_ini = calc_dt(r, m1, star.m, a1, dRHill);
    if(r->integrator == REB_INTEGRATOR_IAS15){
        dt_ini /= 3.;
        printf("dt = %f \n",dt_ini);
        dRHill = -1;
    }
    r->dt = dt_ini;
    
    // Other setup stuff
    n_output = 5000;    //true output ~2x n_output for some reason.
    t_log_output = pow(tmax + 1, 1./(n_output - 1));
    t_output = dt_ini;
    
    //orbiting planetesimal/satellite
    if(p1_satellite_on == 1){
        double x=0.005;
        struct reb_particle pt = {0};
        pt = reb_tools_orbit_to_particle(r->G, p1, planetesimal_mass, x, 0, 0, 0, 0, 0);
        pt.y += p1.y;
        pt.r = 4e-5;            //I think radius of particle is in AU!
        pt.id = r->N;              //1 = planet
        reb_add(r, pt);
    }
    
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
		pt.r 		= 4e-5;
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
    printf("N_outputs total=%d\n",n_o);
    global_free();
}

void heartbeat(struct reb_simulation* r){
    double K = 0, U = 0, min_r = 0, max_val = 0;
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
    
    double E1 = calc_Etot(r,&K,&U);
    //output error stuff - every iteration
    if(fabs((E1 - E0)/E0) > 1e-7 && err_print_msg == 0){
        err_print_msg++;
        fprintf(stderr,"\n\033[1mERROR EXCEEDED for %s\033[0m, t=%f.\n",plntdir,r->t);
    }
    
    if(r->t > 115){
        FILE *xyz_output;
        xyz_output = fopen(xyz_check, "a");
        struct reb_particle* global = r->particles;
        //fprintf(xyz_output, "%.16f,rmin=%.16f,vmax/rmin=%.16f\n",r->t,min_r,max_val);
        int i = 42;
        int j=1;
        //fprintf(xyz_output, "particle %d,x=%.16f,y=%.16f,z=%.16f,%.16f,%.16f,%.16f,%.16f,%.16f,%.16f\n",i,global[i].x,global[i].y,global[i].z,global[i].vx,global[i].vy,global[i].vz,global[i].ax,global[i].ay,global[i].az);
        fprintf(xyz_output, "%.16f,%.16f,%.16f,%.16f\n",r->t,fabs(global[i].x-dxold1),fabs(global[i].y-dyold1),fabs(global[i].z-dzold1));
        //fprintf(xyz_output, "planet %d,x=%.16f,y=%.16f,z=%.16f,%.16f,%.16f,%.16f,%.16f,%.16f,%.16f\n",j,global[j].x,global[j].y,global[j].z,global[j].vx,global[j].vy,global[j].vz,global[j].ax,global[j].ay,global[j].az);
        fprintf(xyz_output, "%.16f,%.16f,%.16f,%.16f\n",r->t,fabs(global[j].x-dxold2),fabs(global[j].y-dyold2),fabs(global[j].z-dzold2));
        dxold1 = global[i].x; dyold1 = global[i].y; dzold1 = global[i].z;
        dxold2 = global[j].x; dyold2 = global[j].y; dzold2 = global[j].z;
        //for(int i=0;i<r->N;i++){
        //    fprintf(xyz_output, "%d,%.16f,%.16f,%.16f,%.16f,%.16f,%.16f,%.16f,%.16f,%.16f\n",i,global[i].x,global[i].y,global[i].z,global[i].vx,global[i].vy,global[i].vz,global[i].ax,global[i].ay,global[i].az);
        //}
        fclose(xyz_output);
    }
    
    //OUTPUT stuff*******
    if(r->t > t_output || r->t <= r->dt){
        FILE *append;
        append = fopen(plntdir, "a");
        fprintf(append, "%.16f,%.16f, %d, %.12f,%.12f,%.16f,%.16f,%.16f,%d,%d\n",r->t,s->t,N_encounters_previous,min_r,max_val,fabs((E1 - E0)/E0),fabs((K - K0)/K0),fabs((U - U0)/U0), r->N,s->N);
        fclose(append);
        
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
        
        reb_output_timing(r, 0);    //output only when outputting values. Saves some time
        t_output *= t_log_output;
        n_o++;
    }
    
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

