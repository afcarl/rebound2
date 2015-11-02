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

double tmax, planetesimal_mass, E0, n_output, dt_ini, t_output, t_log_output, ias_timestep, soft;
int N_encounters = 0, N_encounters_previous, N_encounters_tot = 0, HYBRID_ON, err_print_msg = 0, n_o=0, xyz_output_counter=0, output_error = 0;
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
    //double M_planetesimals = 3e-6;                        //Tot. Mass of all planetesimals (Earth mass, 3e-6)
    //planetesimal_mass = M_planetesimals/N_planetesimals;  //mass of each planetesimal
    planetesimal_mass = 3e-8;                               //each is a moon
    double M_planetesimals = planetesimal_mass*N_planetesimals;
    double ias_epsilon = 1e-8;                              //sets precision of ias15
    double HSR2 = 5;                                        //Transition boundary bet. WHFAST & IAS15. Units of Hill^2
    double dRHill = 0.125;                                   //Sets the timestep - max # Hill radii/timestep.
    soft = 1.6e-4/100.;                                     //gravity softening length scale in AU. R_Neptune/100.
    int seed = atoi(argv[3]);
    
	//Simulation Setup
    r = reb_create_simulation();
	r->integrator	= 1;                                  //REB_INTEGRATOR_IAS15 = 0, WHFAST = 1, WH=3, HYBRID = 5
	r->collision	= REB_COLLISION_NONE;
	r->heartbeat	= heartbeat;
    r->ri_hybrid.switch_ratio = HSR2;
    r->softening = soft;
    if(turn_planetesimal_forces_on==1)r->additional_forces = planetesimal_forces_global;
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
    dt_ini = calc_dt(r, m1, star.m, a1, dRHill, 1);
    
    //planet 2
    double a2=1, m2=5e-5, e2=0.01, inc2=reb_random_normal(0.00001);
    struct reb_particle p2 = {0};
    p2 = reb_tools_orbit_to_particle(r->G, star, m2, a2, e2, inc2, 0, 0, 0);
    p2.r = 0.1;
    p2.id = r->N;
    reb_add(r, p2);
    dt_ini = calc_dt(r, m2, star.m, a2, dRHill, dt_ini);
    
    //calc dt
    if(r->integrator == REB_INTEGRATOR_IAS15){
        dt_ini /= 3.;
        printf("dt = %f \n",dt_ini);
        dRHill = -1;
    }
    printf("timesetep is dt = %f, ri_hybrid.switch_ratio=%f \n",dt_ini,r->ri_hybrid.switch_ratio);
    r->dt = dt_ini;
    
    //N_active and move to COM
    r->N_active = r->N;
    if(r->integrator != REB_INTEGRATOR_WH) reb_move_to_com(r);
    
    //Outputting points
    n_output = 50000;
    t_log_output = pow(tmax + 1, 1./(n_output - 1));
    t_output = dt_ini;
    
    //orbiting planetesimal/satellite
    if(p1_satellite_on == 1){
        double x=0.01;
        struct reb_particle pt = {0};
        //pt = reb_tools_orbit_to_particle(r->G, p1, planetesimal_mass, x, 0, 0, 0, 0, 0);
        //pt.y += p1.y;
        pt = reb_tools_orbit_to_particle(r->G, star, 0, x + a1, 0, 0, 0, 0, 0); //m=planetesimal_mass?
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
        pt = reb_tools_orbit_to_particle(r->G, star, 0, a, 0, inc, Omega, apsis,phi);
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
    E0 = calc_Etot(r, soft);
    
    //Ini mini
    s = reb_create_simulation();    //initialize mini simulation (IAS15)
    double ias_subtime = 3.;        //how much smaller is the ias timestep vs. global?
    ias_timestep = r->dt/ias_subtime;
    ini_mini(r,s,ias_epsilon,turn_planetesimal_forces_on,ias_timestep,soft);
    clock_t t_ini = clock_start();
    
    //Integrate!
    reb_integrate(r, tmax);
    
    //finish
    clock_finish(t_ini,N_encounters_tot,lgnddir);
    printf("N_outputs total=%d\n",n_o);
    global_free();
}

void heartbeat(struct reb_simulation* r){
    double min_r = 1e8, max_val = 1e-8;
    int output_error = 0;
    if(HYBRID_ON == 1){
        if(N_encounters_previous == 0){
            check_for_encounter(r, &N_encounters, N_encounters_previous, &min_r, &max_val, xyz_check);
            if(N_encounters > 0){//1st update in a while, update mini massive bodies, add particles, no int
                s->t = r->t;
                int N_active = s->N_active;
                struct reb_particle* global = r->particles;
                struct reb_particle* mini = s->particles;
                for(int i=0; i<N_active; i++) mini[i] = global[i];
                add_or_subtract_particles(r,s,N_encounters,N_encounters_previous,CEprint);
                update_previous_global_positions(r, N_encounters);
            } //otherwise do nothing.
        } else { //integrate existing mini, update global, add/remove new/old particles.
            reb_integrate(s, r->t);
            update_global(s,r,N_encounters_previous);
            check_for_encounter(r, &N_encounters, N_encounters_previous, &min_r, &max_val, xyz_check);
            add_or_subtract_particles(r,s,N_encounters,N_encounters_previous,CEprint);
            update_previous_global_positions(r, N_encounters);
        }
        update_encounter_indices(&N_encounters, &N_encounters_previous);
    }
    
    double E1 = calc_Etot(r, soft);

    //output error stuff - every iteration
    if(fabs((E1 - E0)/E0) > 1e-6){
        if(err_print_msg == 0){
            err_print_msg++;
            fprintf(stderr,"\n\033[1mERROR EXCEEDED for %s\033[0m, t=%f.\n",plntdir,r->t);
        }
        output_error = 1;
    }
    
    //OUTPUT stuff*******
    if(r->t > t_output || r->t <= r->dt || output_error == 1){
        FILE *append;
        append = fopen(plntdir, "a");
        fprintf(append, "%.16f,%.16f, %d, %.12f,%.12f,%.16f\n",r->t,s->t,r->N,min_r,max_val,fabs((E1 - E0)/E0));
        fclose(append);
        
        reb_output_timing(r, 0);    //output only when outputting values. Saves some time
        t_output = r->t*t_log_output;
        n_o++;
    }
    
    /*
    if(r->t > 112.6 && r->t < 113.2){
        struct reb_particle* global = r->particles;
        struct reb_particle* mini = s->particles;
        char strxyz[50] = {0};
        char temp[7];
        strcat(strxyz,"xyz_temp/outOct28mini_");
        sprintf(temp, "%d",xyz_output_counter);
        strcat(strxyz,temp);
        strcat(strxyz,".txt");
        FILE *xyz_output;
        xyz_output = fopen(strxyz,"w");
        //for(int i=0;i<r->N;i++) fprintf(xyz_output, "%f,%d,%.16f,%.16f,%.16f\n",r->t,global[i].id,global[i].x,global[i].y,global[i].z);
        for(int i=0;i<s->N;i++) fprintf(xyz_output, "%f,%d,%.16f,%.16f,%.16f\n",s->t,mini[i].id,mini[i].x,mini[i].y,mini[i].z);
        fclose(xyz_output);
        xyz_output_counter++;
    }
    */
    /*
    //double t1 = 48.15; double t2 = 48.4;
    double t1 = 201.79; double t2 = 202;
    if(r->t > t1 && r->t < t2){
        int id = 236;
        //int id = 94;
        int pid = 1;
        struct reb_particle* mini = s->particles;
        for(int i=0;i<s->N;i++){
            if(mini[i].id == id){
                double dx = mini[i].x - mini[pid].x;
                double dy = mini[i].y - mini[pid].y;
                double dz = mini[i].z - mini[pid].z;
                double d = sqrt(dx*dx + dy*dy + dz*dz);
                double vx = mini[i].vx - mini[pid].vx;
                double vy = mini[i].vy - mini[pid].vy;
                double vz = mini[i].vz - mini[pid].vz;
                double v = sqrt(vx*vx + vy*vy + vz*vz);
                printf("r->t=%f,s->dt=%.8f,v=%f,d=%f,vx=%.16f,vy=%.16f,vz=%.16f,x=%.16f,y=%.16f,z=%.16f\n",r->t,s->dt,v,d,vx,vy,vz,dx,dy,dz);
            }
        }
    }*/
}

/*
 for(int i=0;i<r->N;i++){
 if(fabs(global[i].ax) + fabs(global[i].ay) + fabs(global[i].az) > 100){
 for(int j=0;j<N_encounters_previous;j++){
 if(global[i].id == previous_encounter_index[j]){
 printf("mini integrated large acc, par %d: ax=%f,ay=%f,az=%f\n",global[i].id,global[i].ax,global[i].ay,global[i].az);
 }
 }
 }
 }*/

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

/*
 if(r->t > 74.5){
 FILE *xyz_output;
 xyz_output = fopen(xyz_check, "a");
 struct reb_particle* global = r->particles;
 //fprintf(xyz_output, "%.16f,rmin=%.16f,vmax/rmin=%.16f\n",r->t,min_r,max_val);
 int i = 20;
 int j=1;
 //fprintf(xyz_output, "particle %d,x=%.16f,y=%.16f,z=%.16f,%.16f,%.16f,%.16f,%.16f,%.16f,%.16f\n",i,global[i].x,global[i].y,global[i].z,global[i].vx,global[i].vy,global[i].vz,global[i].ax,global[i].ay,global[i].az);
 //fprintf(xyz_output, "%.16f,%.16f,%.16f,%.16f,%.16f\n",r->t,s->t,fabs(global[i].x-dxold1),fabs(global[i].y-dyold1),fabs(global[i].z-dzold1));
 //fprintf(xyz_output, "planet %d,x=%.16f,y=%.16f,z=%.16f,%.16f,%.16f,%.16f,%.16f,%.16f,%.16f\n",j,global[j].x,global[j].y,global[j].z,global[j].vx,global[j].vy,global[j].vz,global[j].ax,global[j].ay,global[j].az);
 fprintf(xyz_output, "%.16f,%.16f,%.16f,%.16f,%.16f\n",r->t,s->t,fabs(global[j].x-dxold2),fabs(global[j].y-dyold2),fabs(global[j].z-dzold2));
 dxold1 = global[i].x; dyold1 = global[i].y; dzold1 = global[i].z;
 dxold2 = global[j].x; dyold2 = global[j].y; dzold2 = global[j].z;
 //for(int i=0;i<r->N;i++){
 //    fprintf(xyz_output, "%d,%.16f,%.16f,%.16f,%.16f,%.16f,%.16f,%.16f,%.16f,%.16f\n",i,global[i].x,global[i].y,global[i].z,global[i].vx,global[i].vy,global[i].vz,global[i].ax,global[i].ay,global[i].az);
 //}
 fclose(xyz_output);
 }*/


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

