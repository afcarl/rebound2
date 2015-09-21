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

void heartbeat(struct reb_simulation* r);

struct reb_simulation* s;

double e0;

int main(int argc, char* argv[]){
    struct reb_simulation* r = reb_create_simulation();
    // Setup constants
    r->integrator	= REB_INTEGRATOR_WHFAST;
    r->heartbeat	= heartbeat;
    
    //star
    struct reb_particle star = {0};
    star.m 		= 1;
    star.r		= 0.01;
    reb_add(r, star);
    
    //planet 1
    double a1=0.7, m1=5e-5, e1=0.01;
    struct reb_particle p1 = reb_tools_orbit2d_to_particle(r->G, star, m1, a1, e1, 0., 0.);
    reb_add(r, p1);
    
    //planet 2
    double a2=1, m2=5e-5, e2=0.01;
    struct reb_particle p2 = reb_tools_orbit2d_to_particle(r->G, star, m2, a2, e2, 0., 0.);
    reb_add(r, p2);
    
    
    //move to COM
    reb_move_to_com(r);
    
    //Ini mini
    s = reb_create_simulation();    //initialize mini simulation (IAS15)
    s->integrator = REB_INTEGRATOR_IAS15;
    
    s->exact_finish_time = 1;
    s->dt = r->dt;
    struct reb_particle* restrict const particles = r->particles;
    for(int k=0; k<r->N; k++){
        struct reb_particle p = particles[k];
        reb_add(s,p);
    }
    
    system("rm -vf energy.txt");
    
    //Integrate!
    e0 = reb_tools_energy(r);
    reb_integrate(r, INFINITY);
    
    //********setup calculations and integrate!*********************
}

int mini_on = 0;

void heartbeat(struct reb_simulation* r){
    //Simple algorithm which starts using mini and updating global after 50 years
    if(r->t > 50 && mini_on == 0){
        s->t = r->t;
        struct reb_particle* global = r->particles;
        struct reb_particle* mini = s->particles;
        printf("%d %d \n\n\n",r->N, s->N);
        for(int i=0; i<r->N; i++) mini[i] = global[i];
        
        fprintf(stderr,"\n\033[1m Note!\033[0m Initializing mini simulation (IAS)\n");
        mini_on = 1;
        
        
    } else if(mini_on == 1){
        reb_integrate(s, r->t);
        struct reb_particle* global = r->particles;
        struct reb_particle* mini = s->particles;
        
        //update massive bodies and planetesimal particle
        for(int i=0; i<s->N; i++) global[i] = mini[i];
    }
    
    FILE* f = fopen("energy.txt","a+");
    double e1 = reb_tools_energy(r);
    fprintf(f,"%e %e\n",r->t,fabs((e0-e1)/e0));
    fclose(f); 
    //output stuff
}