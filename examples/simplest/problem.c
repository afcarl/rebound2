/**
 * A very simple test problem
 * 
 * We first create a REBOUND simulation, then we add 
 * two particles and integrate the system for 100 time 
 * units.
 */
#include "rebound.h"
#include <stdio.h>
#include <stdlib.h>

void heartbeat(struct reb_simulation* r){
	//printf("%f\n",r->t);
}

int main(int argc, char* argv[]) {
	struct reb_simulation* r = reb_create_simulation();
	r->dt = 0.1;
	r->heartbeat = heartbeat;
    //r->usleep = 80000;
    r->integrator=REB_INTEGRATOR_WHFAST;
	r->exact_finish_time = 1; // Finish exactly at tmax in reb_integrate(). Default is already 1.

	struct reb_particle p1 = {0};
	p1.m = 1.;
	reb_add(r, p1);
	
	struct reb_particle p2 = {0};
	p2.x = 1;
	p2.vy = 1;
	reb_add(r, p2);

    struct reb_particle* particles = r->particles;
    struct reb_particle out = particles[1];
    printf("values ini: x=%f,y=%f,vx=%f,vy=%f,t=%f\n",out.x,out.y,out.vx,out.vy,r->t);
    reb_integrate(r,50.4);
    /*
    for(int i=0;i<10;i++){
        reb_integrate(r,r->t+5.);
        particles = r->particles;
        out = particles[1];
        printf("values: x=%f,y=%f,vx=%f,vy=%f,t=%f\n",out.x,out.y,out.vx,out.vy,r->t);
        
    }*/
    particles = r->particles;
    out = particles[1];
    printf("values fini: x=%f,y=%f,vx=%f,vy=%f,t=%f\n",out.x,out.y,out.vx,out.vy,r->t);
}

