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
void planetesimal_forces(struct reb_simulation *a);
double calc_Etot(double planetesimal_mass, struct reb_simulation* a);

struct reb_simulation* s;

double e0, planetesimal_mass;
int planetesimal_1, p1_id;

int main(int argc, char* argv[]){
    //********switches******************
    planetesimal_1 = 1;     //=1 to turn on, 0=off
    int planetesimal_2 = 1;
    int turn_planetesimal_forces_on = 1; //=1 to turn on
    planetesimal_mass = 3e-7;
    
    struct reb_simulation* r = reb_create_simulation();
    // Setup constants
    r->integrator	= REB_INTEGRATOR_WHFAST;
    r->heartbeat	= heartbeat;
    r->dt           = 0.005;
    if(turn_planetesimal_forces_on==1)r->additional_forces = planetesimal_forces;
    
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
    //r->N_active = r->N;
    
    //move to COM
    reb_move_to_com(r);
    
    //Ini mini
    s = reb_create_simulation();    //initialize mini simulation (IAS15)
    s->integrator = REB_INTEGRATOR_IAS15;
    s->exact_finish_time = 1;
    s->dt = r->dt;
    if(turn_planetesimal_forces_on==1){
        s->additional_forces = planetesimal_forces;
        fprintf(stderr,"\033[1m Note!\033[0m Planetesimal_forces are on!\n");
    }
    struct reb_particle* restrict const particles = r->particles;
    for(int k=0; k<r->N; k++){
        struct reb_particle p = particles[k];
        reb_add(s,p);
    }
    
    //planetesimal
    if(planetesimal_1 == 1){
        double a=0.5;
        double phi = 5.480386;
        double inc = 0.012694;
        struct reb_particle pt = {0};
        pt = reb_tools_orbit_to_particle(r->G, star, 0, a, 0, inc, 0, 0,phi);
        pt.r 		= 0.04;
        p1_id       = 3;
        pt.id       = p1_id;
        reb_add(r, pt);
        fprintf(stderr, "\033[1m Note!\033[0m planetesimal_1 is on!\n");
        r->N_planetesimals++;
    }
    
    if(planetesimal_2 == 1){
        double a2 = 0.4;
        //a = 2.5;
        double phi2 = 5.380386;
        double inc2 = 0.005;
        struct reb_particle pt2 = {0};
        pt2 = reb_tools_orbit_to_particle(r->G, star, 0, a2, 0, inc2, 0, 0,phi2);
        pt2.r 		= 0.04;
        pt2.id       = 4;
        reb_add(r, pt2);
        fprintf(stderr,"\033[1m Note!\033[0m planetesimal_2 is on!\n");
        r->N_planetesimals++;
    }
    
    system("rm -vf energy.txt");
    
    //Integrate!
    e0 = reb_tools_energy(r);
    //e0 = calc_Etot(planetesimal_mass, r);
    reb_integrate(r, 100);
    
    //********setup calculations and integrate!*********************
}

int mini_on = 0;

void heartbeat(struct reb_simulation* r){
    //Simple algorithm which starts using mini and updating global after 50 years
    if(r->t > 50 && mini_on == 0){
        s->t = r->t;
        struct reb_particle* global = r->particles;
        struct reb_particle* mini = s->particles;
        for(int i=0; i<r->N; i++) mini[i] = global[i];
        
        if(planetesimal_1 == 1){
            for(int i=0;i<r->N;i++){
                if(global[i].id == p1_id){
                    struct reb_particle pt = global[i];
                    reb_add(s,pt);
                    s->N_planetesimals++;
                    fprintf(stderr,"\033[1m Note!\033[0m Added planetesimal_1 to mini simulation!\n");
                    break;
                }
            }
            
        }
        
        printf("%d %d \n\n\n",r->N, s->N);
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
    //double e1 = calc_Etot(planetesimal_mass, r);
    fprintf(f,"%e %e\n",r->t,fabs((e0-e1)/e0));
    fclose(f); 
    //output stuff
    
    reb_output_timing(r, 0);
}

double calc_Etot(double planetesimal_mass, struct reb_simulation* a){
    double m1,m2;
    const int N = a->N;
    //const int N_active = a->N_active;
    const int N_active = a->N - a->N_planetesimals;
    const double G = a->G;
    double L = 0, U = 0, K = 0;
    struct reb_particle* const particles = a->particles;
    for(int i=0;i<N;i++){
        struct reb_particle par = particles[i];
        if(i < N_active) m1 = par.m; else m1 = planetesimal_mass;
        const double dvx = par.vx;
        const double dvy = par.vy;
        const double dvz = par.vz;
        const double dx = par.x;
        const double dy = par.y;
        const double dz = par.z;
        
        //L_tot = m*(r x v)
        const double hx = dy*dvz - dz*dvy;
        const double hy = dz*dvx - dx*dvz;
        const double hz = dx*dvy - dy*dvx;
        L += m1*sqrt ( hx*hx + hy*hy + hz*hz );
        
        //E_tot
        K += 0.5*m1*(dvx*dvx + dvy*dvy + dvz*dvz);
        if(i<N_active){//ignore dE/dx = forces between planetesimals
            for(int j=i+1;j<N;j++){
                struct reb_particle par2 = particles[j];
                if(j < N_active) m2 = par2.m; else m2 = planetesimal_mass;
                double ddx = dx - par2.x;
                double ddy = dy - par2.y;
                double ddz = dz - par2.z;
                U -= G*m1*m2/sqrt(ddx*ddx + ddy*ddy + ddz*ddz);
            }
        }
    }
    
    return K + U;
}

void planetesimal_forces(struct reb_simulation *a){
    const double G = a->G;
    const int N = a->N;
    //const int N_active = a->N_active;
    const int N_active = a->N;// - a->N_planetesimals;
    struct reb_particle* const particles = a->particles;
    
    const double Gm1 = G*planetesimal_mass;
    for(int i=0;i<N_active;i++){
        struct reb_particle* body = &(particles[i]);
        for(int j=N_active;j<N;j++){//add planetesimal forces to massive bodies
            struct reb_particle p = particles[j];
            
            const double dx = body->x - p.x;
            const double dy = body->y - p.y;
            const double dz = body->z - p.z;
            
            const double rijinv = 1.0/sqrt(dx*dx + dy*dy + dz*dz);
            const double ac = -Gm1*rijinv*rijinv*rijinv;  //force/mass = acceleration
            
            body->ax += ac*dx;    //perturbation on planets due to planetesimals.
            body->ay += ac*dy;
            body->az += ac*dz;
        }
    }
}