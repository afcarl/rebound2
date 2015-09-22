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
void calc_ELtot(double* Etot, double* Ktot, double* Utot, double* Ltot, double planetesimal_mass, struct reb_simulation* r);
void calc_ae(double* a, double* e, double* d, struct reb_simulation* r, int i, double t);
void planetesimal_forces(struct reb_simulation *a);

double tmax, planetesimal_mass, CE_exit_time = 0, E_ini, K_ini, U_ini, L_ini, n_output;
int mini_on = 0, planetesimal_1, output_counter = 0, p1_id;
//int *encounter_index; int *previous_encounter_index; double* Hill = NULL;
char plntdir[200] = "output/planet_", lgnddir[200] = "output/planet_";
struct reb_simulation* s;

int main(int argc, char* argv[]){
//********switches******************
    planetesimal_1 = 0;     //=1 to turn on, 0=off
    int planetesimal_2 = 0;
    int turn_planetesimal_forces_on = 0; //=1 to turn on
    
//********setup values******************
    tmax = 100;
    double M_planetesimals = 3e-6/10;
    int N_planetesimals = 2;
    planetesimal_mass = M_planetesimals/N_planetesimals;
    
    struct reb_simulation* r = reb_create_simulation();
	// Setup constants
    r->integrator	= 1;    //REB_INTEGRATOR_IAS15 = 0, WHFAST = 1, WH=3, HYBRID = 5
	r->collision	= REB_COLLISION_NONE;
	r->boundary     = REB_BOUNDARY_OPEN;
	r->heartbeat	= heartbeat;
    if(turn_planetesimal_forces_on==1)r->additional_forces = planetesimal_forces;
    r->ri_hybrid.switch_ratio = 5;     //# hill radii for boundary between switch. Try 3?
    //r->usleep   = 1000; //larger the number, slower OpenGL simulation
	
    // Other constants
    n_output = 10000;
    double boxsize = 5;
	reb_configure_box(r, boxsize, 1, 1, 1);
    
//********bodies*********************
	//star
	struct reb_particle star = {0};
	star.m 		= 1;
	star.r		= 0.01;
    star.id     = 0;        // 0 = star
	reb_add(r, star);
    r->N_active = 1;

    //planet 1
    double a1=0.7, m1=5e-5, e1=0.01;
    struct reb_particle p1 = {0};
    p1 = reb_tools_orbit_to_particle(r->G, star, m1, a1, e1, reb_random_normal(0.0001), 0., 0., 0.);
    p1.id = 1;              //1 = planet
    reb_add(r, p1);
    r->N_active++;
    
    //planet 2
    double a2=1, m2=5e-5, e2=0.01;
    struct reb_particle p2 = {0};
    p2 = reb_tools_orbit_to_particle(r->G, star, m2, a2, e2, reb_random_normal(0.0001), 0., 0., 0.);
    p2.id = 1;              //2 = planet
    reb_add(r, p2);
    r->N_active++;
    
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
    }
    
//********setup calculations and integrate!*********************
    //dt
    r->dt = 0.004;
    
    //move to COM
    if(r->integrator != REB_INTEGRATOR_WH) reb_move_to_com(r);
    
    //Ini mini
    s = reb_create_simulation();    //initialize mini simulation (IAS15)
    s->N_active = r->N_active;
    s->integrator = 0; //REB_INTEGRATOR_IAS15 = 0, WHFAST = 1, WH=3, HYBRID = 5
    if(turn_planetesimal_forces_on==1){
       s->additional_forces = planetesimal_forces;
        fprintf(stderr,"\033[1m Note!\033[0m Planetesimal_forces are on!\n");
    }
    s->exact_finish_time = 1;
    s->dt = r->dt;
    struct reb_particle* restrict const particles = r->particles;
    for(int k=0; k<s->N_active; k++){
        struct reb_particle p = {0};
        p = particles[k];
        reb_add(s,p);
    }
    
    system("rm -v simple_test.txt");
    printf("r->N-1=%d (input for orbits.py)\n",r->N-1);
    //calc_ELtot(&E_ini, &K_ini, &U_ini, &L_ini, 0, r);
    
    //Integrate!
    reb_integrate(r, tmax);

//********setup calculations and integrate!*********************
}

void heartbeat(struct reb_simulation* r){
    //Simple algorithm which starts using mini and updating global after 50 years
    if(r->t > 50 && mini_on == 0){
        s->t = r->t;
        int N_active = r->N_active;
        struct reb_particle* global = r->particles;
        struct reb_particle* mini = s->particles;
        for(int i=0; i<N_active; i++) mini[i] = global[i];
        
        fprintf(stderr,"\n\033[1m Note!\033[0m Initializing mini simulation (IAS)\n");
        mini_on = 1;
        
        if(planetesimal_1 == 1){
            for(int i=0;i<r->N;i++){
                if(global[i].id == p1_id){
                    struct reb_particle pt = global[i];
                    reb_add(s,pt);
                    fprintf(stderr,"\033[1m Note!\033[0m Added planetesimal_1 to mini simulation!\n");
                    break;
                }
            }

        }
        
    } else if(mini_on == 1){
        reb_integrate(s, r->t);
        struct reb_particle* global = r->particles;
        struct reb_particle* mini = s->particles;
        
        //update massive bodies and planetesimal particle
        for(int i=0; i<s->N; i++) global[i] = mini[i];
    }
    
    //output stuff
    if(r->t > output_counter*tmax/n_output){
        output_counter++;
        double E_curr = 0, K_curr = 0, U_curr = 0, L_curr = 0, a_p = 0, d_p = 0, e_p = 0, t = r->t;
        calc_ELtot(&E_curr, &K_curr, &U_curr, &L_curr, planetesimal_mass, r); //calcs Etot all in one go.
        if(r->t <= r->dt){
            printf("r->t=%f,E_ini replacement\n",r->t);
            E_ini = E_curr;
            K_ini = K_curr;
            U_ini = U_curr;
            L_ini = L_curr;
        } else {
            for(int i=1;i<r->N;i++){
                calc_ae(&a_p, &e_p, &d_p, r, i, t);
                
                FILE *append;
                append=fopen("simple_test.txt", "a");
                fprintf(append,"%f,%.8f,%.8f,%.16f,%.16f,%.16f,%.16f,%.16f\n",t,a_p,e_p,fabs((E_ini - E_curr)/E_ini),fabs((K_ini - K_curr)/K_ini), fabs((U_ini - U_curr)/U_ini),fabs((L_ini - L_curr)/L_ini),d_p);
                fclose(append);
                E_curr = 0;
                L_curr = 0;
            }
        }
        
        reb_output_timing(r, 0);
    }
}

void calc_ELtot(double* Etot, double* Ktot, double* Utot, double* Ltot, double planetesimal_mass, struct reb_simulation* a){
    double m1,m2;
    int N_active = a->N_active, N = a->N;
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
    
    *Etot = K + U;
    *Ktot = K;
    *Utot = U;
    *Ltot = L;
}

//Calculates 'a' and 'e' of planet each output.
void calc_ae(double* a, double* e, double* d_out, struct reb_simulation* r, int i, double t){
    struct reb_particle* const particles = r->particles;
    struct reb_particle com = reb_get_com(r);
    struct reb_particle par = particles[i]; //output planets only.
    const double G = r->G;
    const double m = par.m;
    const double mu = G*(com.m + m);
    const double dvx = par.vx-com.vx;
    const double dvy = par.vy-com.vy;
    const double dvz = par.vz-com.vz;
    const double dx = par.x-com.x;
    const double dy = par.y-com.y;
    const double dz = par.z-com.z;
    
    const double vv = dvx*dvx + dvy*dvy + dvz*dvz;
    const double d = sqrt( dx*dx + dy*dy + dz*dz );    //distance
    const double dinv = 1./d;
    const double muinv = 1./mu;
    const double vr = (dx*dvx + dy*dvy + dz*dvz)*dinv;
    const double term1 = vv-mu*dinv;
    const double term2 = d*vr;
    const double ex = muinv*( term1*dx - term2*dvx );
    const double ey = muinv*( term1*dy - term2*dvy );
    const double ez = muinv*( term1*dz - term2*dvz );
    *e = sqrt(ex*ex + ey*ey + ez*ez);   // eccentricity
    *a = -mu/( vv - 2.*mu*dinv );
    *d_out = d;
    
}

void planetesimal_forces(struct reb_simulation *a){
    const double G = a->G;
    const int N = a->N;
    const int N_active = a->N_active;
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
