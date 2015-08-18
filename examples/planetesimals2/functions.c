//
//  functions.c
//  
//
//  Created by Ari Silburt on 2015-06-12.
//
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <sys/time.h>
#include "functions.h"
#include "../../src/rebound.h"

void legend(char* planetdir, char* legenddir, struct reb_simulation* r, double tmax, int N_active, int N, double m_planetesimal, double total_planetesimal_mass, double inner, double outer, double powerlaw, double mp, double a, double e, double Ms, double nrhill, double drh){
    
    system("rm -v output/orbit*.txt");
    
    char* us = "_";
    char* txt = ".txt";
    char* intgrtr;
    
    char str[100] = {0};
    if(r->integrator == REB_INTEGRATOR_HYBRID){
        intgrtr = "HYBRID";
        strcat(str, intgrtr);
        strcat(str, us);
        int hybrid_rint = (int) nrhill;
        char str3[15];
        sprintf(str3, "%d", hybrid_rint);
        strcat(str, str3);
        strcat(str, us);
        int drint = (int) drh;
        char str2[15];
        sprintf(str2, "%d", drint);
        strcat(str, str2);
        
        //char strtime[10];
        //sprintf(strtime, "%d", hybrid_rint);
        
    } else{ //pure IAS15 or WH
        if(r->integrator==REB_INTEGRATOR_IAS15) intgrtr = "IAS15"; else intgrtr = "WHFAST";
        char* teq = "_t=";
        strcat(str, intgrtr);
        strcat(str, teq);
        char dtstr[15];
        sprintf(dtstr, "%.0f", tmax);
        strcat(str, dtstr); //planet directory
    }
    
    strcat(legenddir, str);
    strcat(planetdir, str);
    strcat(planetdir, txt);
    
    char* file = "Properties.txt";
    strcat(legenddir, us);
    strcat(legenddir, file);
    FILE *ff;
    ff=fopen(legenddir, "w");
    fprintf(ff,"General:\ndt, tmax,  N_active, N_Rhill, dRHill, hybrid_switch_ratio, Integrator\n");
    fprintf(ff,"%f,%.1f,%d,%f,%f,%f,%s \n\n",r->dt,tmax,N_active,nrhill,drh,r->ri_hybrid.switch_ratio,intgrtr);
    fprintf(ff,"Planet/Star:\nplanet mass, semi-major axis, e_initial, Stellar Mass\n");
    fprintf(ff,"%f,%f,%f,%f\n\n",mp,a,e,Ms);
    fprintf(ff,"Planetesimal:\nN_planetesimals, Mtot_planetsimal, m_planetesimal, planetesimal boundary conditions: inner/outer edge, powerlaw\n");
    fprintf(ff,"%d,%f,%.11f,%f,%f,%f\n\n",N-N_active,total_planetesimal_mass, m_planetesimal, inner, outer, powerlaw);
    fclose(ff);
    
    //remove any old planet files if there are
    char rmv[100] = {0};
    char* rm_v = "rm -v ";
    strcat(rmv, rm_v);
    strcat(rmv, planetdir);
    system(rmv);
    
}

double calc_dt(struct reb_simulation* r, double mp, double Ms, double a, double N_Rhill, double dRHill){
    if(dRHill > N_Rhill){
        printf("\033[1mWarning!\033[0m dRhill !> N_RHill. Setting dRhill = N_Rhill/2 \n");
        dRHill = 0.5*N_Rhill;
    }
    double e_max = 0.3;  //max hypothesized eccentricity that the planet/esimals could have
    double Hill = a*(1 - e_max)*pow(mp/(3*Ms),1./3.);
    //double r2_E = N_Rhill*N_Rhill*Hill*Hill;
    //r->ri_hybrid.switch_ratio = Ms*r2_E/(mp*a*a*1.21);
    r->ri_hybrid.switch_ratio = N_Rhill;    //looks like ratio is rij squared, so this should be too?
    double vmax = sqrt(r->G*(Ms + mp)*(1 + e_max)/(a*(1 - e_max)));   //peri speed
    double dt = dRHill*Hill/vmax;
    printf("timesetep is dt = %f, hybrid_switch_ratio=%f \n",dt,r->ri_hybrid.switch_ratio);
    
    return dt;
}

void calc_ELtot(double* Etot, double* Ltot, double planetesimal_mass, struct reb_simulation* r){
    //first need to fill arrays
    //struct particle com = particles[0];
    double m1,m2;
    int N_active = r->N_active, N = r->N;
    const double G = r->G;
    double L = 0, E = 0;
    struct reb_particle* const particles = r->particles;
    for(int i=0;i<N;i++){
        struct reb_particle* par = &(particles[i]); //planet occupies first slot.
        if(i > N_active - 1) m1 = planetesimal_mass; else m1 = par->m;
        const double dvx = par->vx;
        const double dvy = par->vy;
        const double dvz = par->vz;
        const double dx = par->x;
        const double dy = par->y;
        const double dz = par->z;
        
        //L_tot = m*(r x v)
        const double hx = dy*dvz - dz*dvy;
        const double hy = dz*dvx - dx*dvz;
        const double hz = dx*dvy - dy*dvx;
        L += m1*sqrt ( hx*hx + hy*hy + hz*hz );
        
        //E_tot
        E += 0.5*m1*(dvx*dvx + dvy*dvy + dvz*dvz);
        for(int j=i+1;j<N;j++){
            struct reb_particle* par2 = &(particles[j]);
            if(j > N_active - 1) m2 = planetesimal_mass; else m2 = par2->m;
            double ddx = dx - par2->x;
            double ddy = dy - par2->y;
            double ddz = dz - par2->z;
            E -= G*m1*m2/sqrt(ddx*ddx + ddy*ddy + ddz*ddz);
        }
    }
    *Etot = E;
    *Ltot = L;
}

//Calculates 'a' and 'e' of planet each output.
void calc_ae(double* a, double* e, struct reb_simulation* r){
    struct reb_particle* const particles = r->particles;
    struct reb_particle com = particles[0];
    struct reb_particle* par = &(particles[1]); //planet occupies first slot.
    const double G = r->G;
    const double m = par->m;
    const double mu = G*(com.m + m);
    const double dvx = par->vx-com.vx;
    const double dvy = par->vy-com.vy;
    const double dvz = par->vz-com.vz;
    const double dx = par->x-com.x;
    const double dy = par->y-com.y;
    const double dz = par->z-com.z;
    
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
}

void planetesimal_forces(struct reb_simulation *a){
    const double G = a->G;
    const int N = a->N;
    const int N_active = a->N_active;
    struct reb_particle* const particles = a->particles;
    struct reb_particle com = particles[0];
    const double Gm1 = G*planetesimal_mass;
    for(int i=1;i<N_active;i++){
        struct reb_particle* planet = &(particles[i]);
        const double x = planet->x-com.x;
        const double y = planet->y-com.y;
        const double z = planet->z-com.z;
        for(int j=N_active;j<N;j++){//add forces to planet
            struct reb_particle* p = &(particles[j]);
            //if(p->id == 4) continue;
            const double xp = p->x-com.x;
            const double yp = p->y-com.y;
            const double zp = p->z-com.z;
            
            const double dx = x - xp;
            const double dy = y - yp;
            const double dz = z - zp;
            const double rinv = 1./sqrt( dx*dx + dy*dy + dz*dz );
            const double ac = Gm1*rinv*rinv*rinv;  //force/mass = acceleration
            
            planet->ax += ac*dx;    //perturbation on planet due to planetesimals
            planet->ay += ac*dy;
            planet->az += ac*dz;
        }
    }
}

int check_for_encounter(struct reb_simulation* const r){
    const int N = r->N;
    const int N_active = r->N_active;
    const int N_var = r->N_var;
    struct reb_particle* restrict const particles = r->particles;
    struct reb_particle p0 = particles[0];
    const int _N_active = ((N_active==-1)?N:N_active)- N_var;
    const int _N_real   = N - N_var;
    int index_of_encounter = 0;
    for (int i=1; i<_N_active; i++){
        struct reb_particle pi = particles[i];
        const double dxi = p0.x - pi.x;
        const double dyi = p0.y - pi.y;
        const double dzi = p0.z - pi.z;
        const double r0i2 = dxi*dxi + dyi*dyi + dzi*dzi;
        const double rhi = r0i2*pow((pi.m/(p0.m*3.)), 2./3.);
        
        for (int j=1; j<_N_real; j++){
            if (i==j) continue;
            
            struct reb_particle pj = particles[j];
            
            const double dx = pi.x - pj.x;
            const double dy = pi.y - pj.y;
            const double dz = pi.z - pj.z;
            const double rij2 = dx*dx + dy*dy + dz*dz;
            const double dxj = p0.x - pj.x;
            const double dyj = p0.y - pj.y;
            const double dzj = p0.z - pj.z;
            const double r0j2 = dxj*dxj + dyj*dyj + dzj*dzj;
            const double rhj = r0j2*pow((pj.m/(p0.m*3.)), 2./3.);
            
            const double ratio = rij2/(rhi+rhj);
            
            if(ratio<r->ri_hybrid.switch_ratio){
                index_of_encounter = j;
                goto outer;
            }
        }
    }
outer:;
    return index_of_encounter;
}

//initialize mini-simulation for close encounters
void ini_mini(struct reb_simulation* const r, struct reb_simulation* s){
    s->N_active = r->N_active;
    s->ri_hybrid.switch_ratio = r->ri_hybrid.switch_ratio;
    s->integrator = REB_INTEGRATOR_IAS15;
    s->additional_forces = planetesimal_forces;
    s->exact_finish_time = 1;
    s->dt = r->dt;
    
    //add massive particles
    struct reb_particle* restrict const particles = r->particles;
    for(int k=0; k<s->N_active; k++){
        struct reb_particle p = particles[k];
        reb_add(s,p);
    }
    //add dummy test particle
    struct reb_particle pt = particles[r->N_active];
    reb_add(s,pt);
    reb_move_to_com(s);
}

void update_mini(struct reb_simulation* const r, struct reb_simulation* s, int encounter_index){
    s->t = r->t;
    
    struct reb_particle* const global = r->particles;
    struct reb_particle* mini = s->particles;
    
    //update massive particles
    for(int i=0; i<s->N_active; i++) mini[i] = global[i];
    
    //update test particle inside Hill Sphere
    mini[s->N_active] = global[encounter_index];
}

void update_global(struct reb_simulation* const s, struct reb_simulation* r, int encounter_index){
    
    struct reb_particle* global = r->particles;
    struct reb_particle* const mini = s->particles;
    
    //update massive particles
    for(int i=0; i<s->N_active; i++) global[i] = mini[i];
        
    //update test particle
    global[encounter_index] = mini[s->N_active];
}

/*
void close_encounter(struct reb_simulation* r){
    double ratio = 0;
    int encounter_index = check_for_encounter(r);
    if(encounter_index != 0){//create new rebound simulation
        struct reb_simulation* s = reb_create_simulation();
        s->N_active = r->N_active;
        s->ri_hybrid.switch_ratio = r->ri_hybrid.switch_ratio;
        s->integrator = REB_INTEGRATOR_IAS15;
        s->exact_finish_time = 1;
        const double timestep = r->dt;
        
        struct reb_particle* restrict const particles = r->particles;
        struct reb_particle p0 = particles[0];
        reb_add(s,p0);
        
        for(int k=1; k<s->N_active; k++){//add massive particles
            struct reb_particle p = particles[k];
            reb_add(s,p);
        }
        //copy contents of particle into new simulation, and then put the particle
        //(from the global sim) temporarily out of harms way.
        int dt_counter=0;   //count # of while loops
        struct reb_particle pt = particles[encounter_index];
        pt.id = 3;
        reb_add(s,pt);
        //printf("\n old params: %f,%f,%f,%f,%f,%f,   ",pt.x,pt.y,pt.vx,pt.vy,pt.ax,pt.ay);
        inactive_particle(&pt,p0.m,r->G);
        //printf("new params: %f,%f,%f,%f,%f,%f \n",pt.x,pt.y,pt.vx,pt.vy,pt.ax,pt.ay);
        
        //different index here, only 3 particles
        while(encounter_index != 0){
            dt_counter++;
            reb_integrate(s, dt_counter*timestep);
            planetesimal_forces(s,r,0); //needs to be changed to 1********************
            encounter_index = check_for_encounter(s,&ratio);
            //struct reb_particle* const p = s->particles;
            //printf("planet1: dt_counter=%d, ax,ay,az = %f,%f,%f, time=%f \n",dt_counter, p[1].x, p[1].y, p[1].z, s->t);
        }
        *CE_exit_time = s->t + r->t;
    }
}
*/

/*
void inactive_particle(struct reb_particle* pt, double Ms, double G){
    double a = 4;   //temporarily put particle way out of harms way.
    
    double phi = 0;
    pt->x 		= a*cos(phi);
    pt->y 		= a*sin(phi);
    pt->z 		= 0;
    double vkep = sqrt(G*Ms/a);
    pt->vx 		= -vkep * sin(phi);
    pt->vy 		= vkep * cos(phi);
    pt->vz       = 0;
    pt->ax       = 0;
    pt->ay       = 0;
    pt->az       = 0;
}*/