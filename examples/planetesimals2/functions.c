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
#include "../../src/integrator_whfast.h"

void legend(char* planetdir, char* legenddir, struct reb_simulation* r, double tmax, double m_planetesimal, double total_planetesimal_mass, double inner, double outer, double powerlaw, double mp, double a, double e, double Ms, double drh){
    
    int N_active = r->N_active, N = r->N;
    
    system("rm -v output/orbit*.txt");
    
    char* us = "_";
    char* txt = ".txt";
    char* intgrtr;
    
    char str[100] = {0};
    if(r->integrator == REB_INTEGRATOR_WHFAST){
        intgrtr = "HYBRID";
        strcat(str, intgrtr);
        strcat(str, us);
        int hybrid_rint = (int) r->ri_hybrid.switch_ratio;
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
        
    } else{ //pure IAS15
        intgrtr = "IAS15";
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
    fprintf(ff,"General:\ndt, tmax, N_active, N_Rhill (ri_hybrid.switch_ratio), dRHill, Integrator\n");
    fprintf(ff,"%f,%.1f,%d,%f,%f,%s \n\n",r->dt,tmax,N_active,r->ri_hybrid.switch_ratio,drh,intgrtr);
    fprintf(ff,"Planet/Star:\nplanet mass, semi-major axis, e_initial, Stellar Mass\n");
    fprintf(ff,"%f,%f,%f,%f\n\n",mp,a,e,Ms);
    fprintf(ff,"Planetesimal:\nN_planetesimals, Mtot_planetsimal, m_planetesimal, planetesimal boundary conditions: inner/outer edge, powerlaw\n");
    fprintf(ff,"%d,%.10f,%.15f,%f,%f,%f\n\n",N-N_active,total_planetesimal_mass, m_planetesimal, inner, outer, powerlaw);
    fclose(ff);
    
    //remove any old planet files if there are
    char rmv[100] = {0};
    char* rm_v = "rm -v ";
    strcat(rmv, rm_v);
    strcat(rmv, planetdir);
    system(rmv);
    
}

double calc_dt(struct reb_simulation* r, double mp, double Ms, double a, double dRHill){
    
    if(dRHill > r->ri_hybrid.switch_ratio){
        printf("\033[1mWarning!\033[0m dRhill !> N_RHill. Setting dRhill = N_Rhill/2 \n");
        dRHill = 0.5*r->ri_hybrid.switch_ratio;
    }
    double e_max = 0.3;  //max hypothesized eccentricity that the planet/esimals could have
    double Hill = a*(1 - e_max)*pow(mp/(3*Ms),1./3.);
    double vmax = sqrt(r->G*(Ms + mp)*(1 + e_max)/(a*(1 - e_max)));   //peri speed
    double dt = dRHill*Hill/vmax;
    printf("timesetep is dt = %f, ri_hybrid.switch_ratio=%f \n",dt,r->ri_hybrid.switch_ratio);
    
    return dt;
}

void calc_ELtot(double* Etot, double* Ltot, double planetesimal_mass, struct reb_simulation* r){
    double m1,m2;
    int N_active = r->N_active, N = r->N;
    const double G = r->G;
    double L = 0, E = 0;
    struct reb_particle* const particles = r->particles;
    for(int i=0;i<N;i++){
        struct reb_particle par = particles[i]; //planet occupies first slot.
        if(i >= N_active) m1 = planetesimal_mass; else m1 = par.m;
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
        E += 0.5*m1*(dvx*dvx + dvy*dvy + dvz*dvz);
        for(int j=i+1;j<N;j++){
            struct reb_particle par2 = particles[j];
            if(j >= N_active) m2 = planetesimal_mass; else m2 = par2.m;
            double ddx = dx - par2.x;
            double ddy = dy - par2.y;
            double ddz = dz - par2.z;
            E -= G*m1*m2/sqrt(ddx*ddx + ddy*ddy + ddz*ddz);
        }
    }
        
    *Etot = E;
    *Ltot = L;
}

//Calculates 'a' and 'e' of planet each output.
void calc_ae(double* a, double* e, double* d_out, struct reb_simulation* r, int i){
    struct reb_particle* const particles = r->particles;
    struct reb_particle com = particles[0];
    struct reb_particle* par = &(particles[i]); //output planets only.
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
        for(int j=N_active;j<N;j++){//add forces to massive bodies
            struct reb_particle* p = &(particles[j]);
            
            const double dx = body->x - p->x;
            const double dy = body->y - p->y;
            const double dz = body->z - p->z;
           
            const double rijinv = 1.0/sqrt(dx*dx + dy*dy + dz*dz);
            const double ac = Gm1*rijinv*rijinv*rijinv;  //force/mass = acceleration
            
            body->ax -= ac*dx;    //perturbation on planets due to planetesimals.
            body->ay -= ac*dy;
            body->az -= ac*dz;
        }
    }
}

void check_for_encounter(struct reb_simulation* const r, int** index_of_encounters, int* N_encounters){
    const int N = r->N;
    const int N_active = r->N_active;
    struct reb_particle* restrict const particles = r->particles;
    struct reb_particle p0 = particles[0];
    int num_encounters = 0;
    for (int i=1; i<N_active; i++){
        struct reb_particle pi = particles[i];
        const double dxi = p0.x - pi.x;
        const double dyi = p0.y - pi.y;
        const double dzi = p0.z - pi.z;
        const double r0i2 = dxi*dxi + dyi*dyi + dzi*dzi;
        const double rhi = r0i2*pow((pi.m/(p0.m*3.)), 2./3.); //can make this faster later
        
        for (int j=1; j<N; j++){
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
            const double rhj = r0j2*pow((pj.m/(p0.m*3.)), 2./3.); //can make this faster later
            
            const double ratio = rij2/(rhi+rhj);
            
            if(ratio<r->ri_hybrid.switch_ratio){
                num_encounters++;
                if(num_encounters == 1){
                    *index_of_encounters = malloc(num_encounters*sizeof(int));
                    *index_of_encounters[0] = j;
                }
                else if(num_encounters > 1){//multiple close encounters
                    int* tmpindex = realloc(*index_of_encounters,num_encounters*sizeof(int));
                    if(!tmpindex){
                        fprintf(stderr,"\n\033[1mAlert!\033[0m Could not reallocate close encounter particle array. Not tracking particle. Exiting.\n");
                        exit(0);
                    } else {
                        tmpindex[num_encounters - 1] = j;
                        *index_of_encounters = tmpindex;
                    }
                }
            }
            
            if(rij2 < 5e-8){
                fprintf(stderr,"\n\033[1mAlert!\033[0m Particle/Planet collision should have happened.\n");
            }
        }
    }
    *N_encounters = num_encounters;
}

//initialize mini-simulation for close encounters
void ini_mini(struct reb_simulation* const r, struct reb_simulation* s){
    s->N_active = r->N_active;
    s->integrator = REB_INTEGRATOR_IAS15;
    s->additional_forces = planetesimal_forces;
    s->exact_finish_time = 1;
    s->dt = r->dt;
    
    struct reb_particle* restrict const particles = r->particles;
    for(int k=0; k<s->N_active; k++){
        struct reb_particle p = particles[k];
        reb_add(s,p);
    }
    
    reb_move_to_com(s);         //before IAS15 simulation starts, move to COM
}

void add_mini(struct reb_simulation* const r, struct reb_simulation* s, int* encounter_index, int N_encounters){
    int N_active = s->N_active;
    struct reb_particle* const global = r->particles;
    struct reb_particle* mini = s->particles;
    
    if(N_encounters == 1){//first update in a while, sync up times, update massive bodies
        s->t = r->t;    
        for(int i=0; i<N_active; i++) mini[i] = global[i];   //massive
    }
    
    //add newest test particle and update it
    int mini_index = N_active-1+N_encounters;       //newest particle crossing Hill radius
    struct reb_particle pt = particles[mini_index];
    reb_add(s,pt);
    mini[mini_index] = global[encounter_index[N_encounters - 1]];            //test particle
}

void update_global(struct reb_simulation* const s, struct reb_simulation* r, int* encounter_index, int N_encounters){
    int N_active = s->N_active;
    struct reb_particle* global = r->particles;
    struct reb_particle* const mini = s->particles;

    //update particles
    for(int i=0; i<N_active; i++) global[i] = mini[i];   //massive
    for(int j=0; j<N_encounters; j++) global[encounter_index[j]] = mini[N_active + j];     //update test particle(s)
    
    //const double dx = mini[1].x - mini[2].x;
    //const double dy = mini[1].y - mini[2].y;
    //const double dz = mini[1].z - mini[2].z;
    //double rij2 = dx*dx + dy*dy + dz*dz;
    //printf("t=%f, particle-planet distance = %.10f\n",r->t, sqrt(rij2));
}

void clock_finish(clock_t timer, int N_encounters, char* legenddir){
    timer = clock() - timer;
    FILE *ff;
    ff=fopen(legenddir, "a");
    double result = ((float)timer)/CLOCKS_PER_SEC;
    fprintf(ff,"Elapsed simulation time is %f s, with %d close encounters.\n",result,N_encounters);
    printf("\n\nSimulation complete. Elapsed simulation time is %f s, with %d close encounters.\n\n",result,N_encounters);
}