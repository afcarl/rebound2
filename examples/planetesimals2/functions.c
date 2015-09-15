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

void legend(char* planetdir, char* legenddir, struct reb_simulation* r, double tmax, double m_planetesimal, double total_planetesimal_mass, double inner, double outer, double powerlaw, double mp, double a, double e, double Ms, double drh, int HYBRID_ON){
    
    int N_active = r->N_active, N = r->N;
    
    system("rm -v output/orbit*.txt");
    
    char* us = "_";
    char* txt = ".txt";
    char* intgrtr;
    
    char str[100] = {0};
    if(r->integrator == REB_INTEGRATOR_WHFAST){
        if(HYBRID_ON == 1)intgrtr = "HYBRID"; else intgrtr = "WHFAST";
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

void calc_ELtot(double* Etot, double* Ktot, double* Utot, double* Ltot, double planetesimal_mass, struct reb_simulation* r){
    double m1,m2;
    int N_active = r->N_active, N = r->N;
    const double G = r->G;
    double L = 0, U = 0, K = 0;
    struct reb_particle* const particles = r->particles;
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
            struct reb_particle p = particles[j];
            
            const double dx = body->x - p.x;
            const double dy = body->y - p.y;
            const double dz = body->z - p.z;
           
            const double rijinv = 1.0/sqrt(dx*dx + dy*dy + dz*dz);
            const double ac = Gm1*rijinv*rijinv*rijinv;  //force/mass = acceleration
            
            body->ax -= ac*dx;    //perturbation on planets due to planetesimals.
            body->ay -= ac*dy;
            body->az -= ac*dz;
        }
    }
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
    
    //reb_move_to_com(s);         //before IAS15 simulation starts, move to COM
    //move_to_com_with_planetesimals(s);
}

//collect the id/array number of all planetesimals involved in a close encounter
void check_for_encounter(struct reb_simulation* const r, int* N_encounters){
    const int N = r->N;
    const int N_active = r->N_active;
    struct reb_particle* restrict const particles = r->particles;
    struct reb_particle p0 = particles[0];
    int num_encounters = 0;
    for (int i=1; i<N_active; i++){
        struct reb_particle body = particles[i];
        const double dxi = p0.x - body.x;
        const double dyi = p0.y - body.y;
        const double dzi = p0.z - body.z;
        const double r0i2 = dxi*dxi + dyi*dyi + dzi*dzi;
        const double rhi = r0i2*Hill[i];
        //const double rhi = r0i2*pow((body.m/(p0.m*3.)), 2./3.); //can make this faster later
        
        for (int j=1; j<N; j++){
            if(i==j) continue;
            struct reb_particle pj = particles[j];
            
            const double dxj = p0.x - pj.x;
            const double dyj = p0.y - pj.y;
            const double dzj = p0.z - pj.z;
            const double r0j2 = dxj*dxj + dyj*dyj + dzj*dzj;
            const double rhj = r0j2*Hill[j];
            //const double rhj = r0j2*pow((pj.m/(p0.m*3.)), 2./3.); //can make this faster later
            
            const double dx = body.x - pj.x;
            const double dy = body.y - pj.y;
            const double dz = body.z - pj.z;
            const double rij2 = dx*dx + dy*dy + dz*dz;
            const double ratio = rij2/(rhi+rhj);
            
            if(ratio<r->ri_hybrid.switch_ratio){
                num_encounters++;
                if(num_encounters == 1) encounter_index[0] = pj.id;
                else if(num_encounters > 1){//multiple close encounters
                    encounter_index = realloc(encounter_index,num_encounters*sizeof(int));
                    encounter_index[num_encounters - 1] = pj.id;
                }
            }
            
            if(rij2 < 5e-8){
                fprintf(stderr,"\n\033[1mAlert!\033[0m Particle/Planet collision should have happened.\n");
            }
        }
    }
    *N_encounters = num_encounters;
}

//Just after mini has been integrated up to r->t, update global.
void update_global(struct reb_simulation* const s, struct reb_simulation* r, int N_encounters_previous, int N_encounters){
    int N_active = s->N_active;
    struct reb_particle* global = r->particles;
    struct reb_particle* const mini = s->particles;
    
    //update massive and planetesimal particles
    for(int i=0; i<N_active; i++) global[i] = mini[i];  //massive particles, always in same order
    for(int j=0; j<N_encounters_previous; j++){
        _Bool particle_update = 0;
        int PEI = previous_encounter_index[j];          //encounter index == global[EI].id
        for(int k=0;particle_update==0 && k<N_encounters_previous;k++){
            int mini_index = N_active + k;
            if(mini[mini_index].id == PEI){
                //printf("mini_index=%d,id=%d, s->N_active=%d\n",mini_index,mini[mini_index].id,N_active);
                global[PEI] = mini[mini_index];
                particle_update = 1;
            }
        }
        
        if(particle_update == 0){
            fprintf(stderr,"\n\033[1mAlert!\033[0m Particle %d couldn't be found in update_global. Exiting. \n",PEI);
            printf("N_encounters=%d,size=%lu,",N_encounters,sizeof(encounter_index)/sizeof(encounter_index[0]));
            for(int i=0;i<N_encounters;i++) printf("EI(%d)=%d,",i,encounter_index[i]);
            printf("\n");
            printf("N_encounters_previous=%d,size=%lu,",N_encounters_previous,sizeof(previous_encounter_index)/sizeof(previous_encounter_index[0]));
            for(int i=0;i<N_encounters_previous;i++) printf("PEI(%d)=%d,",i,previous_encounter_index[i]);
            printf("\n");
            printf("Mini:");
            for(int i=0;i<N_encounters_previous;i++) printf("mini[%d].id=%d,",i,mini[N_active+i].id);
            printf("\n");
            exit(0);
        }
    }
    
}

void add_or_subtract_particles(struct reb_simulation* r, struct reb_simulation* s, int N_encounters, int N_encounters_previous, int dN){
    int N_active = s->N_active;
    struct reb_particle* mini = s->particles;
    struct reb_particle* global = r->particles;
    
    //check to add particles
    for(int i=0;i<N_encounters;i++){
        _Bool index_found = 0;
        int EI = encounter_index[i];
        for(int j=0;index_found == 0 && j<N_encounters_previous;j++){
            if(EI == previous_encounter_index[j]) index_found = 1;
        }
        if(index_found == 0){//couldn't find index, add particle
            struct reb_particle pt = global[EI];
            reb_add(s,pt);
            N_encounters_tot++;
            printf("particle %d added. dN == %d, N_close_encounters=%d\n",EI,dN,N_encounters);
        }
    }
    
    //check to remove particles
    for(int i=0;i<N_encounters_previous;i++){
        _Bool index_found = 0;
        int PEI = previous_encounter_index[i];
        for(int j=0;index_found == 0 && j<N_encounters;j++){
            if(PEI == encounter_index[j]) index_found = 1;
        }
        if(index_found == 0){//couldn't find index, find particle in sim and remove
            int removed_particle = 0;
            for(int k=0;removed_particle==0 && k<N_encounters_previous;k++){
                if(mini[k+N_active].id == PEI){
                    removed_particle = reb_remove(s,k+N_active,1);    //remove particle
                    printf("particle %d leaving. dN == %d, N_close_encounters=%d.\n",PEI,dN,N_encounters);
                }
            }
        }
    }

}

//transfer values from encounter_index to previous_encounter_index
void update_encounter_indices(int* N_encounters, int* N_encounters_previous){
    int size;
    if(*N_encounters == 0) size = 1; else size = *N_encounters;
    
    previous_encounter_index = realloc(previous_encounter_index,size*sizeof(int));
    
    if(*N_encounters == 0) previous_encounter_index[0] = 0; else {
        for(int i=0;i<*N_encounters;i++) previous_encounter_index[i] = encounter_index[i];
    }
    
    //reset encounter index
    encounter_index = realloc(encounter_index,sizeof(int));   //reset to single element
    encounter_index[0] = 0;
    
    //reset encounter counters.
    *N_encounters_previous = *N_encounters;
    *N_encounters = 0;
}

void clock_finish(clock_t timer, int N_encounters, char* legenddir){
    timer = clock() - timer;
    FILE *ff;
    ff=fopen(legenddir, "a");
    double result = ((float)timer)/CLOCKS_PER_SEC;
    fprintf(ff,"Elapsed simulation time is %f s, with %d close encounters.\n",result,N_encounters);
    printf("\n\nSimulation complete. Elapsed simulation time is %f s, with %d close encounters.\n\n",result,N_encounters);
}

struct reb_particle get_com_with_planetesimals(struct reb_simulation* r){
    struct reb_particle com = {.m=0, .x=0, .y=0, .z=0, .vx=0, .vy=0, .vz=0};
    const int N = r->N;
    const int N_active = r->N_active;
    int m;
    struct reb_particle* restrict const particles = r->particles;
    for (int i=0;i<N;i++){
        struct reb_particle p = particles[i];
        if(i>=N_active) m = planetesimal_mass; else m = p.m;
        com.x = com.x*com.m + p.x*m;
        com.y = com.y*com.m + p.y*m;
        com.z = com.z*com.m + p.z*m;
        com.vx = com.vx*com.m + p.vx*m;
        com.vy = com.vy*com.m + p.vy*m;
        com.vz = com.vz*com.m + p.vz*m;
        com.m += m;
        if(com.m > 0){
            com.x /= com.m;
            com.y /= com.m;
            com.z /= com.m;
            com.vx /= com.m;
            com.vy /= com.m;
            com.vz /= com.m;
        }
    }
    return com;
}

void move_to_com_with_planetesimals(struct reb_simulation* const r){
    const int N = r->N;
    struct reb_particle* restrict const particles = r->particles;
    struct reb_particle com = get_com_with_planetesimals(r);
    for (int i=0;i<N;i++){
        particles[i].x  -= com.x;
        particles[i].y  -= com.y;
        particles[i].z  -= com.z;
        particles[i].vx -= com.vx;
        particles[i].vy -= com.vy;
        particles[i].vz -= com.vz;
    }
}