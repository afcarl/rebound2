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

void legend(char* planetdir, char* legenddir, char* xyz_check, char* CEprint, struct reb_simulation* r, double tmax, double m_planetesimal, double total_planetesimal_mass, int N_planetesimals, double inner, double outer, double powerlaw, double mp, double a, double e, double Ms, double drh, double epsilon, int seed, int HYBRID_ON){
    
    int N_active = r->N_active, N = r->N;
    
    //system("rm -v output/orbit*.txt");
    
    char* us = "_";
    char* txt = ".txt";
    char* intgrtr;
    
    char str[110] = {0};
    if(r->integrator == REB_INTEGRATOR_WHFAST){
        if(HYBRID_ON == 1)intgrtr = "HYBRID"; else intgrtr = "WHFAST";
        strcat(str, intgrtr);
        char* teq = "_t";
        strcat(str, teq);
        char dtstr[15];
        sprintf(dtstr, "%.0f", tmax);
        strcat(str, dtstr); //planet directory
        char* Neq = "_Np";
        char Nstr[15];
        sprintf(Nstr, "%d", N_planetesimals);
        strcat(str, Neq);
        strcat(str, Nstr);
        //char* epsln = "_Ep";
        //char eps[15];
        //sprintf(eps,"%.0e",epsilon);
        //strcat(str, epsln);
        //strcat(str, eps);
        
        strcat(str, us);
        char str3[15];
        sprintf(str3, "%.0f", r->ri_hybrid.switch_ratio);
        strcat(str, str3);
        strcat(str, us);
        char str2[15];
        sprintf(str2, "%.2f", drh);
        strcat(str, str2);
        
        //char strtime[10];
        //sprintf(strtime, "%d", hybrid_rint);
        
    } else{ //pure IAS15
        intgrtr = "IAS15";
        char* teq = "_t";
        strcat(str, intgrtr);
        strcat(str, teq);
        char dtstr[15];
        sprintf(dtstr, "%.0f", tmax);
        strcat(str, dtstr); //planet directory
    }
    
    strcat(legenddir, str);
    strcat(planetdir, str);
    strcat(planetdir, txt);
    
    //xyz_check
    strcat(xyz_check,str);
    char* err = "_xyz";
    strcat(xyz_check,err);
    strcat(xyz_check,txt);
    
    //CE print
    strcat(CEprint,str);
    char* CEs = "_CEs";
    strcat(CEprint,CEs);
    strcat(CEprint,txt);
    
    char* file = "Properties.txt";
    strcat(legenddir, us);
    strcat(legenddir, file);
    FILE *ff;
    ff=fopen(legenddir, "w");
    fprintf(ff,"General:\ndt, tmax, N_active, ri_hybrid.switch_ratio, dRHill, ias_epsilon, Seed, Integrator\n");
    fprintf(ff,"%f,%.1f,%d,%f,%f,%.2e,%d,%s \n\n",r->dt,tmax,N_active,r->ri_hybrid.switch_ratio,drh,epsilon,seed,intgrtr);
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
    
    char rmv2[100] = {0};
    strcat(rmv2, rm_v);
    strcat(rmv2, xyz_check);
    system(rmv2);
    
    char rmv3[100] = {0};
    strcat(rmv3, rm_v);
    strcat(rmv3, CEprint);
    system(rmv3);
    
}

double calc_dt(struct reb_simulation* r, double mp, double Ms, double a, double dRHill, double dt_prev){
    if(dRHill > r->ri_hybrid.switch_ratio){
        printf("\033[1mWarning!\033[0m dRhill !> N_RHill. Setting dRhill = N_Rhill/2 \n");
        dRHill = 0.5*r->ri_hybrid.switch_ratio;
    }
    double e_max = 0.3;  //max hypothesized eccentricity that the planet/esimals could have
    double Hill = a*(1 - e_max)*pow(mp/(3*Ms),1./3.);   //Hill radius of massive body
    double vmax = sqrt(r->G*(Ms + planetesimal_mass)*(1 + e_max)/(a*(1 - e_max)));   //perihelion speed of planetesimal
    double dt = dRHill*Hill/vmax;
    
    if(dt_prev < dt) dt = dt_prev;  //make sure I'm taking the smallest dt.
    
    return dt;
}

//Calc Hill Sphere (for speed in check_for_encounter)
void calc_Hill2(struct reb_simulation* r){
    struct reb_particle* restrict const particles = r->particles;
    struct reb_particle p0 = particles[0];
    for(int i=1;i<r->N;i++){
        struct reb_particle body = particles[i];
        double mp;
        if(i>=r->N_active) mp = planetesimal_mass; else mp = body.m;
        Hill2[i] = pow((mp/(p0.m*3.)), 2./3.);
    }
}

double calc_Etot(struct reb_simulation* a, double* K1, double* U1){
    double m1,m2;
    const int N = a->N;
    const int N_active = a->N_active;
    const double G = a->G;
    double L = 0, U = 0, K = 0;
    struct reb_particle* const particles = a->particles;
    for(int i=0;i<N;i++){
        struct reb_particle par = particles[i];
        if(i < N_active) m1 = par.m; else m1 = planetesimal_mass;
        //m1 = par.m;
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
                //m2 = par2.m;
                double ddx = dx - par2.x;
                double ddy = dy - par2.y;
                double ddz = dz - par2.z;
                U -= G*m1*m2/sqrt(ddx*ddx + ddy*ddy + ddz*ddz);
            }
        }
    }
    
    *K1 = K;
    *U1 = U;
    return K + U;
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

void planetesimal_forces_global(struct reb_simulation *r){
    const double G = r->G;
    const int N = r->N;
    const int N_active = r->N_active;
    struct reb_particle* const particles = r->particles;
    const double Gm1 = G*planetesimal_mass;
    
    for(int j=N_active;j<N;j++){//add planetesimal forces to massive bodies
        struct reb_particle p = particles[j];
        _Bool skip = 0;
        for(int k=0;k<N_encounters_previous && skip == 0;k++){
            if(p.id == previous_encounter_index[k]){
                skip++;
            }
        }
            
        if(skip == 0){//If CE don't calculate planetesimal forces in global
            for(int i=0;i<N_active;i++){
                struct reb_particle* body = &(particles[i]);
                const double dx = body->x - p.x;
                const double dy = body->y - p.y;
                const double dz = body->z - p.z;
                
                const double rijinv2 = 1.0/(dx*dx + dy*dy + dz*dz);
                const double ac = -Gm1*rijinv2*sqrt(rijinv2);
                
                body->ax += ac*dx;      //perturbation on planets due to planetesimals.
                body->ay += ac*dy;
                body->az += ac*dz;
            }
        }
    }
}

void planetesimal_forces_mini(struct reb_simulation *s){
    const double G = s->G;
    const int N = s->N;
    const int N_active = s->N_active;
    struct reb_particle* mini = s->particles;
    
    const double Gm1 = G*planetesimal_mass;
    for(int i=0;i<N_active;i++){
        struct reb_particle* body = &(mini[i]);
        for(int j=N_active;j<N;j++){//add planetesimal forces to massive bodies
            struct reb_particle p = mini[j];
            
            const double dx = body->x - p.x;
            const double dy = body->y - p.y;
            const double dz = body->z - p.z;
            
            const double rijinv2 = 1.0/(dx*dx + dy*dy + dz*dz);
            const double ac = -Gm1*rijinv2*sqrt(rijinv2);
            
            body->ax += ac*dx;      //perturbation on planets due to planetesimals.
            body->ay += ac*dy;
            body->az += ac*dz;
        }
    }
    
    //forces from global into mini
    struct reb_particle* const global = r->particles;
    const double timefac = (s->t - t_prev)/(r->t - t_prev);
    int rN_active = r->N_active;
    for(int i=rN_active;i<r->N;i++){    //planetesimals
        if(x_prev[i] != 0){             //find planetesimals with !=0 values, i.e. part of global but not mini
            const double ix = x_prev[i] + timefac*(global[i].x - x_prev[i]); //interpolated values
            const double iy = y_prev[i] + timefac*(global[i].y - y_prev[i]);
            const double iz = z_prev[i] + timefac*(global[i].z - z_prev[i]);
            for(int j=0;j<rN_active;j++){//massive bodies
                struct reb_particle* body = &(mini[j]);
                const double ddx = body->x - ix;
                const double ddy = body->y - iy;
                const double ddz = body->z - iz;
                
                const double rijinv2 = 1.0/(ddx*ddx + ddy*ddy + ddz*ddz);
                const double ac = -Gm1*rijinv2*sqrt(rijinv2);
                
                body->ax += ac*ddx;     //perturbation on planets due to planetesimals.
                body->ay += ac*ddy;
                body->az += ac*ddz;
            }
        }
    }
}

//initialize mini-simulation for close encounters
void ini_mini(struct reb_simulation* const r, struct reb_simulation* s, double ias_epsilon, int turn_planetesimal_forces_on, double ias_timestep){
    s->N_active = r->N_active;
    s->integrator = REB_INTEGRATOR_IAS15;
    if(turn_planetesimal_forces_on==1)s->additional_forces = planetesimal_forces_mini;
    s->exact_finish_time = 1;
    s->ri_ias15.epsilon = ias_epsilon;
    s->dt = ias_timestep;
    
    struct reb_particle* restrict const particles = r->particles;
    for(int k=0; k<s->N_active; k++){
        struct reb_particle p = particles[k];
        reb_add(s,p);
    }
}


//collect the id/array number of all planetesimals involved in a close encounter
void check_for_encounter(struct reb_simulation* const r, struct reb_simulation* const s, int* N_encounters, double* min_r, double* max_val, double ias_timestep){
    const int rN = r->N;
    const int rN_active = r->N_active;
    const int sN = s->N;
    const int sN_active = s->N_active;
    struct reb_particle* const global = r->particles;
    struct reb_particle* const mini = s->particles;
    struct reb_particle p0 = global[0];
    int num_encounters = 0;
    for (int i=1; i<rN_active; i++){
        struct reb_particle body = global[i];
        const double dxi = p0.x - body.x;
        const double dyi = p0.y - body.y;
        const double dzi = p0.z - body.z;
        const double r0i2 = dxi*dxi + dyi*dyi + dzi*dzi;
        const double rhi = r0i2*Hill2[i];
        //const double rhi = r0i2*pow((body.m/(p0.m*3.)), 2./3.); //can make this faster later
        
        for (int j=i+1; j<rN; j++){
            struct reb_particle pj = global[j];
            
            double HSR = r->ri_hybrid.switch_ratio;
            _Bool found_in_mini = 0;
            for(int k=sN_active; k<sN && found_in_mini == 0; k++){
                if(global[j].id == mini[k].id){
                    pj = mini[k];
                    HSR *= 1.07;    //HSR(mini) is a bit bigger so no constant enter/leave
                }
            }
            
            const double dxj = p0.x - pj.x;
            const double dyj = p0.y - pj.y;
            const double dzj = p0.z - pj.z;
            const double r0j2 = dxj*dxj + dyj*dyj + dzj*dzj;
            const double rhj = r0j2*Hill2[j];
            //const double rhj = r0j2*pow((pj.m/(p0.m*3.)), 2./3.); //can make this faster later
            
            const double dx = body.x - pj.x;
            const double dy = body.y - pj.y;
            const double dz = body.z - pj.z;
            const double rij2 = dx*dx + dy*dy + dz*dz;
            const double ratio = rij2/(rhi+rhj);    //(p-p distance/Hill radii)^2
            
            //if(fabs(pj.ax) + fabs(pj.ay) + fabs(pj.az) > 5000) printf("\nlarge acceleration at t=%f for par %d: ax=%f,ay=%f,az=%f",r->t,pj.id,pj.ax,pj.ay,pj.az);

            //int par_id = 120;
            //int CE_yes = 0;
            //if(ratio < HSR) CE_yes = 1;
            //if(pj.id==par_id && body.id==1){
                //printf("t=%f: CE par %d + pl %d,px=%f,py=%f,pz=%f,bx=%f,by=%f,bz=%f,r=%f,CE=%d\n",r->t,pj.id,body.id,pj.x,pj.y,pj.z,body.x,body.y,body.z,sqrt(rij2),CE_yes);
                //printf("t=%f: CE par %d + pl %d, r=%f,ratio=%f,rhi=%f,rhj=%f,dxj=%f,dyj=%f,dzj=%f,CE=%d\n",r->t,pj.id,body.id,sqrt(rij2),ratio,rhi,rhj,pj.x,pj.y,pj.z,CE_yes);
            //}
            
            if(ratio < HSR){
                num_encounters++;
                if(num_encounters == 1) encounter_index[0] = pj.id;
                else if(num_encounters > 1){//multiple close encounters
                    encounter_index = realloc(encounter_index,num_encounters*sizeof(int));
                    encounter_index[num_encounters - 1] = pj.id;
                }
                
                //Super close encounter
                if(rij2 < 1e-7){    //(2*radius of Neptune in AU)^2
                    fprintf(stderr,"\n\033[1mSuper Close Encounter at t=%f!\033[0m Particle %d and Planet %d collision should have happened, r=%f.\n",r->t,pj.id,body.id,sqrt(rij2));
                }
            }
            
            //calculate dt*(vrel/rmin)
            double vx = body.vx - pj.vx;
            double vy = body.vy - pj.vy;
            double vz = body.vz - pj.vz;
            double vrel = sqrt(vx*vx + vy*vy + vz*vz);
            double rr = sqrt(rij2);
            double val = r->dt*vrel/rr;
            if(rr < *min_r) *min_r = rr;
            if(val > *max_val) *max_val = val;
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

void add_or_subtract_particles(struct reb_simulation* r, struct reb_simulation* s, int N_encounters, int N_encounters_previous, int dN, char* CEprint){
    int N_active = s->N_active;
    struct reb_particle* mini = s->particles;
    struct reb_particle* global = r->particles;
    
    //check to add particles
    for(int i=0;i<N_encounters;i++){
        _Bool index_found = 0;
        int EI = encounter_index[i];
        if(EI >= N_active){//don't want to add/remove massive bodies. Already in mini
            for(int j=0;index_found == 0 && j<N_encounters_previous;j++){
                if(EI == previous_encounter_index[j]) index_found = 1;
            }
            if(index_found == 0){//couldn't find index
                struct reb_particle pt = global[EI];
                reb_add(s,pt);
                N_encounters_tot++;
                
                FILE *output;
                output = fopen(CEprint, "a");
                fprintf(output,"t=%f,%f particle %d added. dN == %d, N_close_encounters=%d\n",r->t,s->t,EI,dN,N_encounters);
                for(int i=0;i<N_encounters;i++)fprintf(output,"EI[%d]=%d,",i,encounter_index[i]);
                fprintf(output,"\n");
                for(int i=0;i<N_encounters_previous;i++)fprintf(output,"PEI[%d]=%d,",i,previous_encounter_index[i]);
                fprintf(output,"\n");
                fclose(output);
            }
        }
    }
    
    //check to remove particles
    for(int i=0;i<N_encounters_previous;i++){
        _Bool index_found = 0;
        int PEI = previous_encounter_index[i];
        if(PEI >= N_active){//don't want to add/remove massive bodies. Already in mini
            for(int j=0;index_found == 0 && j<N_encounters;j++){
                if(PEI == encounter_index[j]) index_found = 1;
            }
            if(index_found == 0){//couldn't find index, remove particle
                int removed_particle = 0;
                for(int k=0;removed_particle==0 && k<N_encounters_previous;k++){
                    if(mini[k+N_active].id == PEI){
                        removed_particle = reb_remove(s,k+N_active,1);    //remove particle
                        
                        FILE *output;
                        output = fopen(CEprint, "a");
                        fprintf(output,"t=%f,%f particle %d leaving. dN == %d, N_close_encounters=%d.\n",r->t,s->t,PEI,dN,N_encounters);
                        for(int i=0;i<N_encounters;i++)fprintf(output,"EI[%d]=%d,",i,encounter_index[i]);
                        fprintf(output,"\n");
                        for(int i=0;i<N_encounters_previous;i++)fprintf(output,"PEI[%d]=%d,",i,previous_encounter_index[i]);
                        fprintf(output,"\n");
                        fclose(output);
                    }
                }
            }
        }
    }
}


void update_previous_global_positions(struct reb_simulation* r, int N_encounters){
    struct reb_particle* global = r->particles;
    t_prev = r->t;
    for(int i=r->N_active;i<r->N;i++){
        int ID = global[i].id;
        _Bool found_particle = 0;
        for(int j=0;j<N_encounters && found_particle == 0;j++){
            if(ID == encounter_index[j]){
                found_particle = 1;
                x_prev[i] = 0.;         //reset planetesimals involved in mini to 0
                y_prev[i] = 0.;
                z_prev[i] = 0.;
            }
        }
        if(found_particle == 0){        //planetesimal not involved in mini, update position for later interp.
            x_prev[i] = global[i].x;
            y_prev[i] = global[i].y;
            z_prev[i] = global[i].z;
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

time_t clock_start(){
    char buf[64];
    time_t t_ini = time(NULL);
    struct tm *tmp = gmtime(&t_ini);
    strftime(buf, sizeof(buf), "%j:%H:%M:%S\n", tmp);
    printf("start time (GMT): %s\n",buf);
    
    return t_ini;
}

void clock_finish(clock_t t_ini, int N_encounters, char* legenddir){
    char buf[64];
    time_t t_fini = time(NULL);
    struct tm *tmp = gmtime(&t_fini);
    strftime(buf, sizeof(buf), "%j:%H:%M:%S\n", tmp);
    printf("\nfinish time (GMT): %s\n",buf);
    
    double time = t_fini - t_ini;
    
    FILE *ff;
    ff=fopen(legenddir, "a");
    fprintf(ff,"Elapsed simulation time is %.2f s, with %d close encounters.\n",time,N_encounters);
    printf("\nSimulation complete. Elapsed simulation time is %.2f s, with %d close encounters.\n\n",time,N_encounters);
}

/*
void clock_finish(clock_t timer, int N_encounters, char* legenddir){
    timer = clock() - timer;
    FILE *ff;
    ff=fopen(legenddir, "a");
    double result = ((float)timer)/CLOCKS_PER_SEC;
    fprintf(ff,"Elapsed simulation time is %f s, with %d close encounters.\n",result,N_encounters);
    printf("\n\nSimulation complete. Elapsed simulation time is %f s, with %d close encounters.\n\n",result,N_encounters);
}*/

void global_free(){
    free(encounter_index);
    free(previous_encounter_index);
    free(Hill2);
    free(x_prev);
    free(y_prev);
    free(z_prev);
}
