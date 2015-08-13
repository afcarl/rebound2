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

double calc_dt(struct reb_simulation* r, double mp, double Ms, double a, double e_max, double N_Rhill, double dRHill){
    if(dRHill > N_Rhill){
        printf("\033[1mWarning!\033[0m dRhill !> N_RHill. Setting dRhill = N_Rhill/2 \n");
        dRHill = 0.5*N_Rhill;
    }
    double Hill = a*(1 - e_max)*pow(mp/(3*Ms),1./3.);
    double r2_E = N_Rhill*N_Rhill*Hill*Hill;
    r->ri_hybrid.switch_ratio = Ms*r2_E/(mp*a*a*1.21);
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
