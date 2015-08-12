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
