/**
 * @file 	gravity.c
 * @brief 	Direct gravity calculation, O(N^2).
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 *
 * @details 	This is the crudest implementation of an N-body code
 * which sums up every pair of particles. It is only useful very small 
 * particle numbers (N<~100) as it scales as O(N^2). Note that the MPI
 * implementation is not well tested and only works for very specific
 * problems. This should be resolved in the future. 
 *
 * 
 * @section LICENSE
 * Copyright (c) 2011 Hanno Rein, Shangfei Liu
 *
 * This file is part of rebound.
 *
 * rebound is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * rebound is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with rebound.  If not, see <http://www.gnu.org/licenses/>.
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include "particle.h"
#include "rebound.h"
#include "tree.h"
#include "boundary.h"

#ifdef MPI
#include "communication_mpi.h"
#endif

/**
  * The function loops over all trees to call calculate_forces_for_particle_from_cell() tree to calculate forces for each particle.
  *
  * @param pt Index of the particle the force is calculated for.
  * @param gb Ghostbox plus position of the particle (precalculated). 
  */
void gravity_calculate_acceleration_for_particle(struct Rebound* const r, const int pt, const struct Ghostbox gb);


/**
 * Main Gravity Routine
 */
void gravity_calculate_acceleration(struct Rebound* r){
	struct Particle* const particles = r->particles;
	const int N = r->N;
	const int N_megno = r->N_megno;
	const int N_active = r->N_active;
	const double G = r->G;
	const double softening2 = r->softening*r->softening;
	const unsigned int _gravity_ignore_10 = r->gravity_ignore_10;
	const int _N_start  = (r->integrator==RB_IT_WH?1:0);
	const int _N_active = ((N_active==-1)?N:N_active)- N_megno;
	const int _N_real   = N - N_megno;
	switch (r->gravity){
		case RB_GT_NONE: // Do nothing.
		break;
		case RB_GT_BASIC:
		{
			const int nghostx = r->nghostx;
			const int nghosty = r->nghosty;
			const int nghostz = r->nghostz;
#pragma omp parallel for schedule(guided)
			for (int i=0; i<N; i++){
				particles[i].ax = 0; 
				particles[i].ay = 0; 
				particles[i].az = 0; 
			}
			// Summing over all Ghost Boxes
			for (int gbx=-nghostx; gbx<=nghostx; gbx++){
			for (int gby=-nghosty; gby<=nghosty; gby++){
			for (int gbz=-nghostz; gbz<=nghostz; gbz++){
				struct Ghostbox gb = boundary_get_ghostbox(r, gbx,gby,gbz);
				// Summing over all particle pairs
#pragma omp parallel for schedule(guided)
				for (int i=_N_start; i<_N_real; i++){
				for (int j=_N_start; j<_N_active; j++){
					if (_gravity_ignore_10 && j==1 && i==0 ) continue;
					if (i==j) continue;
					const double dx = (gb.shiftx+particles[i].x) - particles[j].x;
					const double dy = (gb.shifty+particles[i].y) - particles[j].y;
					const double dz = (gb.shiftz+particles[i].z) - particles[j].z;
					const double _r = sqrt(dx*dx + dy*dy + dz*dz + softening2);
					const double prefact = -G/(_r*_r*_r)*particles[j].m;
					
					particles[i].ax    += prefact*dx;
					particles[i].ay    += prefact*dy;
					particles[i].az    += prefact*dz;
				}
				}
			}
			}
			}
		}
		break;
		case RB_GT_COMPENSATED:
		{
			if (r->N_cs<_N_real){
				r->cs = realloc(r->cs,_N_real*sizeof(struct rb_vec3d));
				r->N_cs = _N_real;
			}
			struct rb_vec3d* const cs = r->cs;
#pragma omp parallel for schedule(guided)
			for (int i=0; i<_N_real; i++){
				particles[i].ax = 0.; 
				particles[i].ay = 0.; 
				particles[i].az = 0.; 
				cs[i].x = 0.;
				cs[i].y = 0.;
				cs[i].z = 0.;
			}
			// Summing over all massive particle pairs
#pragma omp parallel for schedule(guided)
			for (int i=_N_start; i<_N_active; i++){
			for (int j=i+1; j<_N_active; j++){
				if (_gravity_ignore_10 && j==1 && i==0 ) continue;
				const double dx = particles[i].x - particles[j].x;
				const double dy = particles[i].y - particles[j].y;
				const double dz = particles[i].z - particles[j].z;
				const double r2 = dx*dx + dy*dy + dz*dz + softening2;
				const double r = sqrt(r2);
				const double prefact  = -G/(r2*r);
				const double prefacti = prefact*particles[i].m;
				const double prefactj = prefact*particles[j].m;
				
				{
				double ax = particles[i].ax;
				cs[i].x  +=	prefactj*dx; 
				particles[i].ax    = ax + cs[i].x;
				cs[i].x  += ax - particles[i].ax; 
				
				double ay = particles[i].ay;
				cs[i].y  +=	prefactj*dy; 
				particles[i].ay    = ay + cs[i].y;
				cs[i].y  += ay - particles[i].ay; 
				
				double az = particles[i].az;
				cs[i].z  +=	prefactj*dz; 
				particles[i].az    = az + cs[i].z;
				cs[i].z  += az - particles[i].az; 
				}
				
				{
				double ax = particles[j].ax;
				cs[j].x  -=	prefacti*dx; 
				particles[j].ax    = ax + cs[j].x;
				cs[j].x  += ax - particles[j].ax; 
				
				double ay = particles[j].ay;
				cs[j].y  -=	prefacti*dy; 
				particles[j].ay    = ay + cs[j].y;
				cs[j].y  += ay - particles[j].ay; 
				
				double az = particles[j].az;
				cs[j].z  -=	prefacti*dz; 
				particles[j].az    = az + cs[j].z;
				cs[j].z  += az - particles[j].az; 
				}
				
			}
			}
			// Testparticles
#pragma omp parallel for schedule(guided)
			for (int i=_N_active; i<_N_real; i++){
			for (int j=_N_start; j<_N_active; j++){
				if (_gravity_ignore_10 && i==1 && j==0 ) continue;
				const double dx = particles[i].x - particles[j].x;
				const double dy = particles[i].y - particles[j].y;
				const double dz = particles[i].z - particles[j].z;
				const double r2 = dx*dx + dy*dy + dz*dz + softening2;
				const double r = sqrt(r2);
				const double prefact = -G/(r2*r)*particles[j].m;
				
				double ax = particles[i].ax;
				cs[i].x  +=	prefact*dx; 
				particles[i].ax    = ax + cs[i].x;
				cs[i].x  += ax - particles[i].ax; 
				
				double ay = particles[i].ay;
				cs[i].y  +=	prefact*dy; 
				particles[i].ay    = ay + cs[i].y;
				cs[i].y  += ay - particles[i].ay; 
				
				double az = particles[i].az;
				cs[i].z  +=	prefact*dz; 
				particles[i].az    = az + cs[i].z;
				cs[i].z  += az - particles[i].az; 
			}
			}
		}
		break;
		case RB_GT_TREE:
		{
#pragma omp parallel for schedule(guided)
			for (int i=0; i<N; i++){
				particles[i].ax = 0; 
				particles[i].ay = 0; 
				particles[i].az = 0; 
			}
			// Summing over all Ghost Boxes
			for (int gbx=-r->nghostx; gbx<=r->nghostx; gbx++){
			for (int gby=-r->nghosty; gby<=r->nghosty; gby++){
			for (int gbz=-r->nghostz; gbz<=r->nghostz; gbz++){
				// Summing over all particle pairs
#pragma omp parallel for schedule(guided)
				for (int i=0; i<N; i++){
					struct Ghostbox gb = boundary_get_ghostbox(r, gbx,gby,gbz);
					// Precalculated shifted position
					gb.shiftx += particles[i].x;
					gb.shifty += particles[i].y;
					gb.shiftz += particles[i].z;
					gravity_calculate_acceleration_for_particle(r, i, gb);
				}
			}
			}
			}
		}
		break;
		default:
			fprintf(stderr,"\n\033[1mError!\033[0m Gravity calculation not yet implemented.\n");
			exit(EXIT_FAILURE);
			break;
		break;
	}
#ifdef MPI
	// Distribute active particles from root to all other nodes.
	// This assures that round-off errors do not accumulate and 
	// the copies of active particles do not diverge. 
	MPI_Bcast(particles, N_active, mpi_particle, 0, MPI_COMM_WORLD); 
#endif

}

void gravity_calculate_variational_acceleration(struct Rebound* r){
	struct Particle* const particles = r->particles;
	const double G = r->G;
	const unsigned int _gravity_ignore_10 = r->gravity_ignore_10;
	const int N = r->N;
	const int N_megno = r->N_megno;
	const int _N_real   = N - N_megno;
#pragma omp parallel for schedule(guided)
	for (int i=_N_real; i<N; i++){
		particles[i].ax = 0; 
		particles[i].ay = 0; 
		particles[i].az = 0; 
	}
	switch (r->gravity){
		case RB_GT_NONE: // Do nothing.
		break;
		case RB_GT_BASIC:
		case RB_GT_COMPENSATED:
#pragma omp parallel for schedule(guided)
			for (int i=_N_real; i<N; i++){
			for (int j=i+1; j<N; j++){
				if (_gravity_ignore_10 && ((i==_N_real+1 && j==_N_real) || (j==_N_real+1 && i==_N_real)) ) continue;
				const double dx = particles[i-N/2].x - particles[j-N/2].x;
				const double dy = particles[i-N/2].y - particles[j-N/2].y;
				const double dz = particles[i-N/2].z - particles[j-N/2].z;
				const double r2 = dx*dx + dy*dy + dz*dz;
				const double _r  = sqrt(r2);
				const double r3inv = 1./(r2*_r);
				const double r5inv = 3.*r3inv/r2;
				const double ddx = particles[i].x - particles[j].x;
				const double ddy = particles[i].y - particles[j].y;
				const double ddz = particles[i].z - particles[j].z;
				const double Gmi = G * particles[i].m;
				const double Gmj = G * particles[j].m;
				
				const double dax =   ddx * ( dx*dx*r5inv - r3inv )
						   + ddy * ( dx*dy*r5inv )
						   + ddz * ( dx*dz*r5inv );
				const double day =   ddx * ( dy*dx*r5inv )
						   + ddy * ( dy*dy*r5inv - r3inv )
						   + ddz * ( dy*dz*r5inv );
				const double daz =   ddx * ( dz*dx*r5inv )
						   + ddy * ( dz*dy*r5inv )
						   + ddz * ( dz*dz*r5inv - r3inv );
				
				particles[i].ax += Gmj * dax;
				particles[i].ay += Gmj * day;
				particles[i].az += Gmj * daz;
				
				particles[j].ax -= Gmi * dax;
				particles[j].ay -= Gmi * day;
				particles[j].az -= Gmi * daz;
			}
			}
			break;
		default:
			fprintf(stderr,"\n\033[1mError!\033[0m Variational equations not implemented.\n");
			exit(EXIT_FAILURE);
			break;
	}

}


// Helper routines for RB_GT_TREE


/**
  * The function calls itself recursively using cell breaking criterion to check whether it can use center of mass (and mass quadrupole tensor) to calculate forces.
  * Calculate the acceleration for a particle from a given cell and all its daughter cells.
  *
  * @param pt Index of the particle the force is calculated for.
  * @param node Pointer to the cell the force is calculated from.
  * @param gb Ghostbox plus position of the particle (precalculated). 
  */
void gravity_calculate_acceleration_for_particle_from_cell(struct Rebound* const r, const int pt, const struct cell *node, const struct Ghostbox gb);

void gravity_calculate_acceleration_for_particle(struct Rebound* const r, const int pt, const struct Ghostbox gb) {
	for(int i=0;i<r->root_n;i++){
		struct cell* node = r->tree_root[i];
		if (node!=NULL){
			gravity_calculate_acceleration_for_particle_from_cell(r, pt, node, gb);
		}
	}
}

void gravity_calculate_acceleration_for_particle_from_cell(struct Rebound* r, const int pt, const struct cell *node, const struct Ghostbox gb) {
	const double G = r->G;
	const double softening2 = r->softening*r->softening;
	struct Particle* const particles = r->particles;
	const double dx = gb.shiftx - node->mx;
	const double dy = gb.shifty - node->my;
	const double dz = gb.shiftz - node->mz;
	const double r2 = dx*dx + dy*dy + dz*dz;
	if ( node->pt < 0 ) { // Not a leaf
		if ( node->w*node->w > r->opening_angle2*r2 ){
			for (int o=0; o<8; o++) {
				if (node->oct[o] != NULL) {
					gravity_calculate_acceleration_for_particle_from_cell(r, pt, node->oct[o], gb);
				}
			}
		} else {
			double _r = sqrt(r2 + softening2);
			double prefact = -G/(_r*_r*_r)*node->m;
#ifdef QUADRUPOLE
			double qprefact = G/(_r*_r*_r*_r*_r);
			particles[pt].ax += qprefact*(dx*node->mxx + dy*node->mxy + dz*node->mxz); 
			particles[pt].ay += qprefact*(dx*node->mxy + dy*node->myy + dz*node->myz); 
			particles[pt].az += qprefact*(dx*node->mxz + dy*node->myz + dz*node->mzz); 
			double mrr 	= dx*dx*node->mxx 	+ dy*dy*node->myy 	+ dz*dz*node->mzz
					+ 2.*dx*dy*node->mxy 	+ 2.*dx*dz*node->mxz 	+ 2.*dy*dz*node->myz; 
			qprefact *= -5.0/(2.0*_r*_r)*mrr;
			particles[pt].ax += (qprefact + prefact) * dx; 
			particles[pt].ay += (qprefact + prefact) * dy; 
			particles[pt].az += (qprefact + prefact) * dz; 
#else
			particles[pt].ax += prefact*dx; 
			particles[pt].ay += prefact*dy; 
			particles[pt].az += prefact*dz; 
#endif
		}
	} else { // It's a leaf node
		if (node->pt == pt) return;
		double _r = sqrt(r2 + softening2);
		double prefact = -G/(_r*_r*_r)*node->m;
		particles[pt].ax += prefact*dx; 
		particles[pt].ay += prefact*dy; 
		particles[pt].az += prefact*dz; 
	}
}

