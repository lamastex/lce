/* 
 * Copyright (C) 2005--2009 Raazesh Sainudiin and Kevin Thornton
 * 
 * This file is part of lce, a C++ class library for lumped coalescent experiments.
 *
 * lce is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */
#ifndef __SFSPARTICLES_HPP__ 
#define __SFSPARTICLES_HPP__ 

#include<ebc_sfsinteg.hpp>


void iterate_sfs_particles_Repel(gsl_rng* rgsl, sfs_particlesV& SFS_Particles, unsigned ParticleWalkLength, 
		unsigned NumberOfParticles, double Temperature,
		std::vector<OneMove> & MovesVector, unsigned NumberOfMoves, const std::valarray<double>& Pis,
		const params& p, const SfsInfo& SFSInfo);


#endif
