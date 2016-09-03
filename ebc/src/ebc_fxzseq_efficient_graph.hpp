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
/*! \file ebc_fxzseq_efficient_graph.hpp
    \brief prototype of BGL methods for unvintaged (sized) n-coalescent and associated shape statistics with vertices as iterators to fstate and zstate sets, ie graph over product space
*/
#ifndef __FXZSEQEFFICIENTGRAPH_HPP__
#define __FXZSEQEFFICIENTGRAPH_HPP__

#include <ebc_sfstypes.hpp>
#include <ebc_graphtypes.hpp>

#include <gsl/gsl_rng.h>

//libsequence headers

bool generateControlledFSequenceWithGraph(f_x_z_efficient_graph& graph,
                                ExtSetItPairVertexMap& extStatesPairVertexMap,
                                state_typeSet& extFStateSet,
                                state_typeSet& extZStateSet,
                                gsl_rng* rgsl,
                                StateSequence& fs,
                                const state_type& z_state,
                                double& fsProposalProb,
                                double& fsKCProb,
                                const size_t nsam,
                                const bool efficientFind = true);


bool generateFSequenceWithGraph(f_x_z_efficient_graph& graph,
                                ExtSetItPairVertexMap& extStatesPairVertexMap,
                                state_typeSet& extFStateSet,
                                state_typeSet& extZStateSet,
                                gsl_rng* rgsl,
                                StateSequence& fs,
                                double& fsProposalProb,
                                double& fsKCProb,
                                const size_t nsam,
                                const bool efficientFind = true);


bool generateCompleteControlledGraph(f_x_z_efficient_graph& graph,
                                    ExtSetItPairVertexMap& extStatesPairVertexMap,
                                    state_typeSet& extFStateSet,
                                    state_typeSet& extZStateSet,
                                    const size_t nsam,
                                    const state_type& z_state);

bool generateCompleteGraph(f_x_z_efficient_graph& graph,
                            ExtSetItPairVertexMap& extStatesPairVertexMap,
                            state_typeSet& extFStateSet,
                            state_typeSet& extZStateSet,
                            const size_t nsam);

fxz_efficient_graph_vertex_descriptor firstVertex(f_x_z_efficient_graph& graph,
                                ExtSetItPairVertexMap& extStatesPairVertexMap,
                                state_typeSet& extFStateSet,
                                state_typeSet& extZStateSet,
                                const size_t nsam,
                                const state_type& z_state);

bool nextFZSequenceState(f_x_z_efficient_graph& graph,
                        ExtSetItPairVertexMap& extStatesPairVertexMap,
                        state_typeSet& extFStateSet,
                        state_typeSet& extZStateSet,
                        gsl_rng* rgsl,
                        StateSequence& fs,
                        double& fsProposalProb,
                        double& fsKCProb,
                        const size_t nsam,
                        const fxz_efficient_graph_vertex_descriptor v,
                        const size_t epoch,
                        bool pathTravelledBefore);

bool allNextFZSequenceStates(f_x_z_efficient_graph& graph,
                            ExtSetItPairVertexMap& extStatesPairVertexMap,
                            state_typeSet& extFStateSet,
                            state_typeSet& extZStateSet,
                            const size_t nsam,
                            const fxz_efficient_graph_vertex_descriptor v,
                            const size_t epoch);


bool expandGraphKingmanForward(f_x_z_efficient_graph& graph,
                                ExtSetItPairVertexMap& extStatesPairVertexMap,
                                state_typeSet& extFStateSet,
                                state_typeSet& extZStateSet,
                                const size_t nsam,
                                const fxz_efficient_graph_vertex_descriptor v,
                                const size_t epoch);


FxzEfficientVertexChoiceTuple chooseNextVertex(f_x_z_efficient_graph& graph,
                                                gsl_rng* rgsl,
                                                const fxz_efficient_graph_vertex_descriptor v,
                                                bool pathTravelledBefore);

FxzEfficientGraphVertexBoolPair vertexAddOrFind(f_x_z_efficient_graph& graph,
                                                ExtSetItPairVertexMap& extStatesPairVertexMap,
                                                const state_typeSetIt& f_stateSetIt,
                                                const state_typeSetIt& z_stateSetIt);





FxzEfficientGraphEdgeBoolPair edgeAddOrFind(f_x_z_efficient_graph& graph,
                                            const fxz_efficient_graph_vertex_descriptor v1,
                                            const fxz_efficient_graph_vertex_descriptor v2,
                                            const double edgeRWprob,
                                            const double edgeKCprob);

//void addStateToF_Sequence(StateSequence& fs, const state_type& f_state,
       // const size_t epoch, const size_t nsam);

//void addStateToZ_Sequence(state_type& zs, const state_type& z_state,
       // const size_t epoch, const size_t nsam);

//string statesToString(const state_type& f_state, const state_type& z_state);

//bool checkFinalState(const state_type& f_state, const size_t nsam);

void print_graph_desc(f_x_z_efficient_graph& graph);


/*
template <typename T>
void stdoutArray(const T& t)
{
    std::cout << "(";
    for (int i = 0; i < t.size(); i++) {
        std::cout << t[i];
        if (i < t.size()-1) std::cout << "\t";
    }
    std::cout << ")";
}
*/

#endif
