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
/*! \file ebc_fxzseq_named_graph.hpp
    \brief prototype of BGL methods for unvintaged (sized) n-coalescent and associated shape statistics using a named graph over the f x z product space which is inefficient on memory but good for images for small nsam
*/
#ifndef __FXZSEQNAMEDGRAPH_HPP__
#define __FXZSEQNAMEDGRAPH_HPP__

#include <ebc_sfstypes.hpp>
#include <ebc_graphtypes.hpp>

#include <gsl/gsl_rng.h>

//libsequence headers

bool generateControlledFSequenceWithGraph(f_x_z_named_graph& graph,
                                        ExtSetItPairVertexMap& extStatesPairVertexMap,
                                        state_typeSet& extFStateSet,
                                        state_typeSet& extZStateSet,
                                        gsl_rng* rgsl,
                                        StateSequence& fs,
                                        const state_type& z_state,
                                        double& fsProposalProb,
                                        double& fsKCProb,
                                        const size_t nsam,
                                        const bool efficientFind = true,
                                        const string addColour = "black");


bool generateFSequenceWithGraph(f_x_z_named_graph& graph,
                                ExtSetItPairVertexMap& extStatesPairVertexMap,
                                state_typeSet& extFStateSet,
                                state_typeSet& extZStateSet,
                                gsl_rng* rgsl,
                                StateSequence& fs,
                                double& fsProposalProb,
                                double& fsKCProb,
                                const size_t nsam,
                                const bool efficientFind = true,
                                const string addColour = "black");

bool generateCompleteControlledGraph(f_x_z_named_graph& graph,
                                    ExtSetItPairVertexMap& extStatesPairVertexMap,
                                    state_typeSet& extFStateSet,
                                    state_typeSet& extZStateSet,
                                    const size_t nsam,
                                    const state_type& z_state);

bool generateCompleteGraph(f_x_z_named_graph& graph,
                            ExtSetItPairVertexMap& extStatesPairVertexMap,
                            state_typeSet& extFStateSet,
                            state_typeSet& extZStateSet,
                            const size_t nsam);

fxz_named_graph_vertex_descriptor firstVertex(f_x_z_named_graph& graph,
                                ExtSetItPairVertexMap& extStatesPairVertexMap,
                                state_typeSet& extFStateSet,
                                state_typeSet& extZStateSet,
                                const size_t nsam, const state_type& z_state);

bool nextFZSequenceState(f_x_z_named_graph& graph,
                        ExtSetItPairVertexMap& extStatesPairVertexMap,
                        state_typeSet& extFStateSet,
                        state_typeSet& extZStateSet,
                        gsl_rng* rgsl,
                        StateSequence& fs,
                        StateSequence& zs,
                        double& fsProposalProb,
                        double& fsKCProb,
                        const size_t nsam,
                        const fxz_efficient_graph_vertex_descriptor v,
                        size_t epoch,
                        bool pathTravelledBefore,
                        const string addColour = "black");

bool allNextFZSequenceStates(f_x_z_named_graph& graph,
                            ExtSetItPairVertexMap& extStatesPairVertexMap,
                            state_typeSet& extFStateSet,
                            state_typeSet& extZStateSet,
                            const size_t nsam,
                            const fxz_efficient_graph_vertex_descriptor v,
                            const size_t epoch);


bool expandGraphKingmanForward(f_x_z_named_graph& graph,
                                ExtSetItPairVertexMap& extStatesPairVertexMap,
                                state_typeSet& extFStateSet,
                                state_typeSet& extZStateSet,
                                const size_t nsam,
                                const fxz_efficient_graph_vertex_descriptor v,
                                const size_t epoch);


FxzNamedVertexChoiceTuple chooseNextVertex(f_x_z_named_graph& graph,
                                            gsl_rng* rgsl,
                                            const fxz_efficient_graph_vertex_descriptor v,
                                            bool pathTravelledBefore);

FxzNamedGraphVertexBoolPair vertexAddOrFind(f_x_z_named_graph& graph,
                                ExtSetItPairVertexMap& extStatesPairVertexMap,
                                const state_typeSetIt& f_stateSetIt,
                                const state_typeSetIt& z_stateSetIt);


FxzEfficientGraphEdgeBoolPair edgeAddOrFind(f_x_z_named_graph& graph,
                                    const fxz_efficient_graph_vertex_descriptor v1,
                                    const fxz_efficient_graph_vertex_descriptor v2,
                                    const double edgeRWprob,
                                    const double edgeKCprob);

void print_graph_desc(f_x_z_named_graph& graph);

void makeDotGraph(f_x_z_named_graph& graph, string filename = "myGraph.dot");


void dijShortestPath(f_x_z_named_graph& graph,
    ExtSetItPairVertexMap& extStatesPairVertexMap,
    state_typeSet& extFStateSet,
    state_typeSet& extZStateSet,
    const state_type& z_state,
    const size_t nsam,
    const string addColour = "black");


void dijLongestPath(f_x_z_named_graph& graph,
    ExtSetItPairVertexMap& extStatesPairVertexMap,
    state_typeSet& extFStateSet,
    state_typeSet& extZStateSet,
    const state_type& z_state,
    const size_t nsam,
    const string addColour = "black");


void mapAverageOuts(f_x_z_named_graph& graph,
    ExtSetItPairVertexMap& extStatesPairVertexMap,
    state_typeSet& extFStateSet,
    state_typeSet& extZStateSet,
    const state_type& z_state,
    const size_t nsam);


pair<fxz_named_graph_vertex_descriptor, fxz_named_graph_vertex_descriptor>
    findSourceAndSink(f_x_z_named_graph& graph,
    ExtSetItPairVertexMap& extStatesPairVertexMap,
    state_typeSet& extFStateSet,
    state_typeSet& extZStateSet,
    const state_type& z_state,
    const size_t nsam);

void colourPath(f_x_z_named_graph& graph,
                FxzEfficientVecVertices& path,
                const string addColour);

void printPath(f_x_z_named_graph& graph,
                FxzEfficientVecVertices& path);


#endif
