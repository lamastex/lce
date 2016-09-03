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
/*! \file ebc_fseq_efficient_graph.hpp
    \brief prototype of BGL methods for unvintaged (sized) n-coalescent and associated shape statistics with vertices as iterators to an fstate set, ie graph over f-space only
*/
#ifndef __FSEQEFFICIENTGRAPH_HPP__
#define __FSEQEFFICIENTGRAPH_HPP__

#include <ebc_sfstypes.hpp>
#include <ebc_graphtypes.hpp>

#include <gsl/gsl_rng.h>

#include <boost/graph/breadth_first_search.hpp>


//libsequence headers

bool generateControlledFSequenceWithGraph(f_efficient_graph& graph,
                                        ExtSetItVertexMap& extStateVertexMap,
                                        state_typeSet& extFStateSet,
                                        state_typeSet& extZStateSet,
                                        gsl_rng* rgsl,
                                        StateSequence& fs,
                                        const state_type& z_state,
                                        double& fsProposalProb,
                                        double& fsKCProb,
                                        const size_t nsam,
                                        const bool efficientFind = false);


bool generateFSequenceWithGraph(f_efficient_graph& graph,
                                ExtSetItVertexMap& extStateVertexMap,
                                state_typeSet& extFStateSet,
                                state_typeSet& extZStateSet,
                                gsl_rng* rgsl,
                                StateSequence& fs,
                                double& fsKCProb,
                                double& fsProposalProb,
                                const size_t nsam,
                                const bool efficientFind = false);


bool generateCompleteControlledGraph(f_efficient_graph& graph,
                                    ExtSetItVertexMap& extStateVertexMap,
                                    state_typeSet& extFStateSet,
                                    state_typeSet& extZStateSet,
                                    const size_t nsam,
                                    const state_type& z_state);

bool generateCompleteGraph(f_efficient_graph& graph,
                            ExtSetItVertexMap& extStateVertexMap,
                            state_typeSet& extFStateSet,
                            state_typeSet& extZStateSet,
                            const size_t nsam);

f_efficient_graph_vertex_descriptor firstVertex(f_efficient_graph& graph,
                                                    ExtSetItVertexMap& extStateVertexMap,
                                                    state_typeSet& extFStateSet,
                                                    state_typeSet& extZStateSet,
                                                    const size_t nsam,
                                                    const state_type& z_state);

bool nextStateSequenceState(f_efficient_graph& graph,
                        ExtSetItVertexMap& extStateVertexMap,
                        state_typeSet& extFStateSet,
                        state_typeSet& extZStateSet,
                        gsl_rng* rgsl,
                        StateSequence& fs,
                        double& fsProposalProb,
                        double& fsKCProb,
                        const size_t nsam,
                        const f_efficient_graph_vertex_descriptor v,
                        const state_type& z_state,
                        const size_t epoch);


bool allNextStateSequenceStates(f_efficient_graph& graph,
                            ExtSetItVertexMap& extStateVertexMap,
                            state_typeSet& extFStateSet,
                            state_typeSet& extZStateSet,
                            const size_t nsam,
                            const f_efficient_graph_vertex_descriptor v,
                            const state_type& z_state,
                            const size_t epoch);


bool expandGraphKingmanForward(f_efficient_graph& graph,
                                ExtSetItVertexMap& extStateVertexMap,
                                state_typeSet& extFStateSet,
                                state_typeSet& extZStateSet,
                                VecSetIt& new_f_statesSetIt,
                                VecSetIt& new_z_statesSetIt,
                                const size_t nsam,
                                const f_efficient_graph_vertex_descriptor v,
                                const state_type& z_state,
                                const size_t epoch);


F_EfficientVertexChoiceTuple chooseNextVertex(f_efficient_graph& graph,
                                                gsl_rng* rgsl,
                                                const size_t nsam,
                                                const f_efficient_graph_vertex_descriptor v,
                                                const state_type& z_state);


void findAllowableNextVertices(f_efficient_graph& graph,
                                const size_t nsam,
                                const f_efficient_graph_vertex_descriptor v,
                                const state_type& z_state,
                                vector<f_efficient_graph_vertex_descriptor>& vecVertices,
                                vector<double>& vecProbs,
                                double& totalAllowedProb);


bool findNextZState(const size_t nsam,
                    const state_type& f_state,
                    const state_type& new_f_state,
                    const state_type& z_state,
                    state_type& new_z_state);


F_EfficientGraphVertexBoolPair vertexAddOrFind(f_efficient_graph& graph,
                                ExtSetItVertexMap& extStateVertexMap,
                                const state_typeSetIt& f_stateSetIt);


F_EfficientGraphEdgeBoolPair edgeAddOrFind(f_efficient_graph& graph,
                                    const f_efficient_graph_vertex_descriptor v1,
                                    const f_efficient_graph_vertex_descriptor v2,
                                    const double edgeKCprob);

//void addStateToF_Sequence(StateSequence& fs, const state_type& f_state,
       // const size_t epoch, const size_t nsam);

//void addStateToZ_Sequence(state_type& zs, const state_type& z_state,
       // const size_t epoch, const size_t nsam);

//string statesToString(const state_type& f_state, const state_type& z_state);

//bool checkFinalState(const state_type& f_state, const size_t nsam);

// print description of the graph
void print_graph_desc(f_efficient_graph& graph);

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

// breadth-first searching
class bfs_accumulate_count_visitor:public default_bfs_visitor {

public:
    bfs_accumulate_count_visitor(size_t* cp) : countPtr(cp), good(true) {};

    void discover_vertex(f_efficient_graph_vertex_descriptor u,
        const f_efficient_graph& g);

    size_t* countPtr; // countPtr is public
    std::map<f_efficient_graph_vertex_descriptor, size_t> accMap;
    bool good;
};

size_t  calcFsequences(f_efficient_graph& graph);



void useGraphToCountFsequences (const size_t max_n,
                                const size_t min_n,
                                bool doXstars,
                                const string& filename);

#endif
