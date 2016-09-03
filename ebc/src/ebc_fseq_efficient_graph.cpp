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
/*! \file ebc_fseq_efficient_graph.cpp
    \brief BGL methods for unvintaged (sized) n-coalescent and associated shape statistics with vertices as iterators to an fstate set, ie graph over f-space only
*/

#include <sstream>  // to be able to manipulate strings as streams

//#include <ebc_output.hpp>
//#include <ebc_params.hpp>

#include <ebc_fseq_efficient_graph.hpp>
#include <ebc_fxzseq_named_graph.hpp>
#include <ebc_fseq_general_graph.hpp>
#include <ebc_output.hpp>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_matrix.h>
//#include <gsl/gsl_rng.h> // is in the hpp file

#include <boost/graph/graphviz.hpp>

#include <numeric> // to use accumulate
#include <limits> // numeric limits

#define SPLITS_DEBUG 0
#define SPLITS_DEBUG1 0
#define TREEOUT_DEBUG 0
#define TIMESOUT_DEBUG 1
#define FSEQOUT_DEBUG 1


/*! generate StateSequences of Kingman's unlabeled n-coalescent (i.e. the unvintaged and sized n-coalescent) forward in time

Using a BGL graph
*/

// generate a controlled f-sequence, using the z state from some data
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
                                        const bool efficientFind)
{
    if(nsam < 2) {
        std::cerr << "\nnsam must be > 1 in GenerateStateSequence\n";
        exit(0);
    }

    // fs is as in Raaz's program, which misses out the 'top' row, epoch 1,
    // where one lineage subtends nsam leaves, and therefore the matrix form
    // is (nsam-1)x(nsam-1)and only has epochs 2..nsam, and (the columns)
    // leaves subtended 1..nsam-1

    //make the the fseq fs is clear
    fs.clear();


    size_t epoch = 1; // starting graph at epoch 1

    // make sure that the first (epoch 1) vertex is in the graph
    f_efficient_graph_vertex_descriptor firstV = firstVertex(graph,
                                            extStateVertexMap,
                                            extFStateSet,
                                            extZStateSet,
                                            nsam,
                                            z_state);

    fsProposalProb = 1.0; // moving to this first state with certainty
    fsKCProb = 1.0;

    // the first state is not shown in the f_seq_type fs

    // take it that we are here, and make the next move
    // when we start we assume that we are on a known path
    bool completed = nextStateSequenceState(graph, extStateVertexMap,
                                extFStateSet, extZStateSet, rgsl,
                                fs, fsProposalProb, fsKCProb, nsam,
                                firstV, z_state, epoch);

    //debug
    if (completed) {
        /*
        std::cout << "\nf-sequence found is " << std::endl;
        printStateSequence(fs);
*/
        // std::cout << "printing map \n" << std::endl;
        //print_graph_desc(graph);
    }
    else {
        std::cout << "failed to find StateSequence" << std::endl;
        std::cout << std::endl;
    }


    return completed;
}



bool generateFSequenceWithGraph(f_efficient_graph& graph,
                                ExtSetItVertexMap& extStateVertexMap,
                                state_typeSet& extFStateSet,
                                state_typeSet& extZStateSet,
                                gsl_rng* rgsl,
                                StateSequence& fs,
                                double& fsProposalProb,
                                double& fsKCProb,
                                const size_t nsam,
                                const bool efficientFind)
{
    size_t initVal = 0;
    //uncontrolled z_state z1
    state_type z1_uncontrolled(nsam-1, initVal); // nsam-1 elements, all 0

    return generateControlledFSequenceWithGraph(graph, extStateVertexMap,
            extFStateSet, extZStateSet, rgsl, fs, z1_uncontrolled,
            fsProposalProb, fsKCProb, nsam);
}


// generate the complete graph with an initial controlling z_state
// graph will only include vertices with f_states compatible with given z_state
bool generateCompleteControlledGraph(f_efficient_graph& graph,
                                    ExtSetItVertexMap& extStateVertexMap,
                                    state_typeSet& extFStateSet,
                                    state_typeSet& extZStateSet,
                                    const size_t nsam,
                                    const state_type& z_state)
{
    if(nsam < 2) {
        std::cerr << "\nnsam must be > 1 in GenerateStateSequence\n";
        exit(0);
    }

    if(nsam > 80) {
        std::cerr << "\nnsam should be <= 80 for this graph type, or you'll use too much RAM\n";
        std::cerr << "Try using the more efficient f-space only graph instead\n" << std::endl;
        exit(0);
    }

    size_t epoch = 1; // starting graph at epoch 1

    // make sure that the first (epoch 1) vertex is in the graph
    f_efficient_graph_vertex_descriptor firstV = firstVertex(graph,
                                            extStateVertexMap,
                                            extFStateSet,
                                            extZStateSet,
                                            nsam,
                                            z_state);

    // take it that we are here, and make the next move
    bool completed = allNextStateSequenceStates(graph, extStateVertexMap,
                                            extFStateSet, extZStateSet, nsam,
                                            firstV, z_state, epoch);

    if (completed) {

    }
    else
        std::cout << "Could not successfully complete graph" << std::endl;

    return completed;
}


// generate the complete graph
// no initial controlling z_state
bool generateCompleteGraph(f_efficient_graph& graph,
                            ExtSetItVertexMap& extStateVertexMap,
                            state_typeSet& extFStateSet,
                            state_typeSet& extZStateSet,
                            const size_t nsam)
{
    size_t initVal = 0;
    //uncontrolled z_state z1
    state_type z1(nsam-1, initVal); // nsam-1 elements, all 0

    return generateCompleteControlledGraph(graph, extStateVertexMap,
                                            extFStateSet, extZStateSet,
                                            nsam, z1);
}


// make sure that the first vertex in the graph is there
// this method makes the first vertex if it is not there already
f_efficient_graph_vertex_descriptor firstVertex(f_efficient_graph& graph,
                                                    ExtSetItVertexMap& extStateVertexMap,
                                                    state_typeSet& extFStateSet,
                                                    state_typeSet& extZStateSet,
                                                    const size_t nsam,
                                                    const state_type& z_state)
{

    // the vertex for EPOCH 1 is the pair (f_1 = (0 ... 0 1), z_state)
    size_t initVal = 0;
    //state_type f_1(initVal, nsam); // nsam elements, all 0
    state_type f_1(nsam, initVal); // nsam elements, all 0
    f_1[nsam-1] = 1; // set nsam-th element to 1

    // make sure the f_state is in the external map
    state_typeSetBool pair_f = extFStateSet.insert(f_1);

    // make sure the z_state is in the external map
    state_typeSetBool pair_z = extZStateSet.insert(z_state);

    // check vertex representing f_state is there in the map or make it if not
    F_EfficientGraphVertexBoolPair addVertex = vertexAddOrFind(graph,
                                    extStateVertexMap,
                                    pair_f.first);

    return addVertex.first;
}


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
                        const size_t epoch)
{
    bool retValue = false;

    // now make sure all the next states are there
    // (add them if not there already)
    bool expanded = false;

    // if the vertex has some out edges,it will have all out-edges, so check this
    fxz_efficient_degree_size_type outs = out_degree(v, graph);
    if (outs > 0) { // v has out edges

        expanded = true;
    }
    else { // try to expand, with a 'dummy' z_state which gives no control

        size_t initVal = 0;
        //uncontrolled z_state z1
        state_type z1(nsam-1, initVal); // nsam-1 elements, all 0

        // set up some containers required by expand subroutine but contents not
        // in fact used here (we may not be able to get to all the f_states
        // given our z_state - rely on chooseNextVertex to find out where we can go
        VecSetIt new_f_statesSetIt;
        VecSetIt new_z_statesSetIt;

        expanded = expandGraphKingmanForward(graph, extStateVertexMap,
                        extFStateSet, extZStateSet,
                        new_f_statesSetIt,
                        new_z_statesSetIt,
                        nsam, v, z1, epoch);
    }

    if (expanded) {

        // choose which next state to go to from our current vertex addVertex.first
        // given our current z_state
        // the reweighting of transition probabilities from that for the vanilla f_state
        // to that for the controlled chain takes place here so that the transition
        // probability returned in the tuple is the reweighted one;
        F_EfficientVertexChoiceTuple choiceTuple = chooseNextVertex(graph, rgsl, nsam,
                                            v, z_state);
        // nextVertex, proposal transition probability, kc transition probability, haveVertex

        bool haveVertex = boost::get<3>(choiceTuple); // false if no vertex could be found

        if (haveVertex) { // if we found a vertex

            f_efficient_graph_vertex_descriptor nextVertex = get<0>(choiceTuple);

            // fold the proposal transition probability into fsProposalProb
            fsProposalProb *= get<1>(choiceTuple);

            // fold the kc transition probability into fsKCProb
            fsKCProb *= get<2>(choiceTuple);

            // we will need the fstate property map
            f_efficient_fstateIt_map_t fstateProp = get(vertex_fstateSetIt, graph);

            // next epoch
            size_t newEpoch = epoch+1;

            // use the internal property maps to get the current and new f states
            state_typeSetIt new_f_stateSetIt = fstateProp[nextVertex];
            state_typeSetIt current_f_stateSetIt = fstateProp[v];

            state_type new_f_state = *new_f_stateSetIt;
            state_type f_state = *current_f_stateSetIt;

            // add the new f_state to the f_seq_container fs
            addStateToF_Sequence(fs, new_f_stateSetIt, newEpoch, nsam);

            // now we have to work out what out new z_state is given the current
            // and next f_states
            state_type new_z_state = z_state;

            bool foundNewZstate = findNextZState(nsam, f_state, new_f_state,
                                                z_state, new_z_state);

            if (foundNewZstate) {

                //make sure the new z state is in the external map
                state_typeSetBool pair_z = extZStateSet.insert(new_z_state);

                // make the next move
                retValue = nextStateSequenceState(graph, extStateVertexMap,
                                                extFStateSet, extZStateSet, rgsl,
                                                fs, fsProposalProb, fsKCProb,
                                                nsam, nextVertex,
                                                new_z_state, newEpoch);
            }
        }
        else { // did not find a vertex to go to
            std::cout << "Problem finding f-sequence ";
            std::cout << "in transition from epoch " << epoch << std::endl;
        }
    }

    else { // graph was not expanded
    // check that is is because we are in the final state

        f_efficient_fstateIt_map_t fstateProp = get(vertex_fstateSetIt, graph);
        state_typeSetIt f_stateSetIt = fstateProp[v];
        state_type f_state = *f_stateSetIt;

        // double check if the current f_state is the final state
        if ((epoch == nsam) && checkFinalState(f_state, nsam))
            retValue = true;
        else
            std::cout << "Error trying to expand graph" << std::endl;
    }

    return retValue;
}



// add all next FZSequence states from the current state
// with transition probabilities reweighted for current z_state
bool allNextStateSequenceStates(f_efficient_graph& graph,
                            ExtSetItVertexMap& extStateVertexMap,
                            state_typeSet& extFStateSet,
                            state_typeSet& extZStateSet,
                            const size_t nsam,
                            const f_efficient_graph_vertex_descriptor v,
                            const state_type& z_state,
                            const size_t epoch)
{

    bool retValue = false;

    // now make sure all the next states are there
    // (add them if not there already)
    bool expanded = false;

    // if the vertex has some out edges,it will have all out-edges, so check this
    f_efficient_degree_size_type outs = out_degree(v, graph);
    if (outs > 0)  { // already expanded - nothing here for us to do
        retValue = true;

    }
    else { // v has no out edges yet

        // set up some variables to be filled by next subroutine
        VecSetIt new_f_statesSetIt;
        VecSetIt new_z_statesSetIt;

        expanded = expandGraphKingmanForward(graph, extStateVertexMap,
                        extFStateSet, extZStateSet,
                        new_f_statesSetIt,
                        new_z_statesSetIt,
                        nsam, v, z_state, epoch);

        if (expanded) {

            // increment the epoch counter
            size_t newEpoch = epoch + 1;

            // iterate through the next states from this vertex and expand them
            // but only if we have not now reached the final state
            if (newEpoch < nsam) {

                VecSetIt_it fit;
                VecSetIt_it zit;

                for (fit = new_f_statesSetIt.begin(), zit = new_z_statesSetIt.begin();
                        fit < new_f_statesSetIt.end(), zit < new_z_statesSetIt.end();
                        fit++, zit++) {

                    retValue = true; // found at least one state to go to

                    //find the vertex
                    f_efficient_graph_vertex_descriptor nextVertex =
                                            extStateVertexMap[(*fit)];

                    state_type new_z_state = *(*zit);

                    // use allNextFStates recursively, using retValue to ensure that
                    // we always get to the final state
                    retValue = (retValue
                                && allNextStateSequenceStates(graph,
                                                    extStateVertexMap,
                                                    extFStateSet,
                                                    extZStateSet,
                                                    nsam,
                                                    nextVertex,
                                                    new_z_state,
                                                    newEpoch));

                    if (!retValue) {
                        std::cout << "Cannot get to final state " << std::endl;
                        std::cout << "epoch is " << newEpoch << std::endl;

                        break; // break out of for loop if failed to reach final state

                    }
                }


            } // end second check on epoch
            else
                retValue = true; // at final epoch

        }
        else { // graph was not expanded
            // check that is is because we are in the final state

            f_efficient_fstateIt_map_t fstateProp = get(vertex_fstateSetIt, graph);
            state_typeSetIt f_stateSetIt = fstateProp[v];
            state_type f_state = *f_stateSetIt;

            // double check if the current f_state is the final state
            if ((epoch == nsam) && checkFinalState(f_state, nsam))
                retValue = true;
            else
                std::cout << "Error trying to expand graph" << std::endl;
        }
    }

    return retValue;
}


// controlled expansion
// find all possible next state pairs given a current state pair
// add them to the graph and the external maps
// return false if we could not add any states, true otherwise
bool expandGraphKingmanForward(f_efficient_graph& graph,
                                ExtSetItVertexMap& extStateVertexMap,
                                state_typeSet& extFStateSet,
                                state_typeSet& extZStateSet,
                                VecSetIt& new_f_statesSetIt,
                                VecSetIt& new_z_statesSetIt,
                                const size_t nsam,
                                const f_efficient_graph_vertex_descriptor v,
                                const state_type& z_state,
                                const size_t epoch)
{
    bool retValue = false;

    // use the internal property map to get the current f state iterator
    f_efficient_fstateIt_map_t fstateProp = get(vertex_fstateSetIt, graph);
    state_typeSetIt f_stateSetIt = fstateProp[v];
    // and hence the f_state
    state_type f_state = *f_stateSetIt;

    // if the current f_state is the final state, we can't do anything
    // and transition probabilities depend on the state we are entering, and
    // we can check the current epoch against the current state
    //if ((f_state.sum() == epoch) && !checkFinalState(f_state, nsam)) {
    if ((epoch < nsam) && !checkFinalState(f_state, nsam)) {

        size_t newEpoch = epoch + 1;

        // set up some variables to be filled by next subroutine
        //VecSetIt new_f_statesSetIt;
        //VecSetIt new_z_statesSetIt;
        VecProbs newProbs;
        double accProbs = 0.0; // use to accumulate probabilities of all transitions
        double rejectedProbs = 0.0; // use to accumulate rejected transitions

        retValue = findNewStates(extFStateSet, extZStateSet, f_state, z_state,
                                    nsam, newEpoch, accProbs, rejectedProbs,
                                    new_f_statesSetIt, new_z_statesSetIt, newProbs);

        // go through the new f states
        if (retValue && !new_f_statesSetIt.empty()) {

            VecSetIt_it fit;
            VecProbsIt dit;

            for (fit = new_f_statesSetIt.begin(), dit = newProbs.begin();
                fit < new_f_statesSetIt.end(), dit < newProbs.end();
                fit++, dit++) {

                // make sure new state is in the graph
                F_EfficientGraphVertexBoolPair addVertex = vertexAddOrFind(graph,
                                                extStateVertexMap,
                                                *fit);

                // make sure edge is in the graph
                // with reweighted probability
                F_EfficientGraphEdgeBoolPair addEdge = edgeAddOrFind(graph,
                                v, addVertex.first, ((*dit)/accProbs));

            }
        }

    } // end if not final state

    return retValue;
}


// choose an edge or next vertex to go to from our current vertex v
// the reweighting of transition probabilities from that for the vanilla f_state
// to that for the controlled chain takes place here so that the transition
// probability returned in the tuple is the reweighted one;
// tuple returned is (vertex, proposal transition probability, kc transition probability, bool)
F_EfficientVertexChoiceTuple chooseNextVertex(f_efficient_graph& graph,
                                                gsl_rng* rgsl,
                                                const size_t nsam,
                                                const f_efficient_graph_vertex_descriptor v,
                                                const state_type& z_state)
{
    bool haveVertex = false;
    double transProposalProb = 0;
    double transKCProb = 0;

    f_efficient_graph_vertex_descriptor nextVertex = v; // initialise to this vertex

    // set up some variables which will be given values by being passed by reference to next subroutine
    F_EfficientVecVertices vertices;
    VecProbs vecProbs;
    double totalAllowedProb = 0;

    // fill in these variables
    // where we can go depends on the z_state
    findAllowableNextVertices(graph, nsam, v, z_state,
                            vertices, vecProbs, totalAllowedProb);

    if (totalAllowedProb > 0 && !vertices.empty()) {

        // choose a vertex as the next vertex
        // choose a random number between 0 and 1
        // use this to choose which vertex to go to
        double rand = gsl_rng_uniform(rgsl);

        // scale it by the total probability of the available paths
        rand *= totalAllowedProb;

        double accProbs = 0;
        // go through the available vertices and choose one
        F_EfficientVecVerticesIt vit;
        VecProbsIt dit;

        for (vit = vertices.begin(), dit = vecProbs.begin();
                vit < vertices.end(), dit < vecProbs.end();
                vit++, dit++) {

            double prob = *dit;
            accProbs += prob;

            if (rand < accProbs) {
                nextVertex = *vit;
                haveVertex = true;
                transProposalProb = prob/totalAllowedProb; // reweighting probabilities
                transKCProb = prob; // pure KC transition probability

                break; // found a target - break out of for loop

                //debug
                //string n1 = get(nameProp, v);
                //string n2 = get(nameProp, nextVertex);
                //std::cout << "Chosen edge from " << n1 << " to "
                        //<< n2 << " with KC proby " << prob << std::endl;

            }
        } // end for loop through available vertices
    }
    // else there is nowhere to go

    F_EfficientVertexChoiceTuple choiceTuple(nextVertex, transProposalProb,
                                        transKCProb, haveVertex);

    return choiceTuple;
}



// find the allowable next vertices from this one given the z_state
// fill in vecVertices, vecProbs, totalAllowedProbs which are passed by reference
void findAllowableNextVertices(f_efficient_graph& graph,
                                const size_t nsam,
                                const f_efficient_graph_vertex_descriptor v,
                                const state_type& z_state,
                                vector<f_efficient_graph_vertex_descriptor>& vertices,
                                vector<double>& vecProbs,
                                double& totalAllowedProb)
{
    // use the internal property map to get the current f state iterator
    f_efficient_fstateIt_map_t fstateProp = get(vertex_fstateSetIt, graph);
    state_typeSetIt f_stateSetIt = fstateProp[v];
    // and hence the f_state
    state_type f_state = *f_stateSetIt;

    // will need the KC edge probabilities map
    f_efficient_kc_prob_map_t edgeKCProbProp = get(edge_kcprob, graph);

    // check whether the z state places any restrictions on where we can go
    size_t init = 0;
    size_t zSum = std::accumulate(z_state.begin(), z_state.end(), init);

    // z tells us about states we must have:  if z has an entry in the
    // position corresponding to *leaves being subtended, then
    // the f sequences must include a state where there is
    // at least one lineage subtending *leaves

    // find the maximum such restriction
    size_t maxRestrictingL = 0;
    size_t ll = z_state.size()-1;

    while (zSum > 0 && ll > 0 && maxRestrictingL == 0) {
        if (z_state[ll] != 0)
            maxRestrictingL = ll+1;
        ll--;
    }

    // go through all the edges from here
    // go through the out-edges of this vertex and add up the probabilities
    // of the available edges
    for (pair<f_efficient_graph_out_edge_it, f_efficient_graph_out_edge_it> ed
                        = out_edges(v, graph);
                    ed.first != ed.second; ++ed.first) {

        // find the target vertex
        f_efficient_graph_vertex_descriptor trg = target(*(ed.first), graph);
        // find the f_state of the target vertex
        state_typeSetIt sit = get(fstateProp, trg);
        state_type nextFstate = *sit;

        // what was the split to get here?
        size_t leavesSubtended = 0;
        size_t split_j = 0;
        size_t split_k = 0;
        size_t i = 0;
        while ((split_j == 0) && (i < nsam)) {
            if ((f_state[i] == nextFstate[i] - 1) || (f_state[i] == nextFstate[i] - 2)) {
                split_j = i+1;
            }
            i++;
        }
        i = split_j;

        while ((split_j > 0) && (leavesSubtended == 0) && (i < nsam)) {
            if (f_state[i] == nextFstate[i] + 1) {
                leavesSubtended = i+1;
                split_k = leavesSubtended - split_j;
            }
            i++;
        }

        if (leavesSubtended > 0) {

            size_t split_jmax = leavesSubtended/2; // integer division to get floor;

            size_t split_jhat = split_jmax;

            // split_jhat is only relevant if there is some restriction from z
            if (maxRestrictingL > 0) {
                split_jhat = leavesSubtended - 1;
                // split_jhat = min(leavesSubtended - 1, maxRestrictingL)
                if (maxRestrictingL < split_jhat)
                    split_jhat = maxRestrictingL;
                if (leavesSubtended % 2 == 0) {
                    // floor and ceiling of leavesSubtended/2 will be same
                    // split_jhat = max(split_jhat, split_jmax)
                    if (split_jhat < split_jmax)
                        split_jhat = split_jmax;
                }
                else {
                    // ceiling of leavesSubtended/2 = floor +1
                    // split_jhat = max(split_jhat, split_jmax+1)
                    if (split_jhat < split_jmax + 1)
                        split_jhat = split_jmax + 1;
                }

            }

            //but is the new state allowed?
            bool allowed = ((maxRestrictingL == 0) ||
                    ((split_j <= split_jhat) && (split_jhat <= split_k)));

            if (allowed) {

                vertices.push_back(trg);
                double prob = get(edgeKCProbProp, *(ed.first));
                vecProbs.push_back(prob);
                totalAllowedProb += prob;

            }

        } // end leavesSubtended > 0
    } // end for loop through edges out of this vertex v
}


// find the z_state that will result from moving from f_state to new_f_state under z_state
// changes new_z_state which is passed by reference
// returns false if no z state could be found
bool findNextZState(const size_t nsam,
                    const state_type& f_state,
                    const state_type& new_f_state,
                    const state_type& z_state,
                    state_type& new_z_state)
{
    bool retValue = false;

    // what was the split to get here?
    size_t leavesSubtended = 0;
    size_t split_j = 0;
    size_t split_k = 0;
    size_t i = 0;
    while ((split_j == 0) && (i < nsam)) {
        if ((f_state[i] == new_f_state[i] - 1) || (f_state[i] == new_f_state[i] - 2)) {
            split_j = i+1;
        }
        i++;
    }
    i = split_j;

    while ((split_j > 0) && (leavesSubtended == 0) && (i < nsam)) {
        if (f_state[i] == new_f_state[i] + 1) {
            leavesSubtended = i+1;
            split_k = leavesSubtended - split_j;
        }
        i++;
    }

    if (leavesSubtended > 0) {
        retValue = true;
        // z is changed by turning any 1 at
        // index[split_j-1[ [split_k-1] to 0
        new_z_state = z_state;
        if (new_z_state[split_j-1] == 1)
            new_z_state[split_j-1] = 0;
        if (new_z_state[split_k-1] == 1)
            new_z_state[split_k-1] = 0;
    }

    return retValue;
}




// add a new vertex or find the vertex description if its already there
// return (f_efficient_graph_vertex_descriptor, bool)
// vertex descriptor is the descriptor for the vertex,
// bool is true if new vertex was inserted, false if it was already there
F_EfficientGraphVertexBoolPair vertexAddOrFind(f_efficient_graph& graph,
                                ExtSetItVertexMap& extStateVertexMap,
                                const state_typeSetIt& f_stateSetIt)
{

    bool inserted = false;

    f_efficient_graph_vertex_descriptor v;

    if (extStateVertexMap.find(f_stateSetIt) == extStateVertexMap.end()) {

        // not found so need to add to graph and maps
        v = add_vertex(graph);
        extStateVertexMap[f_stateSetIt] = v;
        f_efficient_fstateIt_map_t fstateProp = get(vertex_fstateSetIt, graph);
        fstateProp[v] = f_stateSetIt;
        inserted = true;

    }

    else { // already there, just get the description from the external map
        v = extStateVertexMap[f_stateSetIt];

    }

    F_EfficientGraphVertexBoolPair retPair(v, inserted);

    return retPair;
}



// add a new edge or find the edge description if its already there
// the probability is added as a property of the edge
// and also the edge weight property becomes ln(1/probability)
// but if probability = 0, the edge weight = MAX_WEIGHT (defined in ebc_graphtypes.hpp)
// return (f_efficient_graph_edge_descriptor, bool)
// edge descriptor is the descriptor for the edge,
// bool is true if new edge was inserted, false if it was already there
F_EfficientGraphEdgeBoolPair edgeAddOrFind(f_efficient_graph& graph,
                                    const f_efficient_graph_vertex_descriptor v1,
                                    const f_efficient_graph_vertex_descriptor v2,
                                    const double edgeKCprob)
{
    f_efficient_graph_edge_descriptor ed;
    f_efficient_kc_prob_map_t edgeKCProbProp = get(edge_kcprob, graph);
    // we don't use reweighted probabilities
    f_efficient_weight_map_t edgeWeightProp = get(edge_weight, graph);


    F_EfficientGraphEdgeBoolPair p = add_edge(v1, v2, graph);
    if (p.second) { // new edge was inserted

        edgeKCProbProp[p.first] = edgeKCprob;

        if (edgeKCprob > 0.0)
            edgeWeightProp[p.first] = log(1.0/edgeKCprob); // natural log

        // weight is log of 1/probability except if probability == 0.0
        else // MAX_WEIGHT defined in ebc_graphtypes.hpp
            edgeWeightProp[p.first] = 1.0*MAX_WEIGHT;

        // no followed property for f_efficient_graph

    }

    return p;
}


// print description of the graph
void print_graph_desc(f_efficient_graph& graph)
{
    // check the number of vertices
    f_efficient_vertices_size_type numVert = num_vertices (graph);

    if (numVert > 500) {
        std::cout << "\n\nThe graph has more than 500 vertices and will not be printed\n" << std::endl;
    }
    else {

        f_efficient_fstateIt_map_t fstateProp = get(vertex_fstateSetIt, graph);

        f_efficient_kc_prob_map_t edgeKCProbProp = get(edge_kcprob, graph);

        std::cout << "edges(graph) = \n";

        for (pair<f_efficient_graph_edge_it, f_efficient_graph_edge_it> ed = edges(graph);
                            ed.first != ed.second; ++ed.first) {

            f_efficient_graph_vertex_descriptor src = source(*(ed.first), graph);
            f_efficient_graph_vertex_descriptor tgt = target(*(ed.first), graph);
            state_typeSetIt src_f_stateSetIt = get(fstateProp, src);
            state_typeSetIt tgt_f_stateSetIt = get(fstateProp, tgt);
            state_type src_f_state = *src_f_stateSetIt;
            state_type tgt_f_state = *tgt_f_stateSetIt;
            string src_label = stateToString(src_f_state);
            string tgt_label = stateToString(tgt_f_state);

            double prob = edgeKCProbProp[*(ed.first)];

            std::cout << "Edge from " << src_label << " to " << tgt_label << " has kc proby " << prob << std::endl;
        }

    }

}


void bfs_accumulate_count_visitor::discover_vertex(f_efficient_graph_vertex_descriptor u, const f_efficient_graph& g)
{
    if (good) {
        size_t thisAccumulation = 0;
        // if this vertex has no in-edges, give it an accumulation of 1
        if (in_degree(u, g) == 0) {
            accMap.insert(make_pair(u, 1));
        }
        // else (vertex has in-edges)
        else {
            // find the source vertices of all the in edges of this vertex
            for (pair<f_efficient_graph_in_edge_it, f_efficient_graph_in_edge_it> ed
                        = in_edges(u, g);
                    ed.first != ed.second; ++ed.first) {
                size_t acc = accMap[source(*(ed.first), g)];
                // add up their accumulations
                if (thisAccumulation < std::numeric_limits<size_t>::max() - acc) {
                    thisAccumulation += acc;
                }
                else {
                    cout << "\nError, exceeded maximum size_t that can be held ";
                    cout << std::numeric_limits<size_t>::max() << endl;
                    good = false;
                }

            }

            // make this the accumulation for this vertex
            accMap.insert(make_pair(u, thisAccumulation));
        }

        // if this vertex has no out-edges
        if (out_degree(u, g) == 0) {
            // update the value stored at *countPtr
            *countPtr = thisAccumulation;
        }
    }


}


//calculating number of fsequences using bread-first search and bfs_accumulate_count_visitor
size_t  calcFsequences(f_efficient_graph& graph)
{
    size_t finalCount = 0;

    if (num_vertices(graph) > 0) {
        // take the source as the first vertex in the graph
        pair<f_efficient_graph_vertex_it, f_efficient_graph_vertex_it> p = vertices(graph);
        f_efficient_graph_vertex_descriptor sourceVertex = *(p.first);;

        bfs_accumulate_count_visitor vis(&finalCount);

        breadth_first_search(graph, sourceVertex, visitor(vis));

        // vis's countPtr should now have the final count

    }
    else
        cout << "Sorry the graph has no vertices" << endl;

    return finalCount;

}


void useGraphToCountFsequences (const size_t max_n,
                                const size_t min_n,
                                bool doXstars,
                                const string& filename)
{
    const size_t init = 0; // an initialiser for values in state_types

    // map from padded xstars to maps for n
    PaddedXStarToN_FSequenceCountMap paddedXStarMap;

    // count map to copy for insertion
    N_FSequenceCountMap countMap;


    for (size_t n = min_n; n < max_n + 1; n++) {

        cout << "\n************* n = " <<  n << " **************" << endl;

        state_type xstar(n-1, init); // first xstar, all 0's

        state_type padding((max_n - n), init); // make padding of 0's

        // make a container to hold possible xstars for this n
        state_typeSet xstars;

        if (doXstars) {
            // and fill it up
            xstars = makePossibleXStars(xstars, xstar, 2);
        }
        else { // just do this (uncontrolled) xstar
            state_typeSetBool insertBool = xstars.insert(xstar);
        }

        // go through each of the possible x stars
        for (state_typeSetIt it = xstars.begin(); it != xstars.end(); it++) {

            // set up structures we need to create graph
            f_efficient_graph graph; // a graph to count sequences

            state_typeSet extFStateSet; // set of unique f states
            state_typeSet extZStateSet; // set of unique z states

            ExtSetItVertexMap extStateVertexMap; // map from set iterators to vertices

            // make the padded version of the xstar, padding to the size of max_N-1
            state_type padded_xstar = *it;
            padded_xstar.insert(padded_xstar.end(), padding.begin(), padding.end());

            string xstar_label = stateToString(*it);

            cout << "\txstar " <<  xstar_label << endl;

            bool map = generateCompleteControlledGraph(graph, extStateVertexMap,
                        extFStateSet, extZStateSet, n, *it);

            if (map) {
                size_t finalCount = calcFsequences(graph);
                cout << "\tThe f-sequence count is\t" <<  finalCount << "\n\n" << endl;

                // find or make map for this padded xstar in the overall map
                PaddedXStarToN_FSequenceCountMapBool xstarBool
                        = paddedXStarMap.insert(make_pair(padded_xstar, countMap));

                // put the n and finalCount into the count map for this padded xstar
                ((xstarBool.first)->second).insert(make_pair(n, finalCount));
            }
            else {
                cout << "Error creating graph for n" << n
                    << " and xstar " << xstar_label << endl;
            }

        }

    }
    // should now have a map from padded xstars to n's and their f sequence counts
    sendToOutfile(paddedXStarMap, min_n, max_n, filename);

}





