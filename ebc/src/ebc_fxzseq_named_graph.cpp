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
/*! \file ebc_fxzseq__named_graph.cpp
    \brief BGL methods for unvintaged (sized) n-coalescent and associated shape statistics using a named graph over the f x z product space which is inefficient on memory but good for images for small nsam
*/

#include <sstream>  // to be able to manipulate strings as streams
#include <iostream>
#include <iomanip>

//#include <ebc_output.hpp>
//#include <ebc_params.hpp>

#include <ebc_fxzseq_named_graph.hpp>
#include <ebc_fseq_general_graph.hpp>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_matrix.h>
//#include <gsl/gsl_rng.h> // is in the hpp file

#include <boost/graph/graphviz.hpp>


#include <numeric> // to use accumulate

#define SPLITS_DEBUG 0
#define SPLITS_DEBUG1 0
#define TREEOUT_DEBUG 0
#define TIMESOUT_DEBUG 1
#define FSEQOUT_DEBUG 1


/*! generate fsequences of Kingman's unlabeled n-coalescent (i.e. the unvintaged and sized n-coalescent) forward in time

Using a BGL graph
*/

// generate a controlled f-sequence, using the z state from some data
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
                                        const bool efficientFind,
                                        const string addColour)
{
    if(nsam < 2) {
        std::cerr << "\nnsam must be > 1 in GenerateFSequence\n";
        exit(0);
    }

    if(nsam > 30) {
        std::cerr << "\nnsam should be <= 30 for the named graph, or you'll use too much RAM\n";
        std::cerr << "Try using one of the efficient graphs instead\n" << std::endl;
        exit(0);
    }

    // fs is as in Raaz's program, which misses out the 'top' row, epoch 1,
    // where one lineage subtends nsam leaves, and therefore the matrix form
    // is (nsam-1)x(nsam-1)and only has epochs 2..nsam, and (the columns)
    // leaves subtended 1..nsam-1

    //make the the fseq fs is clear
    fs.clear();

    // also make ourselves a container for z_states from epoch 2 onwards
    StateSequence zs;

    size_t epoch = 1; // starting graph at epoch 1

    // make sure that the first (epoch 1) vertex is in the graph
    fxz_named_graph_vertex_descriptor firstV = firstVertex(graph,
                                            extStatesPairVertexMap,
                                            extFStateSet, extZStateSet,
                                            nsam, z_state);

    fsProposalProb = 1.0; // moving to this first state with certainty
    fsKCProb = 1.0;

    // the first state is not shown in the f_seq_type fs

    // take it that we are here, and make the next move
    // when we start we assume that we are on a known path
    bool completed = nextFZSequenceState(graph, extStatesPairVertexMap,
                                extFStateSet, extZStateSet, rgsl,
                                fs, zs, fsProposalProb, fsKCProb, nsam,
                                firstV, epoch, efficientFind, addColour);

    //debug
    if (completed) {
/*
        std::cout << "\nf-sequence found is " << std::endl;

        printStateSequence(fs);

        size_t init = 0; // only print z sequence if there was a control
        if (std::accumulate(z_state.begin(), z_state.end(), init) > 0) {
            std::cout << "\nUsing z-sequence " << std::endl;

            printStateSequence(zs);
        }
*/
/*
       // std::cout << "printing map \n" << std::endl;
        fxz_named_name_map_t nameProp = get(vertex_name, graph);
        std::cout << "printing usng print_vertices" << endl;
        print_vertices(graph, nameProp);
        std::cout << "printing usng print_edges" << endl;
        print_edges(graph, nameProp);
    */
        //print_graph_desc(graph);
    }
    else {
        std::cout << "failed to find fsequence" << std::endl;
        std::cout << std::endl;
    }


    return completed;
}



bool generateFSequenceWithGraph(f_x_z_named_graph& graph,
                                ExtSetItPairVertexMap& extStatesPairVertexMap,
                                state_typeSet& extFStateSet,
                                state_typeSet& extZStateSet,
                                gsl_rng* rgsl,
                                StateSequence& fs,
                                double& fsProposalProb,
                                double& fsKCProb,
                                const size_t nsam,
                                const bool efficientFind,
                                const string addColour)
{
    size_t initVal = 0;
    //uncontrolled z_state z1
    state_type z1_uncontrolled(nsam-1, initVal); // nsam-1 elements, all 0

    return generateControlledFSequenceWithGraph(graph, extStatesPairVertexMap,
            extFStateSet, extZStateSet,
            rgsl, fs, z1_uncontrolled, fsProposalProb, fsKCProb,
            efficientFind, nsam);
}


// generate the complete graph with an initial controlling z_state
bool generateCompleteControlledGraph(f_x_z_named_graph& graph,
                                    ExtSetItPairVertexMap& extStatesPairVertexMap,
                                    state_typeSet& extFStateSet,
                                    state_typeSet& extZStateSet,
                                    const size_t nsam,
                                    const state_type& z_state)
{
    if(nsam < 2) {
        std::cerr << "\nnsam must be > 1 in GenerateFSequence\n";
        exit(0);
    }

    if(nsam > 30) {
        std::cerr << "\nnsam should be <= 30 for the named graph, or you'll use too much RAM\n";
        std::cerr << "Try using one of the efficient graphs instead\n" << std::endl;
        exit(0);
    }

    size_t epoch = 1; // starting graph at epoch 1

    // make sure that the first (epoch 1) vertex is in the graph
    fxz_named_graph_vertex_descriptor firstV = firstVertex(graph,
                                            extStatesPairVertexMap,
                                            extFStateSet, extZStateSet,
                                            nsam, z_state);

    // take it that we are here, and make the next move
    bool completed = allNextFZSequenceStates(graph, extStatesPairVertexMap,
                                            extFStateSet, extZStateSet, nsam,
                                            firstV, epoch);

    if (completed) {
/*
        size_t init = 0; // only print z sequence if there was a control
        if (std::accumulate(z_state.begin(), z_state.end(), init) > 0) {

            std::cout << "\nInitial z_state was ";
            copy (z_state.begin(), z_state.end(), ostream_iterator<size_t>(cout, "\t"));
            cout << "\n";
        }

        std::cout << "printing map \n" << std::endl;
        print_graph_desc(graph);

        // get the number of vertices
        fxz_named_vertices_size_type numVert = num_vertices (graph);
        std::cout << "\n\nThe graph has " << numVert << " vertices\n" << std::endl;

*/

      // get the number of vertices

/*
        std::cout << "vertices(graph) f states= \n";
        fxz_named_fstate_map_t fstateProp = get(vertex_fstate, graph);
        for (pair<fxz_named_graph_vertex_it, fxz_named_graph_vertex_it> p = vertices(graph);
                            p.first != p.second; ++p.first) {

            state_type f_state = get(fstateProp, *(p.first));

            copy (f_state.begin(), f_state.end(), ostream_iterator<size_t>(cout, "\t"));
            cout << "\n";
        }

*/

    /*


        std::cout << "printing usng print_vertices" << endl;
        print_vertices(graph, nameProp);
        std::cout << "printing usng print_edges" << endl;
        print_edges(graph, nameProp);
    */

    }
    else
        std::cout << "Could not successfully complete graph" << std::endl;

    return completed;
}


// generate the complete graph
// no initial controlling z_state
bool generateCompleteGraph(f_x_z_named_graph& graph,
                            ExtSetItPairVertexMap& extStatesPairVertexMap,
                            state_typeSet& extFStateSet,
                            state_typeSet& extZStateSet,
                            const size_t nsam)
{
    size_t initVal = 0;
    //uncontrolled z_state z1
    state_type z1(nsam-1, initVal); // nsam-1 elements, all 0

    return generateCompleteControlledGraph(graph, extStatesPairVertexMap,
                                            extFStateSet, extZStateSet,
                                            nsam, z1);
}


// make sure that the first vertex in the graph is there
// this method makes the first vertex if it is not there already
fxz_named_graph_vertex_descriptor firstVertex(f_x_z_named_graph& graph,
                                ExtSetItPairVertexMap& extStatesPairVertexMap,
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

    // get iterator to this state in the fstates set
    state_typeSetBool pair_f = extFStateSet.insert(f_1);
    // get iterator to this state in the zstates set
    state_typeSetBool pair_z = extZStateSet.insert(z_state);

    // vertex is identified by ordered pair of these iterators
    // check its there in the map or make it if not
    FxzNamedGraphVertexBoolPair addVertex = vertexAddOrFind(graph,
                                    extStatesPairVertexMap,
                                    pair_f.first, pair_z.first);

    return addVertex.first;
}


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
                        fxz_named_graph_vertex_descriptor v,
                        const size_t epoch,
                        bool pathTravelledBefore,
                        const string addColour)
{
    bool retValue = false;

    // make sure this vertex is marked as visited
    fxz_named_vertexcolour_map_t vertexColourProp = get(vertex_colour, graph);
    fxz_named_edgecolour_map_t edgeColourProp = get(edge_colour, graph);

    vertexColourProp[v] = addColour;


    // now make sure all the next states are there
    // (add them if not there already)
    bool expanded = false;

    // if the vertex has some out edges,it will have all out-edges, so check this
    fxz_named_degree_size_type outs = out_degree(v, graph);
    if (outs > 0) { // v has out edges

        expanded = true;
    }
    else { // try to expand

        expanded = expandGraphKingmanForward(graph, extStatesPairVertexMap,
                        extFStateSet, extZStateSet, nsam, v, epoch);
    }

    if (expanded) {

        // choose which next state to go to from our current vertex addVertex.first
        FxzNamedVertexChoiceTuple choiceTuple = chooseNextVertex(graph, rgsl,
                                            v, pathTravelledBefore);
        // (nextVertex, proposal transition probability,
        //          kc transition probability, haveVertex, pathStillTravelledBefore)

        bool haveVertex = boost::get<3>(choiceTuple); // false if no vertex could be found

        if (haveVertex) { // if we found a vertex

            fxz_named_graph_vertex_descriptor nextVertex = get<0>(choiceTuple);
            bool pathStillTravelledBefore = boost::get<4>(choiceTuple);

            // fold the transition proposal probability into fsProposalProb
            fsProposalProb *= get<1>(choiceTuple);

            // fold the transition kc probability into fsKCProb
            fsKCProb *= get<2>(choiceTuple);

            // we will need the fstate property map and zstate property map
            fxz_named_fstateIt_map_t fstateProp = get(vertex_fstateSetIt, graph);
            fxz_named_zstateIt_map_t zstateProp = get(vertex_zstateSetIt, graph);

            // increment the epoch counter
            size_t newEpoch = epoch + 1;

            // use the internal property maps to get the new f state, z state
            state_typeSetIt new_f_stateSetIt = fstateProp[nextVertex];
            state_typeSetIt new_z_stateSetIt = zstateProp[nextVertex];

            // add the new f_state to the f_seq_container fs
            addStateToF_Sequence(fs, new_f_stateSetIt, newEpoch, nsam);

            // add the new z_state to the z_seq_container zs
            addStateToZ_Sequence(zs, new_z_stateSetIt, newEpoch, nsam);

            // make the next move
            retValue = nextFZSequenceState(graph, extStatesPairVertexMap,
                                            extFStateSet, extZStateSet, rgsl,
                                            fs, zs, fsProposalProb, fsKCProb,
                                            nsam, nextVertex,
                                            newEpoch, pathStillTravelledBefore,
                                            addColour);

            vertexColourProp[nextVertex] = addColour;

            FxzNamedGraphEdgeBoolPair p = edge(v, nextVertex, graph);
            edgeColourProp[p.first] = addColour;


        }
        else { // did not find a vertex to go to
            std::cout << "Problem finding f-sequence ";
            std::cout << "in transition from epoch " << epoch << std::endl;
        }
    }

    else { // graph was not expanded
    // check that is is because we are in the final state

        fxz_named_fstateIt_map_t fstateProp = get(vertex_fstateSetIt, graph);
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
bool allNextFZSequenceStates(f_x_z_named_graph& graph,
                            ExtSetItPairVertexMap& extStatesPairVertexMap,
                            state_typeSet& extFStateSet,
                            state_typeSet& extZStateSet,
                            const size_t nsam,
                            fxz_named_graph_vertex_descriptor v,
                            const size_t epoch)
{

    bool retValue = false;

    // now make sure all the next states are there
    // (add them if not there already)
    bool expanded = false;

    // if the vertex has some out edges,it will have all out-edges, so check this
    fxz_named_degree_size_type outs = out_degree(v, graph);
    if (outs > 0)  { // already expanded - nothing here for us to do
        retValue = true;
    }
    else { // v has no out edges yet

        expanded = expandGraphKingmanForward(graph, extStatesPairVertexMap,
                        extFStateSet, extZStateSet, nsam, v, epoch);

        if (expanded) {

            retValue = true;

            // increment the epoch counter
            size_t newEpoch = epoch + 1;

            // iterate through the next states from this vertex and expand them
            // but only if we have not now reached the final state
            if (newEpoch < nsam) {

                // get the targets of the out edges of this vertex v
                // seem to have to do this rather that directly using the out-edge iterator
                // which seems to get muddled if more edges are added in the recursion
                FxzNamedVecVertices targets;
                FxzNamedVecVerticesIt vit;

                for (pair<fxz_named_graph_out_edge_it, fxz_named_graph_out_edge_it> ed
                                    = out_edges(v, graph);
                                ed.first != ed.second; ++ed.first) {

                    targets.push_back(target(*(ed.first), graph));

                }
                // now go through the targets in turn
                for (vit = targets.begin(); vit < targets.end(); vit++) {

                fxz_named_graph_vertex_descriptor nextVertex = *vit;

                    // use allNextFZStates recursively, using retValue to ensure that
                    // we always get to the final state
                    retValue = retValue
                                && (allNextFZSequenceStates(graph,
                                                    extStatesPairVertexMap,
                                                    extFStateSet,
                                                    extZStateSet,
                                                    nsam,
                                                    nextVertex,
                                                    newEpoch));

                    if (!retValue) {
                        std::cout << "Cannot get to final state " << std::endl;
                        std::cout << "epoch is " << newEpoch << std::endl;

                        // get the internal property maps
                        break; // break out of for loop if failed to reach final state
                    }

                }
            } // end second check on epoch

        }
        else { // graph was not expanded
            // check that is is because we are in the final state

            fxz_named_fstateIt_map_t fstateProp = get(vertex_fstateSetIt, graph);
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
bool expandGraphKingmanForward(f_x_z_named_graph& graph,
                                ExtSetItPairVertexMap& extStatesPairVertexMap,
                                state_typeSet& extFStateSet,
                                state_typeSet& extZStateSet,
                                const size_t nsam,
                                fxz_named_graph_vertex_descriptor v,
                                const size_t epoch)
{
    bool retValue = false;

    // use the internal property map to get the current f state iterator
    fxz_named_fstateIt_map_t fstateProp = get(vertex_fstateSetIt, graph);
    state_typeSetIt f_stateSetIt = fstateProp[v];
    // and hence the f_state
    state_type f_state = *f_stateSetIt;

    // if the current f_state is the final state, we can't do anything
    // and transition probabilities depend on the state we are entering, and
    // we can check the current epoch against the current state
    //if ((f_state.sum() == epoch) && !checkFinalState(f_state, nsam)) {
    if ((epoch < nsam) && !checkFinalState(f_state, nsam)) {

        // use the internal property map to get the current z state
        fxz_named_zstateIt_map_t zstateProp = get(vertex_zstateSetIt, graph);
        state_typeSetIt z_stateSetIt = zstateProp[v];
        state_type z_state = *z_stateSetIt;

        size_t newEpoch = epoch + 1;

        // set up some variables to be filled by next subroutine
        VecSetIt new_f_statesSetIt;
        VecSetIt new_z_statesSetIt;
        VecProbs newProbs;
        double accProbs = 0.0; // use to accumulate probabilities of all transitions
        double rejectedProbs = 0.0; // use to accumulate rejected transitions

        retValue = findNewStates(extFStateSet, extZStateSet, f_state, z_state,
                                    nsam, newEpoch, accProbs, rejectedProbs,
                                    new_f_statesSetIt, new_z_statesSetIt, newProbs);
        // go through the new f and z states and add to graph
        if (retValue && !new_f_statesSetIt.empty()) {

            VecSetIt_it fit;
            VecSetIt_it zit;
            VecProbsIt dit;


            for (fit = new_f_statesSetIt.begin(), zit = new_z_statesSetIt.begin(),
                    dit = newProbs.begin();
                fit < new_f_statesSetIt.end(), zit < new_z_statesSetIt.end(),
                    dit < newProbs.end();
                fit++, zit++, dit++) {

                // make sure new state is in the graph
                FxzNamedGraphVertexBoolPair addVertex = vertexAddOrFind(graph,
                                                extStatesPairVertexMap,
                                                *fit,
                                                *zit);

                // make sure edge is in the graph
                // with reweighted probability and pure kc probability
                FxzNamedGraphEdgeBoolPair addEdge = edgeAddOrFind(graph,
                                v, addVertex.first, ((*dit)/accProbs), (*dit));

            }
        }

    } // end if not final state

    return retValue;
}




// choose an edge or next vertex to go to from our current vertex v
// tuple returned is
// (vertex, proposal transition probability, kc transition probability, bool, bool)
FxzNamedVertexChoiceTuple chooseNextVertex(f_x_z_named_graph& graph, gsl_rng* rgsl,
                        const fxz_named_graph_vertex_descriptor v, bool pathTravelledBefore)
{
    // we will need the edge proposal probability property map, kc probability map, and edge followed map
    fxz_named_proposal_prob_map_t edgeProposalProbProp = get(edge_kcrwprob, graph);
    fxz_named_kc_prob_map_t edgeKCProbProp = get(edge_kcprob, graph);
    fxz_named_followed_map_t edgeFollowedProp = get(edge_followed, graph);

    double accProbs = 0; // variable for accumulated probabilities of edges

    bool haveVertex = false;
    bool availablePaths = false;
    bool pathStillTravelledBefore = pathTravelledBefore;
    fxz_named_graph_vertex_descriptor nextVertex = v; // initialise to this vertex

    double transProposalProb = 0.0;  // proposal transition probability
    double transKCProb = 0.0; // pure kc transition probability

    // if we are on a known path we prefer to choose new directions
    if (pathTravelledBefore) {

        fxz_named_graph_edge_descriptor chosenEdge;
        FxzNamedVecEdges availableEdges; // container for available edges
        double availableProb = 0.0;

        //debug
        //std::cout << "Making restricted choices since on known path" << std::endl;

        // go through the out-edges of this vertex and add up the proposal transition probabilities
        // of the available edges
        for (pair<fxz_named_graph_out_edge_it, fxz_named_graph_out_edge_it> ed
                            = out_edges(v, graph);
                        ed.first != ed.second; ++ed.first) {
            bool followed = get(edgeFollowedProp, *(ed.first)); // has edge been followed?
            if (!followed) {
                availableProb += get(edgeProposalProbProp, *(ed.first));
                availableEdges.push_back(*(ed.first));
            }
            //debug
            else {
                //std::cout << "Edge " << *(ed.first)
                        //<< " has already been followed so excluding this edge probability "
                        //<< edgeProbProp[*(ed.first)] << std::endl;
            }
        }

        // did we find any available edges?
        if (!availableEdges.empty()) {

            availablePaths = true;

            if (availableEdges.size() == 1) { // only one choice available
                haveVertex = true;
                chosenEdge = availableEdges[0];
                edgeFollowedProp[chosenEdge] = true;
                nextVertex = target(chosenEdge, graph);
                transProposalProb = availableProb;
                transKCProb = edgeKCProbProp[chosenEdge];
                pathStillTravelledBefore = false; // now on unknown path


            }

            else { // have to make a choice

                // choose a random number between 0 and 1
                // use this to choose which vertex to go to
                double rand = gsl_rng_uniform(rgsl);

                //debug
                //std::cout << "Rand is " << rand << endl;

                // scale it by the total probability of the available paths
                rand *= availableProb;

                //debug
                //std::cout << "Scaled rand is " << rand << endl;

                // go through the available out-edges and choose one
                for (FxzNamedVecEdgesIt eit = availableEdges.begin();
                                        eit < availableEdges.end(); eit++) {

                    double prob = edgeProposalProbProp[*(eit)];

                    accProbs += prob;

                    if (rand < accProbs) {
                        chosenEdge = *(eit);
                        edgeFollowedProp[chosenEdge] = true;
                        nextVertex = target(chosenEdge, graph);
                        haveVertex = true;
                        transProposalProb = prob;
                        transKCProb = edgeKCProbProp[chosenEdge];
                        pathStillTravelledBefore = false; // now on unknown path


                        break; // found a target - break out of for loop

                        //debug
                        //string n1 = get(nameProp, v);
                        //string n2 = get(nameProp, nextVertex);
                        //std::cout << "Chosen edge from " << n1 << " to "
                                //<< n2 << " with proby " << prob << std::endl;

                    }
                } // end for loop through edges if more than one edge available
            }

        }
        //debug
        else {
            //std::cout << "No untravelled paths available" << std::endl;
        }
    }

    // if we are already on an untravelled path or
    // if we are on a travelled path but did not have any untravelled choices
    // we make a free choice from the edges out of this vertex
    if (!pathTravelledBefore || (pathTravelledBefore && !availablePaths)) {

        //debug
        //std::cout << "Making free choices" << std::endl;

        // choose a random number between 0 and 1
        // use this to choose which vertex to go to
        double rand = gsl_rng_uniform(rgsl);

        //debug
        //std::cout << "Rand is " << rand << endl;

        // go through the out-edges of this vertex
        for (pair<fxz_named_graph_out_edge_it, fxz_named_graph_out_edge_it> ed
                            = out_edges(v, graph);
                        ed.first != ed.second; ++ed.first) {
            double prob = edgeProposalProbProp[*(ed.first)];

            accProbs += prob;

            if (rand < accProbs) {
                fxz_named_graph_edge_descriptor chosenEdge = *(ed.first);
                nextVertex = target(chosenEdge, graph);
                haveVertex = true;
                transProposalProb = prob;
                transKCProb = edgeKCProbProp[chosenEdge];
                // don't mark the edge chosen as followed
                // either it is marked already or we got to it by an untravelled path
                pathStillTravelledBefore = (pathTravelledBefore && true);
                break; // found a target - break out of for loop

            }
        }
    }

    // if we could not find a vertex or if were no available vertices, nextVertex = v
    FxzNamedVertexChoiceTuple choiceTuple(nextVertex, transProposalProb,
                                    transKCProb, haveVertex, pathStillTravelledBefore);

    return choiceTuple;
}


// add a new vertex or find the vertex description if its already there
// return (fxz_named_graph_vertex_descriptor, bool)
// vertex descriptor is the descriptor for the vertex,
// bool is true if new vertex was inserted, false if it was already there
FxzNamedGraphVertexBoolPair vertexAddOrFind(f_x_z_named_graph& graph,
                                ExtSetItPairVertexMap& extStatesPairVertexMap,
                                const state_typeSetIt& f_stateSetIt,
                                const state_typeSetIt& z_stateSetIt)
{
    bool inserted = false;

    fxz_named_graph_vertex_descriptor v;
    state_typeSetItPair p(f_stateSetIt, z_stateSetIt);

    if (extStatesPairVertexMap.find(p) == extStatesPairVertexMap.end()) {

        // not found so need to add to graph and maps
        v = add_vertex(graph);
        // get the name label to apply to this vertex in the graph
        string label = statesToString(*f_stateSetIt, *z_stateSetIt);

        extStatesPairVertexMap[p] = v;
        fxz_named_name_map_t nameProp = get(vertex_name, graph);
        fxz_named_fstateIt_map_t fstateProp = get(vertex_fstateSetIt, graph);
        fxz_named_zstateIt_map_t zstateProp = get(vertex_zstateSetIt, graph);
        fxz_named_vertexcolour_map_t vertexColourProp = get(vertex_colour, graph);

        nameProp[v] = label;
        fstateProp[v] = f_stateSetIt;
        zstateProp[v] = z_stateSetIt;
        inserted = true;

        // default colour
        vertexColourProp[v] = "black";

        }

    else { // already there, just get the description from the external pair map
        v = extStatesPairVertexMap[p];

    }

    FxzNamedGraphVertexBoolPair retPair(v, inserted);

    return retPair;
}



// add a new edge or find the edge description if its already there
// the probability is added as a property of the edge
// and also the edge weight property becomes ln(1/probability)
// but if probability = 0, the edge weight = MAX_WEIGHT (defined in ebc_graphtypes.hpp)
// return (fxz_named_graph_edge_descriptor, bool)
// edge descriptor is the descriptor for the edge,
// bool is true if new edge was inserted, false if it was already there
FxzNamedGraphEdgeBoolPair edgeAddOrFind(f_x_z_named_graph& graph,
                                    const fxz_named_graph_vertex_descriptor v1,
                                    const fxz_named_graph_vertex_descriptor v2,
                                    const double edgeRWprob,
                                    const double edgeKCprob)
{
    FxzNamedGraphEdgeBoolPair p = add_edge(v1, v2, graph);

    fxz_named_kc_prob_map_t edgeKCProbProp = get(edge_kcprob, graph);


    if (p.second) { // new edge was inserted

        fxz_named_proposal_prob_map_t edgeProposalProbProp = get(edge_kcrwprob, graph);
        fxz_named_weight_map_t edgeWeightProp = get(edge_weight, graph);
        fxz_named_followed_map_t edgeFollowedProp = get(edge_followed, graph);
        fxz_named_edgecolour_map_t edgeColourProp = get(edge_colour, graph);

        edgeProposalProbProp[p.first] = edgeRWprob;
        edgeKCProbProp[p.first] = edgeKCprob;

        // use the proposal probabilities to calculate edge weights
        if (edgeRWprob > 0.0) {
            edgeWeightProp[p.first] = log(1.0/edgeRWprob); // natural log
        }

        // weight is log of 1/probability except if probability == 0.0
        else  { // MAX_WEIGHT defined in ebc_graphtypes.hpp
            edgeWeightProp[p.first] = 1.0*MAX_WEIGHT;
        }

        // also add the edge followed property, which is false by default
        edgeFollowedProp[p.first] = false;

        //default colour
        edgeColourProp[p.first] = "black";

    }
    else {
        double oldProb = get(edgeKCProbProp, p.first);

        if (oldProb != edgeKCprob)
            std::cout << "existing edge kc probability is  " << oldProb
                    << " whereas new kc probability is " << edgeKCprob << std::endl;
        //else
          //  std::cout << "edge weights match =   " << oldProb << std::endl;
    }

    return p;
}



// print description of the graph
void print_graph_desc(f_x_z_named_graph& graph)
{
    // get the internal property maps
    fxz_named_name_map_t nameProp = get(vertex_name, graph);
    fxz_named_fstateIt_map_t fstateProp = get(vertex_fstateSetIt, graph);
    fxz_named_zstateIt_map_t zstateProp = get(vertex_zstateSetIt, graph);

/*
    std::cout << "vertices(graph) = \n";

    for (pair<fxz_named_graph_vertex_it, fxz_named_graph_vertex_it> p = vertices(graph);
                        p.first != p.second; ++p.first) {

        string n = get(nameProp, *(p.first));
        std::cout << n << std::endl;

        state_type f_state = get(fstateProp, *(p.first));

        cout << "Printing f state from internal property map \n";
            copy (f_state.begin(), f_state.end(), ostream_iterator<size_t>(cout, "\t"));
            cout << "\n";


    }
*/
    fxz_named_proposal_prob_map_t edgeProposalProbProp = get(edge_kcrwprob, graph);
    fxz_named_kc_prob_map_t edgeKCProbProp = get(edge_kcprob, graph);

    std::cout << "edges(graph) = \n";

    for (pair<fxz_named_graph_edge_it, fxz_named_graph_edge_it> ed = edges(graph);
                        ed.first != ed.second; ++ed.first) {

        fxz_named_graph_vertex_descriptor src = source(*(ed.first), graph);
        fxz_named_graph_vertex_descriptor tgt = target(*(ed.first), graph);
        string n1 = nameProp[src];
        string n2 = nameProp[tgt];

        double proposalProb = edgeProposalProbProp[*(ed.first)];
        double kcProb = edgeKCProbProp[*(ed.first)];

        std::cout << "Edge from " << n1 << " to " << n2
                << " has proposal transition probability " << proposalProb
                << " and pure kc transition probability " << kcProb << std::endl;

    }
}


void makeDotGraph(f_x_z_named_graph& graph, string filename)
{
    /*
    fxz_named_name_map_t nameProp = get(vertex_name, graph);
    fxz_named_proposal_prob_map_t edgeProbProp = get(edge_kcrwprob, graph);
    string filename = "myGraph.dot";
    ofstream os(filename.c_str());
    if (os.is_open()) {
        write_graphviz(os, graph, make_label_writer(nameProp),
                        make_label_writer(edgeProbProp));
        os.close();
        std::cout << ".dot file in "
            << filename << std::endl << std::endl;

        }
    else {
        std::cout << "Error: could not open file named "
            << filename << std::endl << std::endl;
    }
    */

    fxz_named_name_map_t nameProp = get(vertex_name, graph);
    fxz_named_proposal_prob_map_t edgeProbProp = get(edge_kcrwprob, graph);
    fxz_named_vertexcolour_map_t vertexColourProp = get(vertex_colour, graph);
    fxz_named_edgecolour_map_t edgeColourProp = get(edge_colour, graph);
    ofstream os(filename.c_str());
    if (os.is_open()) {
        write_graphviz(os, graph, make_coloured_label_writer(nameProp, vertexColourProp),
                        make_coloured_label_writer(edgeProbProp, edgeColourProp));
        os.close();
        std::cout << ".dot file in "
            << filename << std::endl << std::endl;

        }
    else {
        std::cout << "Error: could not open file named "
            << filename << std::endl << std::endl;
    }
}





//dijkstra shortest path
void dijShortestPath(f_x_z_named_graph& graph,
    ExtSetItPairVertexMap& extStatesPairVertexMap,
    state_typeSet& extFStateSet,
    state_typeSet& extZStateSet,
    const state_type& z_state,
    const size_t nsam,
    const string addColour)
{
    fxz_named_graph_vertex_descriptor sourceVertex;
    fxz_named_graph_vertex_descriptor sinkVertex;
    //find the source and sink vertices
    boost::tie(sourceVertex, sinkVertex) =
                    findSourceAndSink(graph, extStatesPairVertexMap,
                    extFStateSet, extZStateSet, z_state, nsam);

    if (sourceVertex != sinkVertex) {

        // edge property maps
        fxz_named_weight_map_t edgeWeightProp = get(edge_weight, graph);
        fxz_named_proposal_prob_map_t edgeKCProbProp = get(edge_kcrwprob, graph);

        FxzEfficientVecVertices shortestPath; // container for the shortest path
                                                    // filled in call to dij(..)

        // call the generalised version of the function to get
        // the highest probability and fill in the shortestPath container
        double proby = dij(graph, edgeWeightProp, edgeKCProbProp,
            shortestPath, sourceVertex, sinkVertex, nsam, addColour);

        std::cout << "The first highest probability route to final epoch found is " << std::endl;

        // print the path
        printPath(graph, shortestPath);

        cout << "\nThe highest route probability is " << proby << "\n\n" << std::endl;

        // colour the path
        colourPath(graph, shortestPath, addColour);


    }
    else {
        std::cout << "Sorry, source and sink vertices are equal";
        std::cout << ": no path produced." << std::endl;
    }

}


//dijkstra longest path
void dijLongestPath(f_x_z_named_graph& graph,
    ExtSetItPairVertexMap& extStatesPairVertexMap,
    state_typeSet& extFStateSet,
    state_typeSet& extZStateSet,
    const state_type& z_state,
    const size_t nsam,
    const string addColour)
{
    fxz_named_graph_vertex_descriptor sourceVertex;
    fxz_named_graph_vertex_descriptor sinkVertex;
    //find the source and sink vertices
    boost::tie(sourceVertex, sinkVertex) =
                    findSourceAndSink(graph, extStatesPairVertexMap,
                    extFStateSet, extZStateSet, z_state, nsam);

    if (sourceVertex != sinkVertex) {

        // edge property map
        fxz_named_proposal_prob_map_t edgeKCProbProp = get(edge_kcrwprob, graph);

        FxzEfficientVecVertices longestPath; // container for the longest path
                                                    // filled in call to dij(..)

        // call the generalised version of the function
        // weights are given as probabilities
        // ie low probability = low weight
        double proby = dij(graph, edgeKCProbProp, edgeKCProbProp,
            longestPath, sourceVertex, sinkVertex, nsam, addColour);

        std::cout << "The first lowest probability route to final epoch found is " << std::endl;

        // print the path
        printPath(graph, longestPath);

        cout << "\nThe lowest route probability is " << proby << "\n\n" << std::endl;

        // colour the path
        colourPath(graph, longestPath, addColour);


    }
    else {
        std::cout << "Sorry, source and sink vertices are equal";
        std::cout << ": no path produced." << std::endl;
    }
}


// find average out edges at each epoch using breadth first search
void mapAverageOuts(f_x_z_named_graph& graph,
    ExtSetItPairVertexMap& extStatesPairVertexMap,
    state_typeSet& extFStateSet,
    state_typeSet& extZStateSet,
    const state_type& z_state,
    const size_t nsam)
{

    fxz_named_graph_vertex_descriptor sourceVertex;
    fxz_named_graph_vertex_descriptor sinkVertex;
    //find the source and sink vertices
    boost::tie(sourceVertex, sinkVertex) =
                    findSourceAndSink(graph, extStatesPairVertexMap,
                    extFStateSet, extZStateSet, z_state, nsam);

    if (sourceVertex != sinkVertex) {

        // call the generalised version of the function to get
        // the highest probability and fill in the shortestPath container
        mapAverageEdges(graph, sourceVertex);
    }
    else {
        std::cout << "Sorry, source and sink vertices are equal";
        std::cout << ": no map constructed." << std::endl;
    }

}



// find the source and sink vertices in a graph
// assumes that the graph is complete, ie the vertices exist
pair<fxz_named_graph_vertex_descriptor, fxz_named_graph_vertex_descriptor>
    findSourceAndSink(f_x_z_named_graph& graph,
    ExtSetItPairVertexMap& extStatesPairVertexMap,
    state_typeSet& extFStateSet,
    state_typeSet& extZStateSet,
    const state_type& z_state,
    const size_t nsam)
{
    // by default, set source and sink to first vertex in map
    pair<fxz_named_graph_vertex_it, fxz_named_graph_vertex_it> p = vertices(graph);
    fxz_named_graph_vertex_descriptor sourceVertex = *(p.first);;
    fxz_named_graph_vertex_descriptor sinkVertex = sourceVertex;

    // the vertex for EPOCH 1 is the pair (f_1 = (0 ... 0 1), z_state)
    size_t initVal = 0;
    state_type f_source(nsam, initVal); // nsam elements, all 0
    f_source[nsam-1] = 1; // set nsam-th element to 1

    state_type f_sink(nsam, initVal); // nsam elements, all 0
    f_sink[0] = nsam; // set 1st element to nsam

    state_type z_sink(nsam-1, initVal); // nsam-1 elements, all 0

    // get iterator to these states in the sets
    state_typeSetBool pair_f_source = extFStateSet.insert(f_source);
    state_typeSetBool pair_z_source = extZStateSet.insert(z_state);
    state_typeSetBool pair_f_sink = extFStateSet.insert(f_sink);
    state_typeSetBool pair_z_sink = extZStateSet.insert(z_sink);

    // find what vertices these are in the graph


    if (pair_f_source.second || pair_z_source.second
                || pair_f_sink.second || pair_z_sink.second) {

        // at least one of them had to be inserted
        // so leave both source and sink set to the first vertex in the graph

    }
    else {

        state_typeSetItPair sourcePair(pair_f_source.first, pair_z_source.first);
        state_typeSetItPair sinkPair(pair_f_sink.first, pair_z_sink.first);

        ExtSetItPairVertexMapIt source_it = extStatesPairVertexMap.find(sourcePair);
        ExtSetItPairVertexMapIt sink_it = extStatesPairVertexMap.find(sinkPair);

        if (source_it != extStatesPairVertexMap.end() &&
            sink_it != extStatesPairVertexMap.end()) { // both found

            sourceVertex = source_it->second;
            sinkVertex = sink_it->second;
        }
        // again if not both source and sink stay set to first vertex
    }
    return make_pair(sourceVertex, sinkVertex);

}


// colour a particular path in a graph
// could be templatised if we want to do this for other graph types
// but any graph type used must have colour maps for vertices and edges
void colourPath(f_x_z_named_graph& graph,
                FxzEfficientVecVertices& path,
                const string addColour)
{

    fxz_named_vertexcolour_map_t vertexColourProp = get(vertex_colour, graph);
    fxz_named_edgecolour_map_t edgeColourProp = get(edge_colour, graph);

    FxzEfficientVecVerticesIt vit1 = path.begin();
    vertexColourProp[*vit1] = addColour;

    // colour the edges and vertices in the shortest path
    if (path.size() > 1) {
        FxzEfficientVecVerticesIt vit2 = path.begin();
        vit2++; // second element

        for ( ; vit2 < path.end(); vit2++) {
            vertexColourProp[*vit2] = addColour;
            FxzNamedGraphEdgeBoolPair pp = edge(*vit1, *vit2, graph);
            edgeColourProp[pp.first] = addColour;
            vit1 = vit2;
        }
    }

}

// print a path
void printPath(f_x_z_named_graph& graph,
                FxzEfficientVecVertices& path)
{

    fxz_named_name_map_t nameProp = get(vertex_name, graph);

    FxzEfficientVecVerticesIt vit;

    for (vit = path.begin() ; vit < path.end(); vit++) {
        std::cout << nameProp[*vit] << "\n";
    }
    std::cout<<"\n" << std::endl;

}


