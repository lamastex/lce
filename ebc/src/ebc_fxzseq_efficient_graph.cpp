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
/*! \file ebc_fxzseq_efficient_graph.cpp
    \brief BGL methods for unvintaged (sized) n-coalescent and associated shape statistics with vertices as iterators to fstate and zstate sets, ie graph over product space
*/

#include <sstream>  // to be able to manipulate strings as streams

//#include <ebc_output.hpp>
//#include <ebc_params.hpp>

#include <ebc_fxzseq_efficient_graph.hpp>
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


/*! generate StateSequences of Kingman's unlabeled n-coalescent (i.e. the unvintaged and sized n-coalescent) forward in time

Using a BGL graph
*/

// generate a controlled f-sequence, using the z state from some data
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
                                const bool efficientFind)
{
    if(nsam < 2) {
        std::cerr << "\nnsam must be > 1 in GenerateStateSequence\n";
        exit(0);
    }

    if(nsam > 70) {
        std::cerr << "\nnsam should be <= 70 for this graph type, or you'll use too much RAM\n";
        std::cerr << "Try using the more efficient f-space only graph instead\n" << std::endl;
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
    fxz_efficient_graph_vertex_descriptor firstV = firstVertex(graph,
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
                                fs, fsProposalProb, fsKCProb,
                                nsam, firstV, epoch, efficientFind);

    //debug
    if (completed) {
        //std::cout << "\nf-sequence found is " << std::endl;
        //printStateSequence(fs);

        // std::cout << "printing map \n" << std::endl;
        //print_graph_desc(graph);
       // std::cout << "printing map \n" << std::endl;
    /*
        fxz_named_name_map_t nameProp = get(vertex_name, graph);
        std::cout << "printing usng print_vertices" << endl;
        print_vertices(graph, nameProp);
        std::cout << "printing usng print_edges" << endl;
        print_edges(graph, nameProp);
    */
        //print_graph_desc(graph);
    }
    else {
        std::cout << "failed to find StateSequence" << std::endl;
        std::cout << std::endl;
    }


    return completed;
}



bool generateFSequenceWithGraph(f_x_z_efficient_graph& graph,
                                ExtSetItPairVertexMap& extStatesPairVertexMap,
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

    return generateControlledFSequenceWithGraph(graph, extStatesPairVertexMap,
            extFStateSet, extZStateSet,
            rgsl, fs, z1_uncontrolled, fsProposalProb, fsKCProb,
            nsam, efficientFind);
}


// generate the complete graph with an initial controlling z_state
bool generateCompleteControlledGraph(f_x_z_efficient_graph& graph,
                                    ExtSetItPairVertexMap& extStatesPairVertexMap,
                                    state_typeSet& extFStateSet,
                                    state_typeSet& extZStateSet,
                                    const size_t nsam,
                                    const state_type& z_state)
{
    if(nsam < 2) {
        std::cerr << "\nnsam must be > 1 in GenerateStateSequence\n";
        exit(0);
    }

    if(nsam > 70) {
        std::cerr << "\nnsam should be <= 70 for this graph type, or you'll use too much RAM\n";
        std::cerr << "Try using the more efficient f-space only graph instead\n" << std::endl;
        exit(0);
    }

    size_t epoch = 1; // starting graph at epoch 1

    // make sure that the first (epoch 1) vertex is in the graph
    fxz_efficient_graph_vertex_descriptor firstV = firstVertex(graph,
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
        fxz_efficient_vertices_size_type numVert = num_vertices (graph);
        std::cout << "\n\nThe graph has " << numVert << " vertices\n" << std::endl;

*/

      // get the number of vertices
        fxz_efficient_vertices_size_type numVert = num_vertices (graph);
        std::cout << "\n\nThe graph for nsam " << nsam << " has " << numVert << " vertices\n" << std::endl;
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

    }
    else
        std::cout << "Could not successfully complete graph" << std::endl;

    return completed;
}


// generate the complete graph
// no initial controlling z_state
bool generateCompleteGraph(f_x_z_efficient_graph& graph,
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
fxz_efficient_graph_vertex_descriptor firstVertex(f_x_z_efficient_graph& graph,
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

    // vertex is identified by ordered pair of f and z states
    // check its there in the map or make it if not
    FxzEfficientGraphVertexBoolPair addVertex = vertexAddOrFind(graph,
                                    extStatesPairVertexMap,
                                    pair_f.first, pair_z.first);

    return addVertex.first;
}


bool nextFZSequenceState(f_x_z_efficient_graph& graph,
                        ExtSetItPairVertexMap& extStatesPairVertexMap,
                        state_typeSet& extFStateSet,
                        state_typeSet& extZStateSet,
                        gsl_rng* rgsl, StateSequence& fs,
                        double& fsProposalProb,
                        double& fsKCProb,
                        const size_t nsam,
                        const fxz_efficient_graph_vertex_descriptor v,
                        const size_t epoch,
                        bool pathTravelledBefore)
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
    else { // try to expand

        expanded = expandGraphKingmanForward(graph, extStatesPairVertexMap,
                        extFStateSet, extZStateSet, nsam, v, epoch);
    }

    if (expanded) {

        // choose which next state to go to from our current vertex addVertex.first

        FxzEfficientVertexChoiceTuple choiceTuple = chooseNextVertex(graph, rgsl,
                                            v, pathTravelledBefore);

        // (nextVertex, proposal transition probability,
        //        kc transition probability, haveVertex, pathStillTravelledBefore)

        bool haveVertex = boost::get<3>(choiceTuple); // false if no vertex could be found

        if (haveVertex) { // if we found a vertex

            fxz_efficient_graph_vertex_descriptor nextVertex = get<0>(choiceTuple);
            bool pathStillTravelledBefore = boost::get<4>(choiceTuple);

            // fold the transition probability into fsProposalProb
            fsProposalProb *= get<1>(choiceTuple);

            // fold the kc probability into fsKCProb
            fsKCProb *= get<2>(choiceTuple);

            // we will need the fstate property map and zstate property map
            fxz_efficient_fstateIt_map_t fstateProp = get(vertex_fstateSetIt, graph);
            fxz_efficient_zstateIt_map_t zstateProp = get(vertex_zstateSetIt, graph);

            // increment the epoch counter
            size_t newEpoch = epoch + 1;

            // use the internal property maps to get the new f state
            state_typeSetIt new_f_stateSetIt = fstateProp[nextVertex];

            // add the new f_state to the f_seq_container fs
            addStateToF_Sequence(fs, new_f_stateSetIt, newEpoch, nsam);

            // make the next move
            retValue = nextFZSequenceState(graph, extStatesPairVertexMap,
                                            extFStateSet, extZStateSet, rgsl,
                                            fs, fsProposalProb, fsKCProb,
                                            nsam, nextVertex,
                                            newEpoch, pathStillTravelledBefore);

        }
        else { // did not find a vertex to go to
            std::cout << "Problem finding f-sequence ";
            std::cout << "in transition from epoch " << epoch << std::endl;
        }
    }

    else { // graph was not expanded
    // check that is is because we are in the final state

        fxz_efficient_fstateIt_map_t fstateProp = get(vertex_fstateSetIt, graph);
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
bool allNextFZSequenceStates(f_x_z_efficient_graph& graph,
                            ExtSetItPairVertexMap& extStatesPairVertexMap,
                            state_typeSet& extFStateSet,
                            state_typeSet& extZStateSet,
                            const size_t nsam,
                            const fxz_efficient_graph_vertex_descriptor v,
                            const size_t epoch)
{

    bool retValue = false;

    // now make sure all the next states are there
    // (add them if not there already)
    bool expanded = false;

    // if the vertex has some out edges,it will have all out-edges, so check this
    fxz_efficient_degree_size_type outs = out_degree(v, graph);
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
                FxzEfficientVecVertices targets;
                FxzEfficientVecVerticesIt vit;

                for (pair<fxz_efficient_graph_out_edge_it, fxz_efficient_graph_out_edge_it> ed
                                    = out_edges(v, graph);
                                ed.first != ed.second; ++ed.first) {

                    targets.push_back(target(*(ed.first), graph));

                }
                // now go through the targets in turn
                for (vit = targets.begin(); vit < targets.end(); vit++) {

                fxz_efficient_graph_vertex_descriptor nextVertex = *vit;

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

            fxz_efficient_fstateIt_map_t fstateProp = get(vertex_fstateSetIt, graph);
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
bool expandGraphKingmanForward(f_x_z_efficient_graph& graph,
                                ExtSetItPairVertexMap& extStatesPairVertexMap,
                                state_typeSet& extFStateSet,
                                state_typeSet& extZStateSet,
                                const size_t nsam,
                                const fxz_efficient_graph_vertex_descriptor v,
                                const size_t epoch)
{
    bool retValue = false;

    // use the internal property map to get the current f state iterator
    fxz_efficient_fstateIt_map_t fstateProp = get(vertex_fstateSetIt, graph);
    state_typeSetIt f_stateSetIt = fstateProp[v];
    // and hence the f_state
    state_type f_state = *f_stateSetIt;

    // if the current f_state is the final state, we can't do anything
    // and transition probabilities depend on the state we are entering, and
    // we can check the current epoch against the current state
    //if ((f_state.sum() == epoch) && !checkFinalState(f_state, nsam)) {
    if ((epoch < nsam) && !checkFinalState(f_state, nsam)) {

        // use the internal property map to get the current z state
        fxz_efficient_zstateIt_map_t zstateProp = get(vertex_zstateSetIt, graph);
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
            FxzEfficientGraphVertexBoolPair addVertex = vertexAddOrFind(graph,
                                            extStatesPairVertexMap,
                                            *fit,
                                            *zit);

            // make sure edge is in the graph
            // with reweighted probability and pure kc probability
            FxzEfficientGraphEdgeBoolPair addEdge = edgeAddOrFind(graph,
                            v, addVertex.first, ((*dit)/accProbs), (*dit));

            }
        }

    } // end if not final state

    return retValue;
}


// choose an edge or next vertex to go to from our current vertex v
// return tuple (vertex, proposal transition probability, kc transition probability, bool, bool)
FxzEfficientVertexChoiceTuple chooseNextVertex(f_x_z_efficient_graph& graph, gsl_rng* rgsl,
                        const fxz_efficient_graph_vertex_descriptor v, bool pathTravelledBefore)
{
    // we will need the edge proposal probability property map, kc probability map and edge followed map
    fxz_efficient_proposal_prob_map_t edgeProposalProbProp = get(edge_kcrwprob, graph);
    fxz_efficient_kc_prob_map_t edgeKCProbProp = get(edge_kcprob, graph);
    fxz_efficient_followed_map_t edgeFollowedProp = get(edge_followed, graph);

    double accProbs = 0; // variable for accumulated probabilities of edges

    bool haveVertex = false;
    bool availablePaths = false;
    bool pathStillTravelledBefore = pathTravelledBefore;
    fxz_efficient_graph_vertex_descriptor nextVertex = v; // initialise to this vertex

    double transProposalProb = 0.0;
    double transKCProb = 0.0;

    // if we are on a known path we prefer to choose new directions
    if (pathTravelledBefore) {

        fxz_efficient_graph_edge_descriptor chosenEdge;
        FxzEfficientVecEdges availableEdges; // container for available edges
        double availableProb = 0.0;

        // go through the out-edges of this vertex and add up the probabilities
        // of the available edges
        for (pair<fxz_efficient_graph_out_edge_it, fxz_efficient_graph_out_edge_it> ed
                            = out_edges(v, graph);
                        ed.first != ed.second; ++ed.first) {
            bool followed = get(edgeFollowedProp, *(ed.first)); // has edge been followed?
            if (!followed) {
                availableProb += get(edgeProposalProbProp, *(ed.first));
                availableEdges.push_back(*(ed.first));
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

                // scale it by the total probability of the available paths
                rand *= availableProb;

                // go through the available out-edges and choose one
                for (FxzEfficientVecEdgesIt eit = availableEdges.begin();
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

    }

    // if we are already on an untravelled path or
    // if we are on a travelled path but did not have any untravelled choices
    // we make a free choice from the edges out of this vertex
    if (!pathTravelledBefore || (pathTravelledBefore && !availablePaths)) {

        // choose a random number between 0 and 1
        // use this to choose which vertex to go to
        double rand = gsl_rng_uniform(rgsl);

        // go through the out-edges of this vertex
        for (pair<fxz_efficient_graph_out_edge_it, fxz_efficient_graph_out_edge_it> ed
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
    FxzEfficientVertexChoiceTuple choiceTuple(nextVertex, transProposalProb,
                                    transKCProb, haveVertex, pathStillTravelledBefore);

    return choiceTuple;
}


// add a new vertex or find the vertex description if its already there
// return (fxz_efficient_graph_vertex_descriptor, bool)
// vertex descriptor is the descriptor for the vertex,
// bool is true if new vertex was inserted, false if it was already there
FxzEfficientGraphVertexBoolPair vertexAddOrFind(f_x_z_efficient_graph& graph,
                                                ExtSetItPairVertexMap& extStatesPairVertexMap,
                                                const state_typeSetIt& f_stateSetIt,
                                                const state_typeSetIt& z_stateSetIt)
{

    bool inserted = false;

    fxz_efficient_graph_vertex_descriptor v;
    state_typeSetItPair p(f_stateSetIt, z_stateSetIt);

    if (extStatesPairVertexMap.find(p) == extStatesPairVertexMap.end()) {

        // not found so need to add to graph and maps
        v = add_vertex(graph);
        extStatesPairVertexMap[p] = v;
        fxz_efficient_fstateIt_map_t fstateProp = get(vertex_fstateSetIt, graph);
        fxz_efficient_zstateIt_map_t zstateProp = get(vertex_zstateSetIt, graph);
        fstateProp[v] = f_stateSetIt;
        zstateProp[v] = z_stateSetIt;
        inserted = true;

    }

    else { // already there, just get the description from the external pair map
        v = extStatesPairVertexMap[p];
    }

    FxzEfficientGraphVertexBoolPair retPair(v, inserted);

    return retPair;
}



// add a new edge or find the edge description if its already there
// the probability is added as a property of the edge
// and also the edge weight property becomes ln(1/probability)
// but if probability = 0, the edge weight = MAX_WEIGHT (defined in ebc_graphtypes.hpp)
// return (fxz_efficient_graph_edge_descriptor, bool)
// edge descriptor is the descriptor for the edge,
// bool is true if new edge was inserted, false if it was already there
FxzEfficientGraphEdgeBoolPair edgeAddOrFind(f_x_z_efficient_graph& graph,
                                    const fxz_efficient_graph_vertex_descriptor v1,
                                    const fxz_efficient_graph_vertex_descriptor v2,
                                    const double edgeRWprob,
                                    const double edgeKCprob)
{
    fxz_efficient_kc_prob_map_t edgeKCProbProp = get(edge_kcprob, graph);

    FxzEfficientGraphEdgeBoolPair p = add_edge(v1, v2, graph);
    if (p.second) { // new edge was inserted

        fxz_efficient_graph_edge_descriptor ed;
        fxz_efficient_proposal_prob_map_t edgeProposalProbProp = get(edge_kcrwprob, graph);
        fxz_efficient_weight_map_t edgeWeightProp = get(edge_weight, graph);
        fxz_efficient_followed_map_t edgeFollowedProp = get(edge_followed, graph);

        edgeProposalProbProp[p.first] = edgeRWprob;
        edgeKCProbProp[p.first] = edgeKCprob;

        // use the proposal probabilities to calculate edge weights
        if (edgeRWprob > 0.0)
            edgeWeightProp[p.first] = log(1.0/edgeRWprob); // natural log

        // weight is log of 1/probability except if probability == 0.0
        else // MAX_WEIGHT defined in ebc_graphtypes.hpp
            edgeWeightProp[p.first] = 1.0*MAX_WEIGHT;

        // also add the edge followed property, which is false by default
        edgeFollowedProp[p.first] = false;
    }

    else {

        double oldProb = get(edgeKCProbProp, p.first);

        if (oldProb != edgeKCprob)
            std::cout << "existing edge kc probability is  " << oldProb
                    << " whereas new kc probability is " << edgeKCprob << std::endl;
        /*
        else
            std::cout << "edge weights match =   " << oldProb << std::endl;
            */
    }


    return p;
}


// print description of the graph
void print_graph_desc(f_x_z_efficient_graph& graph)
{
    // check the number of vertices
    fxz_efficient_vertices_size_type numVert = num_vertices (graph);

    if (numVert > 500) {
        std::cout << "\n\nThe graph has more than 500 vertices and will not be printed\n" << std::endl;
    }
    else {

        fxz_efficient_fstateIt_map_t fstateProp = get(vertex_fstateSetIt, graph);
        fxz_efficient_zstateIt_map_t zstateProp = get(vertex_zstateSetIt, graph);

        fxz_efficient_proposal_prob_map_t edgeProposalProbProp = get(edge_kcrwprob, graph);
        fxz_efficient_kc_prob_map_t edgeKCProbProp = get(edge_kcprob, graph);

        std::cout << "edges(graph) = \n";

        for (pair<fxz_efficient_graph_edge_it, fxz_efficient_graph_edge_it> ed = edges(graph);
                            ed.first != ed.second; ++ed.first) {

            fxz_efficient_graph_vertex_descriptor src = source(*(ed.first), graph);
            fxz_efficient_graph_vertex_descriptor tgt = target(*(ed.first), graph);
            state_typeSetIt src_f_stateSetIt = get(fstateProp, src);
            state_typeSetIt tgt_f_stateSetIt = get(fstateProp, tgt);
            state_typeSetIt src_z_stateSetIt = get(zstateProp, src);
            state_typeSetIt tgt_z_stateSetIt = get(zstateProp, tgt);
            state_type src_f_state = *src_f_stateSetIt;
            state_type tgt_f_state = *tgt_f_stateSetIt;
            state_type src_z_state = *src_z_stateSetIt;
            state_type tgt_z_state = *tgt_z_stateSetIt;
            string src_label = statesToString(src_f_state, src_z_state);
            string tgt_label = statesToString(tgt_f_state, tgt_z_state);

            double proposalProb = edgeProposalProbProp[*(ed.first)];
            double kcProb = edgeKCProbProp[*(ed.first)];

            std::cout << "Edge from " << src_label << " to " << tgt_label
                    << " has proposal transition proby "
                    << proposalProb << " and pure kc transition probability "
                    << kcProb <<std::endl;
        }

    }

}



