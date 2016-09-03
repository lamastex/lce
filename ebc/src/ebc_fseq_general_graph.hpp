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
/*! \file ebc_fseq_general_graph.hpp
    \brief prototypes of general methods associated with use of BGL graphs for unvintaged (sized) n-coalescent and associated shape statistics.
*/
#ifndef __FSEQGENERALGRAPH_HPP__
#define __FSEQGENERALGRAPH_HPP__


#include <ebc_sfstypes.hpp>
#include <ebc_graphtypes.hpp>

#include <gsl/gsl_rng.h>

#include <boost/graph/graphviz.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/breadth_first_search.hpp>

#include <iomanip> // manipulation of io formatting
#include <numeric> // to use accumulate
#include <utility> // to use pairs



/*! \brief Method to try to find an entry in a XStarSequencesMap for an xstar.

This method first looks for an entry in the map for this xstar.  If this is found, the returned pair is the iterator to the (key for xstar, value associated with that key) pair in the map and a bool = true indicating that the value (StateSequenceSetItSummaryMap map) associated with the xstar key is complete.

If the xstar is not already in the map, the method looks for an entry for an equivalent xstar.  An equivalent xstar is one which matches this xstar for all entries indexed [i] for i > 1 (ie, the entries at [0] and [1], corresponding to one and two leaves subtended, are ignored). If an equivalent xstar is found, the value (StateSequenceSetItSummaryMap map) for that equivalent xstar is copied and the copy becomes the value associated with the key for this xstar in the XStarSequencesMap map.  The countGenerated member for each FSeqInfo value in the new copy of the StateSequenceSetItSummaryMap map is set to 0 to reflect the fact that no additional f-sequences had to be generated to get this map.  The returned pair is the iterator to the (key for xstar, value associated with that key) pair in the map (where the value associated with the key is the copied map) and a bool = true indicating that the value (StateSequenceSetItSummaryMap map) associated with the xstar key is complete.

If no equivalent xstar is found, the returned pair is the iterator to the (key for xstar, value associated with that key) pair in the map and a bool = false indicating that the value (StateSequenceSetItSummaryMap map) associated with the xstar key is empty.

\param xstarBool is a variable passed by reference and updated by the method to give a pair (iterator to a (key, value) pair in the map, bool).  The key will correspond to the xstar and the value will be a map for that xstar.  The bool will indicate whether the map is already complete (true) or needs to be filled in (false).
\param xstarToSequencesInfoMap is the map from xstar keys to StateSequenceSetItSummaryMap values.
\param xstar is the xstar which we are trying to find or add as a key to the xstarToSequencesInfoMap map.
\post a key, value pair for the xstar will have been added to the map if that xstar did not already exist as a key in the map.  If the xstar did not exist but an equivalent xstar did, the xstar key will be paired with a value which is a copy of the StateSequenceSetItSummaryMap map for the equivalent xstar but with the countGenerated members for all FSeqInfo values set to 0 to show that no sequences had to be generated.
*/
XStarSequencesMapBool& findExistingMap(XStarSequencesMapBool& xstarBool,
                                        XStarSequencesMap& xstarToSequencesInfoMap,
                                        const state_type& xstar);




void zeroMapCounts(StateSequenceSetItSummaryMap& sssMap);



/*! \brief Method to find new accessible f-state and z-state from current (f, z) state.

This method provides a commmon method for all graph types to find the f-states which can follow from the current f-state given the current z-state.

The method takes the exising f-state and calculates the permitted f-states that can follow (i.e., at the next epoch, moving forwards from epoch 1, one common ancestor) this in the tree given the current z-state, and also the new z-state that would be associated with each new f-state.

This method updates containers and values passed by reference by the calling method.  The set iterators to the acceptable f-states and associated z-states, and the associated Kingman Coalescent forward chain transition probabilities, are stored in containers passed by reference to the method. The values of the total allowable Kingman Coalescent forward chain transition probabilities and the rejected (not allowable) Kingman Coalescent forward chain transition probabilities, both passed by reference to the method, are also updated.

\pre Sets for unique f=states and z-states.  Empty containers for set iterators to allowable new f-states and their associated z-states and transition probabilities.  Variables for total allowable Kingman Coalescent forward chain transition probabilities and not allowable Kingman Coalescent forward chain transition probabilities, initialised to zero.
\param extFStateSet is a reference to a set of unique f-states.  New f-states may be added to this container by the method. The set is also used to find iterators to existing f-states in the set.
\param extZStateSet is a reference to a set of unique z-states.  New z-states may be added to this container by the method. The set is also used to find iterators to existing z-states in the set.
\param f_state is the current f-state.
\param z_state is the current z-state.
\param nsam is the number of individual samples we are working with.
\param newEpoch is the epoch we are moving to, ie the epoch associated with the new states.
\param accProbs is a reference to a value for the total allowable (accepted) Kingman Coalescent forward chain transition probabilities, passed by reference and updated by the method.
\param rejectedProbs is a value for the total not allowable (rejected) Kingman Coalescent forward chain transition probabilities, passed by reference and updated by the method.
\param new_f_statesSetIt is a reference to a container to store iterators into the set of unique f-states.  At the end of the method the container will be filled with iterators representing f-states in the new epoch which can be reached from the current f-state given the current z-state.  The items in this container are effectively paired with the items in the new_z_statesSetIt container, and the newProbs container, ie items with the same index in the three different containers are associated.
\param new_z_statesSetIt is a reference to a container to store iterators into the set of unique z-states.  At the end of the method each item in the container  will representing the z-state in the new epoch associated with the f-state represented by the item in the new_f_statesSetIt container with the same index.
\param newProbs is a reference to a container to store Kingman Coalescent forward transition probabilities.  At the end of the method each item in the container  will be the Kingman Coalescent forward transition probability for the transition from the current f and z state to the f-state and z-state represented by the items in the new_f_statesSetIt and new_z_statesSetIt containers with the same index.
\post Updated sets for unique f=states and z-states.  Containers of set iterators to allowable new f-states and associated z-states and transition probabilities.  Values for total allowable Kingman Coalescent forward chain transition probabilities and not allowable Kingman Coalescent forward chain transition probabilities.
\return true if at least one f-state to move to was found, false otherwise.
*/

bool findNewStates(state_typeSet& extFStateSet,
            state_typeSet& extZStateSet,
            const state_type& f_state,
            const state_type& z_state,
            const size_t nsam,
            const size_t newEpoch,
            double& accProbs,
            double& rejectedProbs,
            VecSetIt& new_f_statesSetIt,
            VecSetIt& new_z_statesSetIt,
            VecProbs& newProbs);



void addStateToF_Sequence(StateSequence& fs,
                            const state_typeSetIt& f_stateIt,
                            const size_t epoch,
                            const size_t nsam);

void addStateToZ_Sequence(StateSequence& zs,
                            const state_typeSetIt& z_stateIt,
                            const size_t epoch,
                            const size_t nsam);

string statesToString(const state_type& f_state,
                        const state_type& z_state);


string stateToString(const state_type& some_state);


bool checkFinalState(const state_type& f_state,
                    const size_t nsam);




void printStateSequence(const StateSequence& sequence);

void printStateSequence(const StateSequenceSetIt& seqSetIt);


state_type& sfsToXstar(state_type& xstar, const sfs_array_type& sfs);

//returns lineageLenItons by reference
std::valarray<double>& epochTimesProdFseq(std::valarray<double>& lineageLenItons,
                                        const std::valarray<double> & epochTimes,
                                        const StateSequence& fs,
                                        const size_t nsam);

// templatised function to generate all fsequences for a given site frequency spectrum ObsSfs
// generate all FSequences (ie, total probability 1) for this SFS
template <typename Graph, typename Map>
pair <XStarSequencesMapIt, bool> generateAllFSeq(Graph& graph,
                                        Map& vertexMap,
                                        state_typeSet& extFStateSet,
                                        state_typeSet& extZStateSet,
                                        gsl_rng* rgsl,
                                        const sfs_array_type& obsSfs,
                                        XStarSequencesMap& xstarToSequencesInfoMap,
                                        StateSequenceSet& fsequenceSet,
                                        const size_t nsam)
{
    // find the xstar for this sfs
    state_type xstar;
    xstar = sfsToXstar(xstar, obsSfs);

    //debug
    //string xstar_label = stateToString(xstar);
    //cout << "\n\nxstar is " << xstar_label << endl;

    // check if we have got fsequences for this xstar already
    XStarSequencesMapBool xstarBool;

    bool gotAllFSequences = false;

    // try to find an existing map for this xstar or its equivalent
    xstarBool = findExistingMap(xstarBool, xstarToSequencesInfoMap, xstar);

    if(xstarBool.second) { // no sequences and no equivalent matches so far for this xstar

        double totalProposalProbs = 0.0;

        size_t iii = 0;

        //debug
        //std::cout << "new xstar - generating sequences " << endl;

        // xstarBool.first is an iterator pointing to a pair (xstar, fsequenceMap)
        // generate the sequences for this xstar and summarise them in the map
        //for(iii = 0; (iii< 100000000) && (totalProposalProbs < 0.999999999); iii++) {
        for(iii = 0; (iii < 1000000000) && (totalProposalProbs < 0.999999999); iii++) {

            StateSequence fsNew;
            double fsProposalProb = 0.0;
            double fsKCProb = 0.0;

            bool success = generateControlledFSequenceWithGraph(graph,
                                                                vertexMap,
                                                                extFStateSet,
                                                                extZStateSet,
                                                                rgsl,
                                                                fsNew,
                                                                ((xstarBool.first)->first),
                                                                fsProposalProb,
                                                                fsKCProb,
                                                                nsam);

            if (success) {

                // FseqInfo to be copied, already given a generated count and imp sample count of 1
                FseqInfo fseqInfo(1, 1, 0.0, 0.0);

                //try inserting the f_sequence into the set
                StateSequenceSetItBool fsNewItBool = fsequenceSet.insert(fsNew);
                // fsNewItBool.first will point to the iterator to the state sequence newFs in the set

                // now see if this iterator is already in the map for this particular xstar
                StateSequenceSetItSummaryMapBool sssMapBool
                        = (xstarBool.first)->second.insert(make_pair(fsNewItBool.first, fseqInfo));

                if(sssMapBool.second) { // new insertion into this xstar's map
                    // sssMapBool.first is an iterator pointing to new map element
                    // which is a pair<fseqSetIt, FseqInfo>
                    (sssMapBool.first)->second.proposalProb = fsProposalProb;    // update prop_prob
                    (sssMapBool.first)->second.kcProb = fsKCProb;    // update kc probability

                    //increment the total proposal prob of all f-seqns visited
                    totalProposalProbs += fsProposalProb;

                }

                else // already there so increment the generated and importance sample counts
                        // for this particular f-seqn

                    (sssMapBool.first)->second.countImpSample += 1;
                    (sssMapBool.first)->second.countGenerated += 1;

            } // end if successful generation of fsequence

            gotAllFSequences = !(totalProposalProbs < 0.999999999);

        } // end of loop trying to get fsequences covering most of the probability for this xstar

    } // end of section generating fsequences for a previously unseen xstar

    else { // xstar is already there so use the sequences currently there
        gotAllFSequences = true;
            // xstarBool.first is an iterator pointing to a pair (xstar, fsequenceMap)
        //debug
        //cout << xstar_label << " has already been done" << std::endl;
    }

    // xstarBool.first is an iterator pointing to a pair (xstar, fsequenceMap)
    // gotAllSequences says if we got all the sequence to make proposal probability add to 1

    return make_pair(xstarBool.first, gotAllFSequences);
}



// templatised function to generate an importance sample of FSequences for this SFS
template <typename Graph, typename Map>
pair <XStarSequencesMapIt, size_t> generateImportanceSample(Graph& graph,
                                        Map& vertexMap,
                                        state_typeSet& extFStateSet,
                                        state_typeSet& extZStateSet,
                                        gsl_rng* rgsl,
                                        const sfs_array_type& obsSfs,
                                        XStarSequencesMap& xstarToSequencesInfoMap,
                                        StateSequenceSet& fsequenceSet,
                                        const size_t nsam,
                                        const size_t importanceSampleSize)
{
    // find the xstar for this sfs
    state_type xstar;
    xstar = sfsToXstar(xstar, obsSfs);

    //debug
    //string xstar_label = stateToString(xstar);
    //cout << "\n\nxstar is " << xstar_label << endl;

    // check if we have got fsequences for this xstar already
    XStarSequencesMapBool xstarBool;

    // try to find an existing map for this xstar or its equivalent
    xstarBool = findExistingMap(xstarBool, xstarToSequencesInfoMap, xstar);

    size_t totalFSImpSample = 0;

    if(xstarBool.second) { // no equivalent sequence done so far for this xstar

        double totalProposalProbs = 0.0;

        // xstarBool.first is an iterator pointing to a pair (xstar, fsequenceMap)
        // generate the sequences for this xstar and summarise them in the map
        while (totalFSImpSample < importanceSampleSize) {

            StateSequence fsNew;
            double fsProposalProb = 0.0;
            double fsKCProb = 0.0;

            // make sure generation is not trying to use efficient generation of all fsequences
            bool success = generateControlledFSequenceWithGraph(graph, vertexMap,
                    extFStateSet, extZStateSet,
                    rgsl, fsNew, ((xstarBool.first)->first), fsProposalProb, fsKCProb, nsam, false);

            if (success) {

                // increment count of total f-sequences
                totalFSImpSample ++;

                FseqInfo fseqInfo(1, 1, 0.0, 0.0);

                //try inserting the f_sequence into the set
                StateSequenceSetItBool fsNewItBool = fsequenceSet.insert(fsNew);
                // fsNewItBool.first will point to the iterator to the state sequence newFs in the set

                // now see if this iterator is already in the map for this particular xstar
                StateSequenceSetItSummaryMapBool sssMapBool
                        = (xstarBool.first)->second.insert(make_pair(fsNewItBool.first, fseqInfo));

                if(sssMapBool.second) { // new insertion into this xstar's map
                    // sssMapBool.first is an iterator pointing to new map element
                    // which is a pair<fseqSetIt, FseqInfo>
                    (sssMapBool.first)->second.proposalProb = fsProposalProb;    // update prop_prob
                    (sssMapBool.first)->second.kcProb = fsKCProb;    // update kc probability

                    // increment the total proposal prob of all f-seqns visited
                    totalProposalProbs += fsProposalProb;
                }

                else // already there so increment the number of visits to this particular f-seqn

                    (sssMapBool.first)->second.countImpSample += 1;
                    (sssMapBool.first)->second.countGenerated += 1;

            } // end if successful generation of fsequence

        } // end of for loop trying to get a set number of fsequences for this xstar

    } // end of section generating fsequences for a previously unseen xstar

    else { // xstar is already there or equivalent found
            // so use some sequences currently there
            // xstarBool.first is an iterator pointing to a pair (xstar, fsequenceMap)

        // need to find out how many fsequences there are:
        StateSequenceSetItSummaryMapIt mapIt;

        for (mapIt = (xstarBool.first->second).begin();
                    mapIt != (xstarBool.first->second).end(); mapIt++) {

            totalFSImpSample += (mapIt->second.countImpSample);
        }

        //debug
        //cout << xstar_label << " has already been done" << std::endl;
    }

    // xstarBool.first is an iterator pointing to a pair (xstar, fsequenceMap)

    return make_pair(xstarBool.first, totalFSImpSample);
}




// templatised function to generate or add to an importance sample of FSequences for this SFS
// return tuple (iterator, importance sample size, number of f sequences generated)
template <typename Graph, typename Map>
boost::tuple <XStarSequencesMapIt, size_t, size_t> generateOrAddToImportanceSample(Graph& graph,
                                        Map& vertexMap,
                                        state_typeSet& extFStateSet,
                                        state_typeSet& extZStateSet,
                                        gsl_rng* rgsl,
                                        const sfs_array_type& obsSfs,
                                        XStarSequencesMap& xstarToSequencesInfoMap,
                                        StateSequenceSet& fsequenceSet,
                                        const size_t nsam,
                                        const size_t importanceSampleSize,
                                        const size_t isIncr)
{
    // find the xstar for this sfs
    state_type xstar;
    xstar = sfsToXstar(xstar, obsSfs);

    // map to be copied
    StateSequenceSetItSummaryMap fsequenceMap;

    //debug
    //string xstar_label = stateToString(xstar);
    //cout << "\n\nxstar is " << xstar_label << endl;

    // if importanceSampleSize is 0, find if we have already done this one

    // check if we have got fsequences for this xstar already
    XStarSequencesMapBool xstarBool;

    // try to find an existing map for this xstar or its equivalent
    xstarBool = findExistingMap(xstarBool, xstarToSequencesInfoMap, xstar);

    size_t totalFSGenerated = 0;
    size_t totalFSImpSample = 0;

    if(((importanceSampleSize == 0) && xstarBool.second) // no sequences so far for this xstar
            || ((importanceSampleSize > 0) && !xstarBool.second)) { // exists, adding to importance sample

        // xstarBool.first is an iterator pointing to a pair (xstar, fsequenceMap)
        // generate the sequences for this xstar and summarise them in the map
        while (totalFSImpSample < isIncr) {

            StateSequence fsNew;
            double fsProposalProb = 0.0;
            double fsKCProb = 0.0;

            // make sure generation is not trying to use efficient generation of all fsequences
            bool success = generateControlledFSequenceWithGraph(graph, vertexMap,
                    extFStateSet, extZStateSet,
                    rgsl, fsNew, ((xstarBool.first)->first), fsProposalProb, fsKCProb, nsam, false);

            if (success) {

                // increment count of total f-sequences generated and in importance sample
                totalFSGenerated ++;
                totalFSImpSample ++;

                FseqInfo fseqInfo(1, 1, 0.0, 0.0);

                //try inserting the f_sequence into the set
                StateSequenceSetItBool fsNewItBool = fsequenceSet.insert(fsNew);
                // fsNewItBool.first will point to the iterator to the state sequence newFs in the set

                // now see if this iterator is already in the map for this particular xstar
                StateSequenceSetItSummaryMapBool sssMapBool
                        = (xstarBool.first)->second.insert(make_pair(fsNewItBool.first, fseqInfo));

                if(sssMapBool.second) { // new insertion into this xstar's map
                    // sssMapBool.first is an iterator pointing to new map element
                    // which is a pair<fseqSetIt, FseqInfo>
                    (sssMapBool.first)->second.proposalProb = fsProposalProb;    // update prop_prob
                    (sssMapBool.first)->second.kcProb = fsKCProb;    // update kc probability
                }

                else // already there so increment the number of visits to this particular f-seqn

                    (sssMapBool.first)->second.countImpSample += 1;
                    (sssMapBool.first)->second.countGenerated += 1;

                } // end if successful generation of fsequence

        } // end of for loop trying to get a set number of fsequences for this xstar

    } // end of section generating fsequences for a previously unseen xstar

    if((importanceSampleSize == 0) && !xstarBool.second) {
        // xstar is already there so use the sequences currently there
            // xstarBool.first is an iterator pointing to a pair (xstar, fsequenceMap)

        // there have been no sequences actually generated
        // need to find out how many fsequences there are:
        StateSequenceSetItSummaryMapIt mapIt;

        for (mapIt = (xstarBool.first->second).begin();
                    mapIt != (xstarBool.first->second).end(); mapIt++) {

            totalFSImpSample += (mapIt->second.countImpSample);
        }

        //debug
        //cout << xstar_label << " has already been done" << std::endl;
    } // end if importanceSampleSize == 0

    // xstarBool.first is an iterator pointing to a pair (xstar, fsequenceMap)

    boost::tuple <XStarSequencesMapIt, size_t, size_t>
            retTuple(xstarBool.first, totalFSImpSample, totalFSGenerated);

    return retTuple;
}





// templatised function to make dot with colour graphs
// type supplied for graph must have colour maps for vertex and edge
template <typename Graph>
void makeColouredUnnamedDotGraph(Graph& graph, string filename= "myGraph.dot")
{
    typedef typename property_map<Graph, vertex_colour_t>::type vertexcolour_map_t;
    typedef typename property_map<Graph, edge_colour_t>::type edgecolour_map_t;

    vertexcolour_map_t vertexColourProp = get(vertex_colour, graph);
    edgecolour_map_t edgeColourProp = get(edge_colour, graph);
    ofstream os(filename.c_str());
    if (os.is_open()) {
        write_graphviz(os, graph, make_unnamed_coloured_label_writer(vertexColourProp),
                        make_unnamed_coloured_label_writer(edgeColourProp));
        os.close();
        std::cout << ".dot file in "
            << filename << std::endl << std::endl;

        }
    else {
        std::cout << "Error: could not open file named "
            << filename << std::endl << std::endl;
    }
}

// templatised function to make black and white (uncoloured) dot graphs
template <typename Graph>
void makeBandWUnnamedDotGraph(Graph& graph, string filename= "myGraph.dot")
{
    ofstream os(filename.c_str());
    if (os.is_open()) {
        write_graphviz(os, graph, make_unnamed_black_and_white_label_writer(),
                        make_unnamed_black_and_white_label_writer());
        os.close();
        std::cout << ".dot file in "
            << filename << std::endl << std::endl;

        }
    else {
        std::cout << "Error: could not open file named "
            << filename << std::endl << std::endl;
    }
}




template <typename T>
void stdoutArray(const T& t)
{
    std::cout << "(";
    for (size_t i = 0; i < t.size(); i++) {
        std::cout << t[i];
        if (i < t.size()-1) std::cout << "\t";
    }
    std::cout << ")";
}



//dijkstra shortest path
template <typename Graph, typename WeightProp, typename ProbPropMap,
            typename Container, typename vertex_descriptor>
double  dij(Graph& graph,
                WeightProp& weightMap,
                ProbPropMap& edgeProbProp,
                Container& shortestPath,
                vertex_descriptor sourceVertex,
                vertex_descriptor sinkVertex,
                size_t nsam,
                const string addColour="black")
{
    std::cout << std::scientific << std::showpoint << std::setprecision(6);

    // typedefs
    typedef typename boost::graph_traits<Graph>::edge_descriptor edge_descriptor;
    typedef typename boost::graph_traits<Graph>::vertex_iterator vertex_iterator;
    typedef typename boost::graph_traits<Graph>::edge_iterator edge_iterator;

    // container for parents
    vector<vertex_descriptor> parent(num_vertices(graph));

    // container for distances
    vector<double> shortestDistances(num_vertices(graph), 0.0);  // all initialised to 0.0

    //start with each vertex as its own parent
    for (pair<vertex_iterator, vertex_iterator> p = vertices(graph);
                            p.first != p.second; ++p.first) {
        parent[*(p.first)] = *(p.first);
    }

    // call algorithm
    dijkstra_shortest_paths(graph, sourceVertex,
                    weight_map(weightMap)
                    .predecessor_map(&parent[0])
                    .distance_map(&shortestDistances[0]));

    /*
    fxz_named_name_map_t nameProp = get(vertex_name, graph);


    cout << "highest probabilities and parents" << endl;
    for (pair<fxz_named_graph_vertex_it, fxz_named_graph_vertex_it> p = vertices(graph);
                            p.first != p.second; ++p.first) {

        double d = shortestDistances[*(p.first)];
        double pr = 1.0/(exp(d));

        cout << "distance to " << nameProp[*(p.first)]
                << " = " << pr
                << ", ";
        cout << "parent of " << nameProp[*(p.first)]
                << " = " << nameProp[parent[*(p.first)]]
                << ", " << endl;

    }
    */

    // container for shortest path
    shortestPath.clear();
    shortestPath.resize(nsam);

    // colour property maps
    fxz_named_vertexcolour_map_t vertexColourProp = get(vertex_colour, graph);
    fxz_named_edgecolour_map_t edgeColourProp = get(edge_colour, graph);

    vertex_descriptor predecessor = sinkVertex;
    size_t index = nsam-1;
    shortestPath[index] = predecessor;

    vertexColourProp[predecessor] = addColour;

    double proby = 1.0;

    typedef typename property_map<Graph, vertex_colour_t>::type vertexcolour_map_t;
    typedef typename property_map<Graph, edge_colour_t>::type edgecolour_map_t;

    while((predecessor != sourceVertex) && (index > 0)) {
        index--;
        vertex_descriptor antecessor = predecessor;
        predecessor = parent[antecessor];
        std::pair<edge_descriptor, bool> pp = edge(predecessor, antecessor, graph);
        shortestPath[index] = predecessor;
        proby *= get(edgeProbProp, pp.first); // accumulate probabilities
    }

    return proby;
}


// breadth-first searching
class bfs_average_edges_visitor:public default_bfs_visitor {

public:
    bfs_average_edges_visitor(EpochVerticesInfoMap* emp) : epochMapPtr(emp) {};

    template < typename Vertex, typename Graph >
    void discover_vertex(Vertex u, const Graph & g)
    {
        BFSInfo info(1, 0, 0);
        size_t init = 0;
        //find the state of this vertex
        //typedef typename property_map<Graph, vertex_fstateSetIt_t>::type fstateIt_map_t;

        state_type fs = *(get(get(vertex_fstateSetIt, g), u));
        size_t epoch = accumulate(fs.begin(), fs.end(), init);

        EpochVerticesInfoMapBool epochBool =
            (*epochMapPtr).insert(make_pair(epoch, info));

        if(!epochBool.second) { // already there
            epochBool.first->second.countVert++;

        }
        // add the outedges
        epochBool.first->second.countAllOutEdge += (boost::out_degree(u, g));
        // number of in-edges must be the epoch above's out edges
        if (epoch > 1) {
            EpochVerticesInfoMapIt it = (*epochMapPtr).find(epoch-1);
            epochBool.first->second.countAllInEdge = it->second.countAllOutEdge;
        }

    }

    EpochVerticesInfoMap* epochMapPtr; // map ptr is public
};

//bread-first search using a bfs_averageouts_visitor
template <typename Graph, typename vertex_descriptor>
void  mapAverageEdges(Graph& graph, vertex_descriptor sourceVertex)
{
    EpochVerticesInfoMap epochMap;
    bfs_average_edges_visitor vis(&epochMap);

    breadth_first_search(graph, sourceVertex, visitor(vis));

    // vis's epochMap should now have data from the vertices



    EpochVerticesInfoMapIt it;

    for (it = epochMap.begin(); it != epochMap.end(); it++) {
        cout << "\nepoch " << (it->first) << " has " << (it->second).countVert << " vertices" << endl;
        cout << "and the total number of out-edges for this epoch is "
                <<  (it->second).countAllOutEdge << endl;
        cout << "so the average number of out-edges per vertex for this epoch is "
                <<  ((it->second).countAllOutEdge/(it->second).countVert) << endl;
        cout << "and the total number of in-edges for this epoch is "
                <<  (it->second).countAllInEdge << endl;
        cout << "so the average number of in-edges per vertex for this epoch is "
                <<  ((it->second).countAllInEdge/(it->second).countVert) << endl;
    }
}


// make possible xstars from a given xstar by making all further index postions 0 or 1
// called recursively to fill up the container xstars
state_typeSet& makePossibleXStars(state_typeSet& xstars, state_type xstar, const size_t index);

// get the Cherry coalescences statistic from an fsequence
size_t getCherryCoalescentStat(const StateSequence& fs);

// get the uneven coalescences statistic from an fsequence
size_t getUnevenCoalescentStat(const StateSequence& fs);

#endif
