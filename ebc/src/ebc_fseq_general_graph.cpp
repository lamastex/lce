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
/*! \file ebc_fseq_general_graph.cpp
    \brief general methods associated with use of BGL graphs for unvintaged (sized) n-coalescent and associated shape statistics.

An f-sequence in these files is a 'forward' f-sequence (the final state is the
'current' state representing the tree in the epoch where it has n leaves).
The fsequence does not include the ultimate common ancestor state with one lineage.
*/


#include <sstream>  // to be able to manipulate strings as streams

//#include <ebc_output.hpp>
//#include <ebc_params.hpp>

#include <ebc_fseq_general_graph.hpp>

//#include <gsl/gsl_randist.h>
//#include <gsl/gsl_sf_gamma.h>
//#include <gsl/gsl_matrix.h>
//#include <gsl/gsl_rng.h> // is in the hpp file

//#include <boost/graph/graphviz.hpp>
//#include <boost/graph/dijkstra_shortest_paths.hpp>

#include <numeric> // to use accumulate





// tries to find this xstar or its equivalent in the map
// if it finds it, xstarBool.second is set to false
//      and xstarBool.first->second is a map of fsequences
// if it does not find it, xstarBool.second stays true
//      and xstarBool.first->second is an empty map
XStarSequencesMapBool& findExistingMap(XStarSequencesMapBool& xstarBool,
                                        XStarSequencesMap& xstarToSequencesInfoMap,
                                        const state_type& xstar)
{
    // map to be copied
    StateSequenceSetItSummaryMap fsequenceMap;

    // check if we have got fsequences for this xstar already
    xstarBool =
        xstarToSequencesInfoMap.insert(make_pair(xstar, fsequenceMap));

    if(xstarBool.second) { // no sequences so far for this xstar

        // but we could have sequences for a simpler version
        XStarSequencesMapIt alt_it = xstarToSequencesInfoMap.end();
        if (xstar.size() > 1) {
            state_type xstarAlt1 = xstar;
            xstarAlt1[0] == 0 ? xstarAlt1[0] = 1: xstarAlt1[0] = 0;
            alt_it = xstarToSequencesInfoMap.find(xstarAlt1);
            if (alt_it == xstarToSequencesInfoMap.end()) { // did not find alt1
                state_type xstarAlt2 = xstar;
                xstarAlt2[1] == 0 ? xstarAlt2[1] = 1: xstarAlt2[1] = 0;
                alt_it = xstarToSequencesInfoMap.find(xstarAlt2);
                if (alt_it == xstarToSequencesInfoMap.end()) { // did not find alt2
                    state_type xstarAlt3 = xstarAlt1;
                    xstarAlt3[1] == 0 ? xstarAlt3[1] = 1: xstarAlt3[1] = 0;
                    alt_it = xstarToSequencesInfoMap.find(xstarAlt3);
                }
            }
        }

        if (alt_it != xstarToSequencesInfoMap.end()) { // found an equivalent existing xstar

            //debug
            //state_type alt_xstar = alt_it->first;
            //string label = stateToString(alt_xstar);
            //cout << "found an alternative " << label << endl;

            xstarBool.first->second = alt_it->second; // copy that map
            zeroMapCounts(xstarBool.first->second); // zero the counts
            xstarBool.second = false;
        }
    }
    //else {
      //  cout << "found with xstar" << endl;
    //}
    return xstarBool; // return by reference
}


// zero the newly generated counts in a given map
void zeroMapCounts(StateSequenceSetItSummaryMap& sssMap)
{
    StateSequenceSetItSummaryMapIt mapIt;

    for (mapIt = sssMap.begin(); mapIt != sssMap.end(); mapIt++) {

        mapIt->second.countGenerated = 0;
    }

}


// find new states to go to from current state given current z_state
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
            VecProbs& newProbs)
{
    bool retValue = false;

    // a value to use for checking accumulated probabilities
    // accumulated probabilities have to be < accProbTolerance different to 1.0
    double accProbTolerance = 0.00001;

    // z tells us about states we must have:  if z has an entry in the
    // position corresponding to *leaves being subtended, then
    // the f sequences must include a state where there is
    // at least one lineage subtending *leaves

    // find the maximum such restriction
    // check whether the z state places any restrictions on where we can go
    size_t init = 0;
    size_t zSum = std::accumulate(z_state.begin(), z_state.end(), init);

    size_t maxRestrictingL = 0;
    size_t ll = z_state.size()-1;

    while (zSum > 0 && ll > 0 && maxRestrictingL == 0) {
        if (z_state[ll] != 0)
            maxRestrictingL = ll+1;
        ll--;
    }

    // element at index i, i = 0..nsam-1, is number of lineages subtending i+1 leaves
    // we can only split up a lineage subtending > 1 leaves,
    // ie elements indexed i = 1..nsam-1
    for (size_t i = 1; i < f_state.size(); i++) {

        // can only split a lineage if there are any of the lineage to split
        if (f_state[i] > 0) {

            size_t leavesSubtended = i+1; // number of leaves subtended
            size_t lineages = f_state[i]; // number of lineages

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

            // if we ignore the restrictions imposed by the z state,
            // we can split leavesSubtended into the following unordered pairs
            // (1,leaves-1), (2,leaves-2).. (split_jmax, leaves-split_jmax)
            for (size_t split_j = 1; split_j <= split_jmax; split_j++) {

                size_t split_k = leavesSubtended-split_j;

                // calculate the probability of the new state
                // this depends on what epoch we are going into
                // which we can tell from the new f_state

                // start by saying transition probability is
                // lineages/(nsam-newEpoch-1)
                double uncontrolledTransProb =
                        lineages/static_cast<double>(nsam - newEpoch + 1);
                // double if split_j != split_k
                if (split_j != split_k)
                    uncontrolledTransProb *= 2.0;

                //but is the new state allowed?
                bool allowed = ((maxRestrictingL == 0) ||
                        ((split_j <= split_jhat) && (split_jhat <= split_k)));

                if (allowed) {

                    retValue = true;

                    // the new state is given by removing 1 from the element
                    // at f_state[i], ie the number of lineages subtending i+1 leaves
                    // and adding 1 to each of the elements at [split_j-1], [split_k-1]
                    // ie adding 1 to each of the number of lineages subtending
                    // split_j, split_k leaves

                    // find the new f and z states
                    state_type new_f_state = f_state;
                    new_f_state[leavesSubtended-1] -= 1;
                    new_f_state[split_j-1] += 1;
                    new_f_state[split_k-1] += 1;

                    // z is changed by turning any 1 at
                    // index[split_j-1[ [split_k-1] to 0
                    state_type new_z_state = z_state;
                    if (new_z_state[split_j-1] == 1)
                        new_z_state[split_j-1] = 0;
                    if (new_z_state[split_k-1] == 1)
                        new_z_state[split_k-1] = 0;

                    // get iterator to this state in the fstates set
                    state_typeSetBool pair_f = extFStateSet.insert(new_f_state);
                    // get iterator to this state in the zstates set
                    state_typeSetBool pair_z = extZStateSet.insert(new_z_state);

                    new_f_statesSetIt.push_back(pair_f.first);
                    new_z_statesSetIt.push_back(pair_z.first);
                    newProbs.push_back(uncontrolledTransProb);

                    // accumulate probabilities of edges from this vertex v
                    accProbs += uncontrolledTransProb;

                }
                else {

                    rejectedProbs += uncontrolledTransProb;
                }


            } // end loop through possible splits of number of leaves subtended chosen
        }// end if there is a lineage to split
    }// end loop through possible numbers of leaves subtended

    if ((1.0 - (accProbs + rejectedProbs) > accProbTolerance)
            || (accProbs + rejectedProbs - 1.0 > accProbTolerance)) {

                std::cout << "Warning: All unadjusted accumulated probabilities are "
                                                << accProbs << std::endl;

    }

    return retValue;

}



// add an f_state to the fsequence container
void addStateToF_Sequence(StateSequence& fs, const state_typeSetIt& f_stateIt,
        const size_t epoch, const size_t nsam)
{
    // the fs only has epochs 2.. onwards,
    // ie does not show number of lineages subtending nsam leaves

    if (epoch > 1 && epoch <= nsam) {
        fs.push_back(f_stateIt);
    }
}



// add a z_state to the zsequence container
void addStateToZ_Sequence(StateSequence& zs, const state_typeSetIt& z_stateIt,
        const size_t epoch, const size_t nsam)
{
    // the zs only has epochs 2.. onwards

    if (epoch > 1 && epoch <= nsam) {
        zs.push_back(z_stateIt);
    }
}




// convert a pair of state_types to a single string as a name or label or identifer
string statesToString(const state_type& f_state, const state_type& z_state)
{
    ostringstream stm;

    stm << "((";
    for (size_t f = 0; f < f_state.size(); f++) {
        stm << f_state[f];
        if (f < f_state.size()-1) stm << " ";
    }
    stm << "),(";
    for (size_t z = 0; z < z_state.size(); z++) {
        stm << z_state[z];
        if (z < z_state.size()-1) stm << " ";
    }
    stm << "))";

    return stm.str();
}


// convert a state_type to a single string as a name or label or identifer
string stateToString(const state_type& some_state)
{
    ostringstream stm;

    stm << "(";
    for (size_t f = 0; f < some_state.size(); f++) {
        stm << some_state[f];
        if (f < some_state.size()-1) stm << " ";
    }
    stm << ")";

    return stm.str();
}


// in the final state, all the elements should be in [0], ie subtending 1 lineage
bool checkFinalState(const state_type& f_state, const size_t nsam)
{
    size_t init = 0;

    return ((f_state.size() > 0) && (f_state[0] == nsam)
            && (std::accumulate(f_state.begin(), f_state.end(), init) == nsam));
}





/*

void makeColouredUnnamedDotGraph(f_x_z_named_graph& graph, string filename)
{

    fxz_named_vertexcolour_map_t vertexColourProp = get(vertex_colour, graph);
    fxz_named_edgecolour_map_t edgeColourProp = get(edge_colour, graph);
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
*/

void printStateSequence(const StateSequence& sequence)
{

    StateSequenceConstIt it;
    for (it = sequence.begin(); it < sequence.end(); it++) {
        // the actual state is *(*it)
        state_type s = *(*it);
        std::cout << "( ";
        copy (s.begin(), s.end(), ostream_iterator<size_t>(std::cout, " "));
        std::cout << ")\n";
    }
    std::cout << endl;

}

void printStateSequence(const StateSequenceSetIt& seqSetIt)
{
    StateSequence fseq = *seqSetIt;

    printStateSequence(fseq);

}




/*! This returns the shape statistic Raaz refers to as x* that gives the topology information in the observed sfs

x* gives all information available in the observed sfs regarding the equivalence classes of binary coalescent trees upto f-sequences.

When the i-th element of the sfs is > 0, then we know that there had to have been a lineage in the underlying coalescent tree (or equivalently in the corresponding f-sequence) that lead to the
i-ton mutation(s) observed in sfs[i].

If sfs[i] > 0 then the shape statistic x*[i] = 1 else x*[i]
xstar will only have nsam-1 elements, ie we are not interested in nsam mutations on the same site
This is the same size as the sfs data we are given
*/
state_type& sfsToXstar(state_type& xstar, const sfs_array_type& sfs)
{
    xstar.clear();
    if (sfs.size() > 0) {
        size_t init = 0;
        xstar.resize(sfs.size(), init); // one smaller than the sfs, all 0's

        // can't use transform because of valarray for sfs_array_type?
        for (size_t i = 0; i < xstar.size(); i++) {
            if (sfs[i] > 0)
                xstar[i] = 1;

        }
    }
    return xstar;
}

/* return (EpochTimes^Transpose * fseq), i.e. the lengths of lineages leading to singleton, doubleton,...,(n-1)-ton mutations

    Input 1: the EpochTimes, where EpochTimes[nsam-2] = t_n = Time to 1st Coal event, ..., EpochTimes[0] = t_2 = Time to last Coal event.
    Input 2: the f-sequence fs under unvintaged and sized (Kingman's unlabeled) n-coalescent.
        Input 3: nsam is the sample size n in the n-coalescent.
        Output: LineageLenItons = (EpochTimes^Transpose * fseq)
*/
std::valarray<double>& epochTimesProdFseq(std::valarray<double>& lineageLenItons,
                                        const std::valarray<double> & epochTimes,
                                        const StateSequence& fs,
                                        const size_t nsam)
{

    if (epochTimes.size() != nsam-1) {
        std::cout << "Problem in epochTimesProdFseq: epochTimes is the wrong size\n";
        std::cout << "epochTimes is" << std::endl;
        stdoutArray(epochTimes);
        exit(1);
    }

    // we now know that epochTimes has nsam-1 elements

    // fsSetIt is an iterator to a StateSequence in a set
    // a StateSequence is a vector of state_types
    // my vectors go from state at epoch 2 to state at epoch nsam,
    // ie exclude the state at epoch 1 when one lineage subtends all nsam leaves
    lineageLenItons.resize(nsam-1, 0.0);

    for (size_t j = 0; j < nsam - 1; j++) {
        double l_j = 0.0;
        for (size_t i = 0; i < nsam - j - 1; i++) {

            l_j += epochTimes[i] * static_cast<double>((*(fs[i]))[j]);
        }
        lineageLenItons[j] = l_j;
    }
    //debug
    //stdoutArray(epochTimes);
    //cout << endl;
    //printStateSequence(fs);
    //cout << endl;
    //stdoutArray(lineageLenItons);

    return lineageLenItons; // return by reference
}


// make possible xstars from a given xstar by making all further index postions 0 or 1
// called recursively to fill up the container xstars
state_typeSet& makePossibleXStars(state_typeSet& xstars,
                                    state_type xstar,
                                    const size_t index)
{

    if ((index > 1) && (index < xstar.size())) {
        // take the given x star, make two copies, put both into container, and send off again
        state_typeSetBool insertBool = xstars.insert(xstar);
        makePossibleXStars(xstars, xstar, index + 1);

        xstar[index] = 1;
        insertBool = xstars.insert(xstar);
        makePossibleXStars(xstars, xstar, index + 1);
    }
    //else do nothing


    return xstars; // return by reference
}


// get the Cherry coalescences statistic from an fsequence
// cherry stat is the number of times in the tree that there is a coalescence
// of a pair of leaves
size_t getCherryCoalescentStat(const StateSequence& fs)
{
    size_t cherryStat = 0;

    size_t fsLength = fs.size();

    size_t before = 0; // will be okay for all nsam > 2

    if (fsLength == 1) {
        // case where nsam = 2, but we don't have (0 1) state in the fsequence
        before = 1;
    }

        for (size_t i = 0; i < fsLength; i++) {
            // second element has index [1]
            size_t after = (*(fs[i]))[1];
            cherryStat += ((before - after) == 1 ? 1 : 0);
            before = after;
            }
    return cherryStat;
}

// get the uneven coalescences statistic from an fsequence
// cherry stat is the number of times in the tree that there is a coalescence
// of two sets of distinct size
size_t getUnevenCoalescentStat(const StateSequence& fs)
{
    size_t unevenStat = 0;

    size_t fsLength = fs.size();

    // make the starting state, for epoch 1, which is not in the fsequence
    state_type before(fsLength + 1, 0);
    before[fsLength] = 1; // last element in before will be 1

    for (size_t i = 0; i < fsLength; i++) {
        int maxDiff = 0;
        state_type after = *(fs[i]);
        size_t stateLength = fsLength + 1;
        for (size_t j = 0; j < stateLength; j++) {
            int thisDiff = after[j] - before[j];
            if (thisDiff > maxDiff) maxDiff = thisDiff;
        }
        before = after;

        unevenStat += (maxDiff == 1 ? 1 : 0);
    }
    return unevenStat;
}
