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

/*! \file ebc_sfsgraphtypes.hpp
    \brief typedefs for BGL graphs and associated structures
*/

#ifndef __SFSGRAPHTYPES_HPP__
#define __SFSGRAPHTYPES_HPP__


#include <boost/graph/graph_utility.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/tuple/tuple.hpp> // to use tuples


// maximum weight a graph edge can have
// used for shortest path routines
#define MAX_WEIGHT 100


using namespace std;
using namespace boost;

/** @name Fsequence types.
Types related to fsequences. */

//@{

/*! \brief Types for f and z states.

An fstate is a representation of the state of a tree between two coalescent events,
the value at the ith position showing the number of lineages subtending i leaves.

A zstate is associated with a sequence of fstates and an sfs and is a representation
of the xstar which would apply to the sfs in a particular epoch given the sequence of
fstates.
*/
typedef vector<size_t> state_type;
typedef state_type::iterator state_type_itr;



/*! Class for sorting state_type lexicographically.
*/
class StateTypeSorting
{
  public:
    bool operator() (const state_type& s1, const state_type& s2) const {
        return lexicographical_compare(s1.begin(), s1.end(), s2.begin(), s2.end());
    }
};

/*! Types for sets of state_types with lexicographic sorting
*/
typedef set<state_type, StateTypeSorting > state_typeSet;
typedef state_typeSet::iterator state_typeSetIt;
typedef pair<state_typeSetIt, bool> state_typeSetBool;



/*! A class for sorting state_typesSetIts by lexicographic ordering of the state_type pointed to.
*/
class SetItSorting
{
  public:
    bool operator() (const state_typeSetIt& s1, const state_typeSetIt& s2) const
    {
        StateTypeSorting lSort;
        return (lSort((*s1), (*s2)));
    }
};


/*! Type for a pair of iterators into sets of state_types
*/
typedef pair<state_typeSetIt, state_typeSetIt> state_typeSetItPair;


/*!  A class for sorting state_typesSetItPairs by lexicographic ordering of the state_types pointed to.
*/
class SetItPairSorting
{
  public:
    bool operator() (const state_typeSetItPair& p1, const state_typeSetItPair& p2) const
    {
        bool retValue = false;
        StateTypeSorting lSort;
        if (lSort((*(p1.first)), (*(p2.first)))) {// p1.1 < p2.1
            retValue = true;
        }
        else if (!lSort((*(p2.first)), (*(p1.first)))) {   //  p1.1 >= p2.1 && p2.1 >= p1.1
            if (lSort((*(p1.second)), (*(p2.second)))) // p1.2 < p2.2
                retValue = true;
        }

        return retValue;
    }
};


/*!  Types to represent an f-sequence as a sequence of iterators into a set of f-states.
*/
typedef std::vector < state_typeSetIt > StateSequence;
typedef StateSequence::iterator StateSequenceIt;
typedef StateSequence::const_iterator StateSequenceConstIt;

/*!  A class for sorting state sequences by lexicographic ordering of the state_types pointed to.
*/
class StateSequenceSorting
{
    public:
        bool operator() (const StateSequence& s1, const StateSequence& s2) const
        {
            SetItSorting sitSort;
            StateSequenceConstIt it1;
            StateSequenceConstIt it2;
            bool retValue = false;
            for (it1 = s1.begin(), it2 = s2.begin();
                    it1 < s1.end(), it2 < s2.end();
                    it1++, it2++) {
                // make a local comparison of pairs of states
                if (sitSort(*it1, *it2)) {
                    retValue = true;
                    break;          // break out of the for loop

                }
                if (sitSort(*it2, *it1)) {
                    break;          // break out of the loop with retValue false
                }
                // if we are still here, the local comparison must be equal
            }
            // will return false if equal, which is what we want for strict weak ordering
            return retValue;
        }
};


/*! Types for a set of StateSequences with lexicographic ordering.
*/
typedef std::set< StateSequence, StateSequenceSorting > StateSequenceSet;
typedef StateSequenceSet::iterator StateSequenceSetIt;
typedef StateSequenceSet::const_iterator StateSequenceSetConstIt;
typedef pair<StateSequenceSetIt, bool> StateSequenceSetItBool;


/*!  A class for sorting StateSequenceSetIts by lexicographic ordering of the StateSequences pointed to.
*/
class StateSequenceItSorting
{
    public:
        bool operator() (const StateSequenceSetIt& ss_it1,
                            const StateSequenceSetIt& ss_it2) const
        {
            StateSequenceSorting ssSort;
            return (ssSort(*ss_it1, *ss_it2));
        }
};
/*! Types for a container of state_typeSetIts.
*/
typedef vector<state_typeSetIt> VecSetIt;
typedef VecSetIt::iterator VecSetIt_it;

/*! Types for a container of probabilities as doubles.
*/
typedef vector<double> VecProbs;
typedef VecProbs::iterator VecProbsIt;

//@}


/** @name Graph types

Types related to graphs. */
//@{

/** @ name Custom property tags. */

//@{
/*! Tags for iterators to f_states and z_states.
*/
enum vertex_fstateSetIt_t { vertex_fstateSetIt = 151 };
enum vertex_zstateSetIt_t { vertex_zstateSetIt = 161 };
/*! Tag for edge pure kingman coalescent probability.
*/
enum edge_kcprob_t { edge_kcprob = 111 };
/*! Tag for edge reweighted kingman coalescent probability.
*/
enum edge_kcrwprob_t { edge_kcrwprob = 121 };
/*! Tag for whether we have followed an edge.
*/
enum edge_followed_t { edge_followed = 131 };
/*! Tags for vertex colour and edge colour.
*/
enum vertex_colour_t { vertex_colour = 211 };
enum edge_colour_t { edge_colour = 221 };
/*! Tag for accumulating an f-sequence count.
*/
enum vertex_accumulation_t { vertex_accumulation = 611 };
/*! Put the custom properties into the namespace.
*/
namespace boost {
    BOOST_INSTALL_PROPERTY(vertex, fstateSetIt);
    BOOST_INSTALL_PROPERTY(vertex, zstateSetIt);
    BOOST_INSTALL_PROPERTY(edge, kcrwprob);
    BOOST_INSTALL_PROPERTY(edge, kcprob);
    BOOST_INSTALL_PROPERTY(edge, followed);
    BOOST_INSTALL_PROPERTY(vertex, colour);
    BOOST_INSTALL_PROPERTY(edge, colour);
}

//@}


/*! Define vertex properties:  name, fstateSetIt, zstateSetIt, colour.
*/
typedef property<vertex_name_t, std::string,
    property<vertex_fstateSetIt_t, state_typeSetIt,
    property<vertex_zstateSetIt_t, state_typeSetIt,
    property<vertex_colour_t, string > > > > f_x_z_named_graph_vertex_properties;

/*! Define vertex properties:  fstateSetIt, zstateSetIt.
*/
typedef property<vertex_fstateSetIt_t, state_typeSetIt,
    property<vertex_zstateSetIt_t, state_typeSetIt > >
    f_x_z_efficient_graph_vertex_properties;

/*! Define vertex properties:  fstateSetIt.
*/
typedef property<vertex_fstateSetIt_t, state_typeSetIt >
    f_efficient_graph_vertex_properties;


/*! Define edge properties for all graphs.
*/
typedef property<edge_kcrwprob_t, double,
    property<edge_kcprob_t, double,
    property<edge_weight_t, double,
    property<edge_followed_t, bool,
    property<edge_colour_t, string > > > > >fs_graph_edge_properties;

/*! Type for the BGL graph without properties.
*/
typedef boost::adjacency_list<
    setS,       // store out-edges in a stl::set
    vecS,       // store vertices in an stl::vector
    bidirectionalS   // bidirectional graph
    > fseq_graph;

/** @name Graph types.
*/
//@{

/*! \brief Type for the BGL graph on product space with named vertices.

This graph is on the product space of f-state and z-state.  ie, every ordered pair (f, z) has a separate vertex.

Each vertex has a corresponding name 'label'.

Both the Kingman Coalescent transition probability and the (z-state dependent) proposal probabilities are stored as edge properties.

This is an inefficient (more memory hungry) graph but is a convenient type for small graphs to be drawn as dot diagrams with vertex labels indicating the f-state and z-state.

This graph type can also have coloured vertices and edges, again designed for pretty representation in dot graphs.

Like the f_x_z_efficient_graph type, for small n (< about 13)this graph type can also be used to produce fsequences whose proposal probability sums to 1 more quickly than the f_efficient_graph, ie fewer fsequences need be produced in total to cover the whole proposal probability.  For very small n this difference makes little difference to the noticeable speed of the process but for n between around 8 and 13 it could make a more significant difference.
*/
typedef boost::adjacency_list<
    setS,       // store out-edges in a stl::set
    vecS,       // store vertices in an stl::vector
    bidirectionalS,   // bidirectional graph
    f_x_z_named_graph_vertex_properties,
    fs_graph_edge_properties
    > f_x_z_named_graph;

/*! \brief Type for the BGL graph on product space  with unnamed verices.

This graph is on the product space of f-state and z-state.  ie, every ordered pair (f, z) has a separate vertex.

There are no vertex name 'labels'.

Both the Kingman Coalescent transition probability and the (z-state dependent) proposal probabilities are stored as edge properties.

The efficiency of this graph depends on the criteria.  If the graph is used to produce a whole graph (all the states, given a starting (f, z) pair) then this type will take more memory than the f_efficient_graph (there are more vertices and memory is also needed to hold the vertex -> z-state property map), but the graph will not be produced more efficiently than the f_efficient_graph.  However, these memory considerations are only important for very large n (around 70 and above).

If the graph is being used to produce fsequences, the larger number of vertices will consume more memory, but the actual process of finding an fsequence is less computationally intensive than the f_efficient_graph because the z-states for each vertex are explicit and less computation is needed to find 'permissable' edges from a current (f,z) pair and their proposal probabilties.

For small n (< about 13), this graph type can also be used to produce fsequences whose proposal probability sums to 1 more quickly than the f_efficient_graph, ie fewer fsequences need be produced in total to cover the whole proposal probability.  For very small n this difference makes little difference to the noticeable speed of the process but for n between around 8 and 13 it could make a more significant difference.
*/
typedef boost::adjacency_list<
    setS,       // store out-edges in a stl::set
    vecS,       // store vertices in an stl::vector
    bidirectionalS,   // bidirectional graph
    f_x_z_efficient_graph_vertex_properties,
    fs_graph_edge_properties
    > f_x_z_efficient_graph;

/*! \brief Type for the BGL graph using iterators into a set for vertex f state for efficient implementation of fstate space graph.

This graph is only on the f-state space. There are no vertex name 'labels'.

Only the Kingman Coalescent transition probability is stored as an edge property.

The efficiency of this graph depends on the criteria.  If the graph is used to produce a whole graph (all the states, given a starting f state and z state then this type will less memory than the f_x_z_efficient_graph (because there are fewer vertices and no need for the vertex -> z-state property map) but will not be produced more efficiently than the f_efficient_graph.  For very large n - greater than about 70 - this is probably the best graph to use.

If the graph is being used to produce fsequences, the larger number of vertices will consume less memory, but the actual process of finding an fsequence is more computationally intensive than the f_x_z_efficient_graph because the permissable edges from each vertex, and their proposal probabilities, must be recalculated at each move.

This graph type will take longer than the f_x_z_efficient_graph to find all the unique fewer fsequences which cover the whole proposal probability.  This only matters for smaller values of n (being most noticeable for n between around 8 and 13).  For larger values of n it is too computionally lenghtly with either graph type to try to produce all the fsequences (proposal probabilities summing to 1) and this graph type's smaller size, in memory, makes it a better type for the Importance Sampling methods which will then be needed.
*/
typedef boost::adjacency_list<
    setS,       // store out-edges in a stl::set
    vecS,       // store vertices in an stl::vector
    bidirectionalS,   // bidirectional graph
    f_efficient_graph_vertex_properties,
    fs_graph_edge_properties
    > f_efficient_graph;


//@}

/** @name Graph trait and property types.
*/
//@{

/** @name Types for fseq_graph. */

//@{
/*! Vertex and edge descriptors. */
typedef boost::graph_traits<fseq_graph>::vertex_descriptor fseq_graph_vertex_descriptor;
typedef boost::graph_traits<fseq_graph>::edge_descriptor fseq_graph_edge_descriptor;

/*! Size types. */
typedef graph_traits<fseq_graph>::degree_size_type degree_size_type_t;
typedef graph_traits<fseq_graph>::vertices_size_type vertices_size_type_t;
typedef graph_traits<fseq_graph>::edges_size_type edges_size_type_t;

//@}

/** @name Types associated with for f_x_y_named_graph. */

//@{
/*! Iterators. */
typedef graph_traits<f_x_z_named_graph>::vertex_iterator fxz_named_graph_vertex_it;
typedef graph_traits<f_x_z_named_graph>::edge_iterator fxz_named_graph_edge_it;
typedef graph_traits<f_x_z_named_graph>::out_edge_iterator fxz_named_graph_out_edge_it;
typedef graph_traits<f_x_z_named_graph>::in_edge_iterator fxz_named_graph_in_edge_it;
typedef graph_traits<f_x_z_named_graph>::adjacency_iterator fxz_named_graph_adjacency_it;

/*! Property maps. */
typedef property_map<f_x_z_named_graph, vertex_name_t>::type fxz_named_name_map_t;
typedef property_map<f_x_z_named_graph, vertex_fstateSetIt_t>::type fxz_named_fstateIt_map_t;
typedef property_map<f_x_z_named_graph, vertex_zstateSetIt_t>::type fxz_named_zstateIt_map_t;
typedef property_map<f_x_z_named_graph, vertex_colour_t>::type fxz_named_vertexcolour_map_t;
typedef property_map<f_x_z_named_graph, edge_weight_t>::type fxz_named_weight_map_t;
typedef property_map<f_x_z_named_graph, edge_kcrwprob_t>::type fxz_named_proposal_prob_map_t;
typedef property_map<f_x_z_named_graph, edge_kcprob_t>::type fxz_named_kc_prob_map_t;
typedef property_map<f_x_z_named_graph, edge_followed_t>::type fxz_named_followed_map_t;
typedef property_map<f_x_z_named_graph, edge_colour_t>::type fxz_named_edgecolour_map_t;

/*! Vertex and edge descriptors. */
typedef boost::graph_traits<f_x_z_named_graph>::vertex_descriptor fxz_named_graph_vertex_descriptor;
typedef boost::graph_traits<f_x_z_named_graph>::edge_descriptor fxz_named_graph_edge_descriptor;


/*! Size types.*/
typedef graph_traits<f_x_z_named_graph>::degree_size_type fxz_named_degree_size_type;
typedef graph_traits<f_x_z_named_graph>::vertices_size_type fxz_named_vertices_size_type;
typedef graph_traits<f_x_z_named_graph>::edges_size_type fxz_named_edges_size_type;

//@}

/** @name Types for f_x_y_efficient_graph. */

//@{
/*! Iterators. */
typedef graph_traits<f_x_z_efficient_graph>::vertex_iterator fxz_efficient_graph_vertex_it;
typedef graph_traits<f_x_z_efficient_graph>::edge_iterator fxz_efficient_graph_edge_it;
typedef graph_traits<f_x_z_efficient_graph>::out_edge_iterator fxz_efficient_graph_out_edge_it;
typedef graph_traits<f_x_z_efficient_graph>::in_edge_iterator fxz_efficient_graph_in_edge_it;
typedef graph_traits<f_x_z_efficient_graph>::adjacency_iterator fxz_efficient_graph_adjacency_it;

/*! Property maps. */
typedef property_map<f_x_z_efficient_graph, vertex_name_t>::type fxz_efficient_name_map_t;
typedef property_map<f_x_z_efficient_graph, vertex_fstateSetIt_t>::type fxz_efficient_fstateIt_map_t;
typedef property_map<f_x_z_efficient_graph, vertex_zstateSetIt_t>::type fxz_efficient_zstateIt_map_t;
typedef property_map<f_x_z_efficient_graph, edge_kcrwprob_t>::type fxz_efficient_proposal_prob_map_t;
typedef property_map<f_x_z_efficient_graph, edge_kcprob_t>::type fxz_efficient_kc_prob_map_t;
typedef property_map<f_x_z_efficient_graph, edge_weight_t>::type fxz_efficient_weight_map_t;
typedef property_map<f_x_z_efficient_graph, edge_followed_t>::type fxz_efficient_followed_map_t;


/*! Vertex and edge descriptors. */
typedef boost::graph_traits<f_x_z_efficient_graph>::vertex_descriptor fxz_efficient_graph_vertex_descriptor;
typedef boost::graph_traits<f_x_z_efficient_graph>::edge_descriptor fxz_efficient_graph_edge_descriptor;

/*! Size types. */
typedef graph_traits<f_x_z_efficient_graph>::degree_size_type fxz_efficient_degree_size_type;
typedef graph_traits<f_x_z_efficient_graph>::vertices_size_type fxz_efficient_vertices_size_type;
typedef graph_traits<f_x_z_efficient_graph>::edges_size_type fxz_efficient_edges_size_type;

//@}

/** @name Types for f_efficient_graph. */

//@{
/*! Iterators. */
typedef graph_traits<f_efficient_graph>::vertex_iterator f_efficient_graph_vertex_it;
typedef graph_traits<f_efficient_graph>::edge_iterator f_efficient_graph_edge_it;
typedef graph_traits<f_efficient_graph>::out_edge_iterator f_efficient_graph_out_edge_it;
typedef graph_traits<f_efficient_graph>::in_edge_iterator f_efficient_graph_in_edge_it;
typedef graph_traits<f_efficient_graph>::adjacency_iterator f_efficient_graph_adjacency_it;

/*! Property maps. */
typedef property_map<f_efficient_graph, vertex_name_t>::type f_efficient_name_map_t;
typedef property_map<f_efficient_graph, vertex_fstateSetIt_t>::type f_efficient_fstateIt_map_t;
/*! \note  f_efficient does not use the kcrwprob property. */
typedef property_map<f_efficient_graph, edge_kcprob_t>::type f_efficient_kc_prob_map_t;
typedef property_map<f_efficient_graph, edge_weight_t>::type f_efficient_weight_map_t;
/*! \note f_efficient does not use the followed property. */


/*! Vertex and edge descriptors. */
typedef boost::graph_traits<f_efficient_graph>::vertex_descriptor f_efficient_graph_vertex_descriptor;
typedef boost::graph_traits<f_efficient_graph>::edge_descriptor f_efficient_graph_edge_descriptor;

/*! Size types.*/
typedef graph_traits<f_efficient_graph>::degree_size_type f_efficient_degree_size_type;
typedef graph_traits<f_efficient_graph>::vertices_size_type f_efficient_vertices_size_type;
typedef graph_traits<f_efficient_graph>::edges_size_type f_efficient_edges_size_type;

//@}



//@}

/*! \brief Pairs for descriptors and bools. */
typedef std::pair<fxz_named_graph_vertex_descriptor, bool> FxzNamedGraphVertexBoolPair;
typedef std::pair<fxz_named_graph_edge_descriptor, bool> FxzNamedGraphEdgeBoolPair;

typedef std::pair<fxz_efficient_graph_vertex_descriptor, bool> FxzEfficientGraphVertexBoolPair;
typedef std::pair<fxz_efficient_graph_edge_descriptor, bool> FxzEfficientGraphEdgeBoolPair;

typedef std::pair<f_efficient_graph_vertex_descriptor, bool> F_EfficientGraphVertexBoolPair;
typedef std::pair<f_efficient_graph_edge_descriptor, bool> F_EfficientGraphEdgeBoolPair;


/*! \brief Tuple for fxz_named_graph.

(vertex, proposal transition probability, kc transition probability, bool, bool).
*/
typedef boost::tuple<fxz_named_graph_vertex_descriptor, double, double, bool, bool> FxzNamedVertexChoiceTuple;

/*! \brief Tuple for fxz_efficient_graph.

(vertex, proposal transition probability, kc transition probability, bool, bool)
*/
typedef boost::tuple<fxz_efficient_graph_vertex_descriptor, double, double, bool, bool> FxzEfficientVertexChoiceTuple;

/*! \brief Tuple for f_efficient_graph.

(vertex, proposal transition probability, kc transition probability, bool)
*/
typedef boost::tuple<f_efficient_graph_vertex_descriptor, double, double, bool> F_EfficientVertexChoiceTuple;


/*! \brief  Pair of state_type. */
typedef pair<state_type, state_type> fz_pair;

/*! \brief Container of fxz_named_graph edge_descriptors.
*/
typedef vector<fxz_named_graph_edge_descriptor> FxzNamedVecEdges;
typedef FxzNamedVecEdges::iterator FxzNamedVecEdgesIt;

/*! \brief Container of fxz_named_graph vertex_descriptors.
*/
typedef vector<fxz_named_graph_vertex_descriptor> FxzNamedVecVertices;
typedef FxzNamedVecVertices::iterator FxzNamedVecVerticesIt;

/*! \brief Container of fxz_efficient_graph edge descriptors.
*/
typedef vector<fxz_efficient_graph_edge_descriptor> FxzEfficientVecEdges;
typedef FxzEfficientVecEdges::iterator FxzEfficientVecEdgesIt;

/*! \brief Container of fxz_efficient_graph vertex_descriptors.
*/
typedef vector<fxz_efficient_graph_vertex_descriptor> FxzEfficientVecVertices;
typedef FxzEfficientVecVertices::iterator FxzEfficientVecVerticesIt;

/*! \brief Container of f_efficient_graph edge descriptors.
*/
typedef vector<f_efficient_graph_edge_descriptor> F_EfficientVecEdges;
typedef F_EfficientVecEdges::iterator F_EfficientVecEdgesIt;

/*! \brief Container of f_efficient_graph vertex_descriptors.
*/
typedef vector<f_efficient_graph_vertex_descriptor> F_EfficientVecVertices;
typedef F_EfficientVecVertices::iterator F_EfficientVecVerticesIt;

//@}


/** @name Mapping types used with or for graphs. */
//@{

/*! \brief Map from pairs of state_typeSetIts to fseq_graph_vertex_descriptor. */
typedef map< state_typeSetItPair, fseq_graph_vertex_descriptor, SetItPairSorting> ExtSetItPairVertexMap;
typedef ExtSetItPairVertexMap::iterator ExtSetItPairVertexMapIt;
typedef pair<ExtSetItPairVertexMapIt, bool>  ExtSetItPairVertexMapBool;

/*! \brief Map from a state_typeSetIt to fseq_graph_vertex_descriptor. */
typedef map< state_typeSetIt, fseq_graph_vertex_descriptor, SetItSorting> ExtSetItVertexMap;
typedef ExtSetItVertexMap::iterator ExtSetItVertexMapIt;
typedef pair<ExtSetItVertexMapIt, bool>  ExtSetItVertexMapBool;

/*! \brief Struct for tracking counts and probabilities of f-sequences. */
struct FseqInfo
{
    size_t countImpSample;     /*!< count for importance sampling. */
    size_t countGenerated;     /*!< count of newly generated sfs. */
    double proposalProb;
    double kcProb;         /*! probability/weight */

    FseqInfo(size_t ic, size_t gc, double pp, double kcp)
        : countImpSample(ic), countGenerated(gc), proposalProb(pp), kcProb(kcp) {} ; /*! constructor*/
};

/*! \brief  Map from a sequence set iterator to information about that sequence.
*/
typedef std::map< StateSequenceSetIt, FseqInfo, StateSequenceItSorting > StateSequenceSetItSummaryMap;
typedef StateSequenceSetItSummaryMap::iterator StateSequenceSetItSummaryMapIt;
typedef std::pair<StateSequenceSetItSummaryMapIt, bool>  StateSequenceSetItSummaryMapBool;


/*! \brief Map from a state (eg xstar) to a map of f-sequence set iterators and their information.

Remember that the same f-sequence will have different probabilities depending on xstar.
*/
typedef std::map < state_type, StateSequenceSetItSummaryMap, StateTypeSorting >
    XStarSequencesMap;
typedef XStarSequencesMap::iterator XStarSequencesMapIt;
typedef pair<XStarSequencesMapIt, bool>  XStarSequencesMapBool;

/** @name Mapping types used for counting possible fsequences with graphs. */
//@{

/*! Map from an n to a count of fsequences.
*/
typedef map<size_t, size_t> N_FSequenceCountMap;
typedef N_FSequenceCountMap::iterator N_FSequenceCountMapIt;
typedef pair<N_FSequenceCountMapIt, bool> N_FSequenceCountMapBool;

/*! Map from a padded xstar to a map from n to count of fsequences.
*/
typedef map < state_type, N_FSequenceCountMap, StateTypeSorting > PaddedXStarToN_FSequenceCountMap;
typedef PaddedXStarToN_FSequenceCountMap::iterator PaddedXStarToN_FSequenceCountMapIt;
typedef pair<PaddedXStarToN_FSequenceCountMapIt, bool> PaddedXStarToN_FSequenceCountMapBool;

//@}

//@}

/** @name Graphviz types.

Types to help with the production .dot graphs with graphviz.
*/
//@{

/*! \brief  Label writer for coloured graphs with name labels. */
template <class Name, class Colour>
class coloured_label_writer {
    public:
        coloured_label_writer(Name _name, Colour _colour) : name(_name), colour(_colour) {}
        template <class VertexOrEdge>
        void operator()(std::ostream& out, const VertexOrEdge& v) const {
          out << "[label=\"" << name[v] << "\", color= " << colour[v] << "]";
        }
    private:
        Name name;
        Colour colour;
};

/*! \brief  Function to make label writer for coloured graphs with name labels. */
template <class Name, class Colour>
inline coloured_label_writer<Name, Colour>
    make_coloured_label_writer(Name n, Colour c) {
return coloured_label_writer<Name, Colour>(n, c);
}

/*! \brief  Label writer for coloured graphs with no name labels. */
template <class Colour>
class unnamed_coloured_label_writer {
    public:
        unnamed_coloured_label_writer(Colour _colour) : colour(_colour) {}
        template <class VertexOrEdge>
        void operator()(std::ostream& out, const VertexOrEdge& v) const {
          out << "[shape=point,arrowhead=none,label=\"\", color= " << colour[v] << "]";
        }
    private:
        Colour colour;
};

template <class Colour>
/*! \brief  Function to make label writer for coloured graphs with no name labels. */
inline unnamed_coloured_label_writer<Colour>
    make_unnamed_coloured_label_writer(Colour c) {
return unnamed_coloured_label_writer<Colour>(c);
}


/*! \brief  Label writer for black and white  graphs with no name labels. */
class unnamed_black_and_white_label_writer {
    public:
        unnamed_black_and_white_label_writer() {}
        template <class VertexOrEdge>
        void operator()(std::ostream& out, const VertexOrEdge& v) const {
          out << "[shape=point,arrowhead=none,label=\"\"]";
        }
};

/*! \brief  Function to make label writer for black and white graphs with no name labels. */
inline unnamed_black_and_white_label_writer
    make_unnamed_black_and_white_label_writer() {
return unnamed_black_and_white_label_writer();
}


//@}

/*! \brief Struct for tracking counts of vertices and out-edges in a breadth-first search.
*/
struct BFSInfo
{
    vertices_size_type_t countVert;
    degree_size_type_t countAllOutEdge;
    degree_size_type_t countAllInEdge;

    BFSInfo(vertices_size_type_t cv, degree_size_type_t coe, degree_size_type_t cie)
            : countVert(cv), countAllOutEdge(coe), countAllInEdge(cie) {};
};


/*! \brief Map from an epoch to information about vertices and edges corresponding to that epoch. */
typedef std::map< size_t, BFSInfo > EpochVerticesInfoMap;
typedef EpochVerticesInfoMap::iterator EpochVerticesInfoMapIt;
typedef std::pair<EpochVerticesInfoMapIt, bool>  EpochVerticesInfoMapBool;



#endif
