/* -*- mode: C++; indent-tabs-mode: nil; -*-
 *
 * This file is a part of LEMON, a generic C++ optimization library.
 *
 * Copyright (C) 2017
 * Egervary Jeno Kombinatorikus Optimalizalasi Kutatocsoport
 * (Egervary Research Group on Combinatorial Optimization, EGRES).
 *
 * Permission to use, modify and distribute this software is granted
 * provided that this copyright notice appears in all copies. For
 * precise terms see the accompanying LICENSE file.
 *
 * This software is provided "AS IS" with no warranty of any kind,
 * express or implied, and with no claim as to its suitability for any
 * purpose.
 *
 */

#ifndef LEMON_COMPACT_GRAPH_H
#define LEMON_COMPACT_GRAPH_H

///\ingroup graphs
///\file
///\brief CompactDigraph class.

#include <lemon/core.h>
#include <lemon/bits/graph_extender.h>

#include <algorithm>

namespace lemon {

  class CompactDigraphBase {

  public:

    CompactDigraphBase()
      : built(false), node_num(0), arc_num(0),
        node_first_out(NULL),
        arc_target(NULL) {}

    ~CompactDigraphBase() {
      if (built) {
        delete[] node_first_out;
        delete[] arc_target;
      }
    }

    class Node {
      friend class CompactDigraphBase;
    protected:
      int id;
      Node(int _id) : id(_id) {}
    public:
      Node() {}
      Node (Invalid) : id(-1) {}
      bool operator==(const Node& node) const { return id == node.id; }
      bool operator!=(const Node& node) const { return id != node.id; }
      bool operator<(const Node& node) const { return id < node.id; }
    };

    class Arc {
      friend class CompactDigraphBase;
    protected:
      int id;
      int source;
      Arc(int _id, int _source) : id(_id), source(_source) {}
    public:
      Arc() { }
      Arc (Invalid) : id(-1), source(-1) {}
      bool operator==(const Arc& arc) const { return id == arc.id; }
      bool operator!=(const Arc& arc) const { return id != arc.id; }
      bool operator<(const Arc& arc) const { return id < arc.id; }
    };

    Node source(const Arc& e) const { return Node(e.source); }
    Node target(const Arc& e) const { return Node(arc_target[e.id]); }

    void first(Node& n) const { n.id = node_num - 1; }
    static void next(Node& n) { --n.id; }

  private:

    void nextSource(Arc& e) const {
      if (e.id == -1) return;
      int last = node_first_out[e.source] - 1;
      while (e.id == last) {
        --e.source;
        last = node_first_out[e.source] - 1;
      }
    }

  public:

    void first(Arc& e) const {
      e.id = arc_num - 1;
      e.source = node_num - 1;
      nextSource(e);
    }
    void next(Arc& e) const {
      --e.id;
      nextSource(e);
    }

    void firstOut(Arc& e, const Node& n) const {
      e.source = n.id;
      e.id = node_first_out[n.id];
      if (e.id == node_first_out[n.id + 1]) e = INVALID;
    }
    void nextOut(Arc& e) const {
      ++e.id;
      if (e.id == node_first_out[e.source + 1]) e = INVALID;
    }

    void firstIn(Arc& e, const Node& n) const {
      first(e);
      while(e != INVALID && target(e) != n) {
        next(e);
      }
    }
    void nextIn(Arc& e) const {
      Node arcTarget = target(e);
      do {
        next(e);
      } while(e != INVALID && target(e) != arcTarget);
    }

    static int id(const Node& n) { return n.id; }
    static Node nodeFromId(int id) { return Node(id); }
    int maxNodeId() const { return node_num - 1; }

    static int id(const Arc& e) { return e.id; }
    Arc arcFromId(int id) const {
      int *l = std::upper_bound(node_first_out, node_first_out + node_num, id) - 1;
      int src = l - node_first_out;
      return Arc(id, src);
    }
    int maxArcId() const { return arc_num - 1; }

    typedef True NodeNumTag;
    typedef True ArcNumTag;

    int nodeNum() const { return node_num; }
    int arcNum() const { return arc_num; }

  private:

    template <typename Digraph, typename NodeRefMap>
    class ArcLess {
    public:
      typedef typename Digraph::Arc Arc;

      ArcLess(const Digraph &_graph, const NodeRefMap& _nodeRef)
        : digraph(_graph), nodeRef(_nodeRef) {}

      bool operator()(const Arc& left, const Arc& right) const {
        return nodeRef[digraph.target(left)] < nodeRef[digraph.target(right)];
      }
    private:
      const Digraph& digraph;
      const NodeRefMap& nodeRef;
    };

  public:

    typedef True BuildTag;

    void clear() {
      if (built) {
        delete[] node_first_out;
        delete[] arc_target;
      }
      built = false;
      node_num = 0;
      arc_num = 0;
    }

    template <typename Digraph, typename NodeRefMap, typename ArcRefMap>
    void build(const Digraph& digraph, NodeRefMap& nodeRef, ArcRefMap& arcRef) {
      typedef typename Digraph::Node GNode;
      typedef typename Digraph::Arc GArc;

      built = true;

      node_num = countNodes(digraph);
      arc_num = countArcs(digraph);

      node_first_out = new int[node_num + 1];

      arc_target = new int[arc_num];

      int node_index = 0;
      for (typename Digraph::NodeIt n(digraph); n != INVALID; ++n) {
        nodeRef[n] = Node(node_index);
        ++node_index;
      }

      ArcLess<Digraph, NodeRefMap> arcLess(digraph, nodeRef);

      int arc_index = 0;
      for (typename Digraph::NodeIt n(digraph); n != INVALID; ++n) {
        int source = nodeRef[n].id;
        std::vector<GArc> arcs;
        for (typename Digraph::OutArcIt e(digraph, n); e != INVALID; ++e) {
          arcs.push_back(e);
        }
        if (!arcs.empty()) {
          node_first_out[source] = arc_index;
          std::sort(arcs.begin(), arcs.end(), arcLess);
          for (typename std::vector<GArc>::iterator it = arcs.begin();
               it != arcs.end(); ++it) {
            int target = nodeRef[digraph.target(*it)].id;
            arcRef[*it] = Arc(arc_index, source);
            arc_target[arc_index] = target;
            ++arc_index;
          }
        } else {
          node_first_out[source] = arc_index;
        }
      }
      node_first_out[node_num] = arc_num;
    }

    template <typename ArcListIterator>
    void build(int n, ArcListIterator first, ArcListIterator last) {
      built = true;

      node_num = n;
      arc_num = static_cast<int>(std::distance(first, last));

      node_first_out = new int[node_num + 1];

      arc_target = new int[arc_num];

      int arc_index = 0;
      for (int i = 0; i != node_num; ++i) {
        node_first_out[i] = arc_index;
        for ( ; first != last && (*first).first == i; ++first) {
          int j = (*first).second;
          LEMON_ASSERT(j >= 0 && j < node_num,
            "Wrong arc list for CompactDigraph::build()");
          arc_target[arc_index] = j;
          ++arc_index;
        }
      }
      LEMON_ASSERT(first == last,
        "Wrong arc list for CompactDigraph::build()");
      node_first_out[node_num] = arc_num;
    }

  protected:
    bool built;
    int node_num;
    int arc_num;
    int *node_first_out;
    int *arc_target;
  };

  typedef DigraphExtender<CompactDigraphBase> ExtendedCompactDigraphBase;


  /// \ingroup graphs
  ///
  /// \brief A static directed graph class.
  ///
  /// \ref CompactDigraph is a highly efficient digraph implementation
  /// similar to \ref StaticDigraph. It is more memory efficient but does
  /// not provide efficient iteration over incoming arcs.
  ///
  /// It stores only one \c int values for each node and one \c int value
  /// for each arc. Its \ref InArcIt implementation is inefficient and
  /// provided only for compatibility with the \ref concepts::Digraph "Digraph concept".
  ///
  /// This type fully conforms to the \ref concepts::Digraph "Digraph concept".
  /// Most of its member functions and nested classes are documented
  /// only in the concept class.
  ///
  /// \sa concepts::Digraph
  class CompactDigraph : public ExtendedCompactDigraphBase {

  private:
    /// Graphs are \e not copy constructible. Use DigraphCopy instead.
    CompactDigraph(const CompactDigraph &) : ExtendedCompactDigraphBase() {};
    /// \brief Assignment of a graph to another one is \e not allowed.
    /// Use DigraphCopy instead.
    void operator=(const CompactDigraph&) {}

  public:

    typedef ExtendedCompactDigraphBase Parent;

  public:

    /// \brief Constructor
    ///
    /// Default constructor.
    CompactDigraph() : Parent() {}

    /// \brief The node with the given index.
    ///
    /// This function returns the node with the given index.
    /// \sa index()
    static Node node(int ix) { return Parent::nodeFromId(ix); }

    /// \brief The arc with the given index.
    ///
    /// This function returns the arc with the given index.
    /// \sa index()
    Arc arc(int ix) { return arcFromId(ix); }

    /// \brief The index of the given node.
    ///
    /// This function returns the index of the the given node.
    /// \sa node()
    static int index(Node node) { return Parent::id(node); }

    /// \brief The index of the given arc.
    ///
    /// This function returns the index of the the given arc.
    /// \sa arc()
    static int index(Arc arc) { return Parent::id(arc); }

    /// \brief Number of nodes.
    ///
    /// This function returns the number of nodes.
    int nodeNum() const { return node_num; }

    /// \brief Number of arcs.
    ///
    /// This function returns the number of arcs.
    int arcNum() const { return arc_num; }

    /// \brief Build the digraph copying another digraph.
    ///
    /// This function builds the digraph copying another digraph of any
    /// kind. It can be called more than once, but in such case, the whole
    /// structure and all maps will be cleared and rebuilt.
    ///
    /// This method also makes possible to copy a digraph to a CompactDigraph
    /// structure using \ref DigraphCopy.
    ///
    /// \param digraph An existing digraph to be copied.
    /// \param nodeRef The node references will be copied into this map.
    /// Its key type must be \c Digraph::Node and its value type must be
    /// \c CompactDigraph::Node.
    /// It must conform to the \ref concepts::ReadWriteMap "ReadWriteMap"
    /// concept.
    /// \param arcRef The arc references will be copied into this map.
    /// Its key type must be \c Digraph::Arc and its value type must be
    /// \c CompactDigraph::Arc.
    /// It must conform to the \ref concepts::WriteMap "WriteMap" concept.
    ///
    /// \note If you do not need the arc references, then you could use
    /// \ref NullMap for the last parameter. However the node references
    /// are required by the function itself, thus they must be readable
    /// from the map.
    template <typename Digraph, typename NodeRefMap, typename ArcRefMap>
    void build(const Digraph& digraph, NodeRefMap& nodeRef, ArcRefMap& arcRef) {
      if (built) Parent::clear();
      Parent::build(digraph, nodeRef, arcRef);
    }

    /// \brief Build the digraph from an arc list.
    ///
    /// This function builds the digraph from the given arc list.
    /// It can be called more than once, but in such case, the whole
    /// structure and all maps will be cleared and rebuilt.
    ///
    /// The list of the arcs must be given in the range <tt>[begin, end)</tt>
    /// specified by STL compatible itartors whose \c value_type must be
    /// <tt>std::pair<int,int></tt>.
    /// Each arc must be specified by a pair of integer indices
    /// from the range <tt>[0..n-1]</tt>. <i>The pairs must be in a
    /// non-decreasing order with respect to their first values.</i>
    /// If the k-th pair in the list is <tt>(i,j)</tt>, then
    /// <tt>arc(k-1)</tt> will connect <tt>node(i)</tt> to <tt>node(j)</tt>.
    ///
    /// \param n The number of nodes.
    /// \param begin An iterator pointing to the beginning of the arc list.
    /// \param end An iterator pointing to the end of the arc list.
    ///
    /// For example, a simple digraph can be constructed like this.
    /// \code
    ///   std::vector<std::pair<int,int> > arcs;
    ///   arcs.push_back(std::make_pair(0,1));
    ///   arcs.push_back(std::make_pair(0,2));
    ///   arcs.push_back(std::make_pair(1,3));
    ///   arcs.push_back(std::make_pair(1,2));
    ///   arcs.push_back(std::make_pair(3,0));
    ///   CompactDigraph gr;
    ///   gr.build(4, arcs.begin(), arcs.end());
    /// \endcode
    template <typename ArcListIterator>
    void build(int n, ArcListIterator begin, ArcListIterator end) {
      if (built) Parent::clear();
      CompactDigraphBase::build(n, begin, end);
      notifier(Node()).build();
      notifier(Arc()).build();
    }

    /// \brief Clear the digraph.
    ///
    /// This function erases all nodes and arcs from the digraph.
    void clear() {
      Parent::clear();
    }

  public:

    Node baseNode(const OutArcIt &arc) const {
      return Parent::source(static_cast<const Arc&>(arc));
    }

    Node runningNode(const OutArcIt &arc) const {
      return Parent::target(static_cast<const Arc&>(arc));
    }

    Node baseNode(const InArcIt &arc) const {
      return Parent::target(static_cast<const Arc&>(arc));
    }

    Node runningNode(const InArcIt &arc) const {
      return Parent::source(static_cast<const Arc&>(arc));
    }

  };

}

#endif
