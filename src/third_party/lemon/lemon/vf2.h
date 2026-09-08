/* -*- mode: C++; indent-tabs-mode: nil; -*-
 *
 * This file is a part of LEMON, a generic C++ optimization library.
 *
 * Copyright (C) 2015-2017
 * EMAXA Kutato-fejleszto Kft. (EMAXA Research Ltd.)
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

#ifndef LEMON_VF2_H
#define LEMON_VF2_H

///\ingroup graph_properties
///\file
///\brief VF2 algorithm \cite cordella2004sub.

#include <lemon/core.h>
#include <lemon/concepts/graph.h>
#include <lemon/bfs.h>
#include <lemon/bits/vf2_internals.h>

#include <vector>

namespace lemon {
  namespace bits {
    namespace vf2 {

      class AlwaysEq {
      public:
        template<class T1, class T2>
        bool operator()(T1&, T2&) const {
          return true;
        }
      };

      template<class M1, class M2>
      class MapEq{
        const M1 &_m1;
        const M2 &_m2;
      public:
        MapEq(const M1 &m1, const M2 &m2) : _m1(m1), _m2(m2) { }
        bool operator()(typename M1::Key k1, typename M2::Key k2) const {
          return _m1[k1] == _m2[k2];
        }
      };

      template <class G>
      class BfsLeaveOrder : public BfsVisitor<G> {
        int i;
        const G &_g;
        std::vector<typename G::Node> &_order;
      public:
        BfsLeaveOrder(const G &g, std::vector<typename G::Node> &order)
          : i(0), _g(g), _order(order){
        }
        void process(const typename G::Node &node) {
          _order[i++]=node;
        }
      };
    }
  }


  ///%VF2 algorithm class.

  ///\ingroup graph_isomorphism This class provides an efficient
  ///implementation of the %VF2 algorithm \cite cordella2004sub
  ///for variants of the (Sub)graph Isomorphism problem.
  ///
  ///There is also a \ref vf2() "function-type interface" called \ref vf2()
  ///for the %VF2 algorithm, which is probably more convenient in most
  ///use cases.
  ///
  ///\tparam G1 The type of the graph to be embedded.
  ///The default type is \ref ListGraph.
  ///\tparam G2 The type of the graph g1 will be embedded into.
  ///The default type is \ref ListGraph.
  ///\tparam M The type of the NodeMap storing the mapping.
  ///By default, it is G1::NodeMap<G2::Node>
  ///\tparam NEQ A bool-valued binary functor determinining whether a node is
  ///mappable to another. By default, it is an always-true operator.
  ///
  ///\sa vf2()
#ifdef DOXYGEN
  template<class G1, class G2, class M, class NEQ >
#else
  template<class G1 = ListGraph,
           class G2 = ListGraph,
           class M = typename G1::template NodeMap<G2::Node>,
           class NEQ = bits::vf2::AlwaysEq >
#endif
  class Vf2 {
    //The graph to be embedded
    const G1 &_g1;

    //The graph into which g1 is to be embedded
    const G2 &_g2;

    //Functor with bool operator()(G1::Node,G2::Node), which returns 1
    //if and only if the two nodes are equivalent
    NEQ _nEq;

    //Current depth in the search tree
    int _depth;

    //The current mapping. _mapping[v1]=v2 iff v1 has been mapped to v2,
    //where v1 is a node of G1 and v2 is a node of G2
    M &_mapping;

    //_order[i] is the node of g1 for which a pair is searched if depth=i
    std::vector<typename G1::Node> _order;
 
    //_conn[v2] = number of covered neighbours of v2
    typename G2::template NodeMap<int> _conn;

    //_currEdgeIts[i] is the last used edge iterator in the i-th
    //depth to find a pair for node _order[i]
    std::vector<typename G2::IncEdgeIt> _currEdgeIts;

    //lookup tables for cutting the searchtree
    typename G1::template NodeMap<int> _rNew1t, _rInOut1t;

    MappingType _mapping_type;

    bool _deallocMappingAfterUse;

    //cut the search tree
    template<MappingType MT>
    bool cut(const typename G1::Node n1,const typename G2::Node n2) const {
      int rNew2=0,rInOut2=0;
      for(typename G2::IncEdgeIt e2(_g2,n2); e2!=INVALID; ++e2) {
        const typename G2::Node currNode=_g2.oppositeNode(n2,e2);
        if(_conn[currNode]>0)
          ++rInOut2;
        else if(MT!=SUBGRAPH&&_conn[currNode]==0)
          ++rNew2;
      }
      switch(MT) {
      case INDUCED:
        return _rInOut1t[n1]<=rInOut2&&_rNew1t[n1]<=rNew2;
      case SUBGRAPH:
        return _rInOut1t[n1]<=rInOut2;
      case ISOMORPH:
        return _rInOut1t[n1]==rInOut2&&_rNew1t[n1]==rNew2;
      default:
        return false;
      }
    }

    template<MappingType MT>
    bool feas(const typename G1::Node n1,const typename G2::Node n2) {
      if(!_nEq(n1,n2))
        return 0;

      for(typename G1::IncEdgeIt e1(_g1,n1); e1!=INVALID; ++e1) {
        const typename G1::Node& currNode=_g1.oppositeNode(n1,e1);
        if(_mapping[currNode]!=INVALID)
          --_conn[_mapping[currNode]];
      }
      bool isIso=1;
      for(typename G2::IncEdgeIt e2(_g2,n2); e2!=INVALID; ++e2) {
        int& connCurrNode = _conn[_g2.oppositeNode(n2,e2)];
        if(connCurrNode<-1)
          ++connCurrNode;
        else if(MT!=SUBGRAPH&&connCurrNode==-1) {
          isIso=0;
          break;
        }
      }

      for(typename G1::IncEdgeIt e1(_g1,n1); e1!=INVALID; ++e1) {
        const typename G2::Node& currNodePair=_mapping[_g1.oppositeNode(n1,e1)];
        int& connCurrNodePair=_conn[currNodePair];
        if(currNodePair!=INVALID&&connCurrNodePair!=-1) {
          switch(MT) {
          case INDUCED:
          case ISOMORPH:
            isIso=0;
            break;
          case SUBGRAPH:
            if(connCurrNodePair<-1)
              isIso=0;
            break;
          }
          connCurrNodePair=-1;
        }
      }
      return isIso&&cut<MT>(n1,n2);
    }

    void addPair(const typename G1::Node n1,const typename G2::Node n2) {
      _conn[n2]=-1;
      _mapping.set(n1,n2);
      for(typename G2::IncEdgeIt e2(_g2,n2); e2!=INVALID; ++e2) {
        int& currConn = _conn[_g2.oppositeNode(n2,e2)];
        if(currConn!=-1)
          ++currConn;
      }
    }

    void subPair(const typename G1::Node n1,const typename G2::Node n2) {
      _conn[n2]=0;
      _mapping.set(n1,INVALID);
      for(typename G2::IncEdgeIt e2(_g2,n2); e2!=INVALID; ++e2) {
        int& currConn = _conn[_g2.oppositeNode(n2,e2)];
        if(currConn>0)
          --currConn;
        else if(currConn==-1)
          ++_conn[n2];
      }
    }

    void initOrder() {
      //determine the order in which we will find pairs for the nodes of g1
      //BFS order is more efficient in practice than DFS
      bits::vf2::BfsLeaveOrder<G1> v(_g1,_order);
      BfsVisit<G1,bits::vf2::BfsLeaveOrder<G1> > bfs(_g1, v);
      bfs.run();
    }

    template<MappingType MT>
    bool extMatch() {
      while(_depth>=0) {
        if(_depth==static_cast<int>(_order.size())) {
          //all nodes of g1 are mapped to nodes of g2
          --_depth;
          return true;
        }
        typename G1::Node& nodeOfDepth = _order[_depth];
        const typename G2::Node& pairOfNodeOfDepth = _mapping[nodeOfDepth];
        typename G2::IncEdgeIt &edgeItOfDepth = _currEdgeIts[_depth];
        //the node of g2 whose neighbours are the candidates for
        //the pair of nodeOfDepth
        typename G2::Node currPNode;
        if(edgeItOfDepth==INVALID) {
          typename G1::IncEdgeIt fstMatchedE(_g1,nodeOfDepth);
          //if pairOfNodeOfDepth!=INVALID, we don't use fstMatchedE
          if(pairOfNodeOfDepth==INVALID) {
            for(; fstMatchedE!=INVALID &&
                  _mapping[_g1.oppositeNode(nodeOfDepth,
                                            fstMatchedE)]==INVALID;
                ++fstMatchedE) ; //find fstMatchedE
          }
          if(fstMatchedE==INVALID||pairOfNodeOfDepth!=INVALID) {
            //We found no covered neighbours, this means that
            //the graph is not connected (or _depth==0). Each
            //uncovered (and there are some other properties due
            //to the spec. problem types) node of g2 is
            //candidate. We can read the iterator of the last
            //tried node from the match if it is not the first
            //try (match[nodeOfDepth]!=INVALID)
            typename G2::NodeIt n2(_g2);
            //if it's not the first try
            if(pairOfNodeOfDepth!=INVALID) {
              n2=++typename G2::NodeIt(_g2,pairOfNodeOfDepth);
              subPair(nodeOfDepth,pairOfNodeOfDepth);
            }
            for(; n2!=INVALID; ++n2)
              if(MT!=SUBGRAPH) {
                if(_conn[n2]==0&&feas<MT>(nodeOfDepth,n2))
                  break;
              }
              else if(_conn[n2]>=0&&feas<MT>(nodeOfDepth,n2))
                break;
            // n2 is the next candidate
            if(n2!=INVALID){
              addPair(nodeOfDepth,n2);
              ++_depth;
            }
            else // there are no more candidates
              --_depth;
            continue;
          }
          else {
            currPNode=_mapping[_g1.oppositeNode(nodeOfDepth,
                                                fstMatchedE)];
            edgeItOfDepth=typename G2::IncEdgeIt(_g2,currPNode);
          }
        }
        else {
          currPNode=_g2.oppositeNode(pairOfNodeOfDepth,
                                     edgeItOfDepth);
          subPair(nodeOfDepth,pairOfNodeOfDepth);
          ++edgeItOfDepth;
        }
        for(; edgeItOfDepth!=INVALID; ++edgeItOfDepth) {
          const typename G2::Node currNode =
            _g2.oppositeNode(currPNode, edgeItOfDepth);
          if(_conn[currNode]>0&&feas<MT>(nodeOfDepth,currNode)) {
            addPair(nodeOfDepth,currNode);
            break;
          }
        }
        edgeItOfDepth==INVALID?--_depth:++_depth;
      }
      return false;
    }

    //calculate the lookup table for cutting the search tree
    void initRNew1tRInOut1t() {
      typename G1::template NodeMap<int> tmp(_g1,0);
      for(unsigned int i=0; i<_order.size(); ++i) {
        const typename G1::Node& orderI = _order[i];
        tmp[orderI]=-1;
        for(typename G1::IncEdgeIt e1(_g1,orderI); e1!=INVALID; ++e1) {
          const typename G1::Node currNode=_g1.oppositeNode(orderI,e1);
          if(tmp[currNode]>0)
            ++_rInOut1t[orderI];
          else if(tmp[currNode]==0)
            ++_rNew1t[orderI];
        }
        for(typename G1::IncEdgeIt e1(_g1,orderI); e1!=INVALID; ++e1) {
          const typename G1::Node currNode=_g1.oppositeNode(orderI,e1);
          if(tmp[currNode]!=-1)
            ++tmp[currNode];
        }
      }
    }
  public:
    ///Constructor

    ///Constructor

    ///\param g1 The graph to be embedded into \e g2.
    ///\param g2 The graph \e g1 will be embedded into.
    ///\param m \ref concepts::ReadWriteMap "read-write" NodeMap
    ///storing the found mapping.
    ///\param neq A bool-valued binary functor determining whether a node is
    ///mappable to another. By default it is an always true operator.
    Vf2(const G1 &g1, const G2 &g2, M &m, const NEQ &neq = NEQ() ) :
      _g1(g1), _g2(g2), _nEq(neq), _depth(0), _mapping(m),
      _order(countNodes(g1)), _conn(g2,0),
      _currEdgeIts(countNodes(g1),INVALID), _rNew1t(g1,0), _rInOut1t(g1,0),
      _mapping_type(SUBGRAPH), _deallocMappingAfterUse(0)
    {
      initOrder();
      initRNew1tRInOut1t();
      for(typename G1::NodeIt n(g1);n!=INVALID;++n)
        m[n]=INVALID;
    }

    ///Destructor

    ///Destructor.
    ///

    ~Vf2(){
      if(_deallocMappingAfterUse)
        delete &_mapping;
    }

    ///Returns the current mapping type

    ///Returns the current mapping type
    ///
    MappingType mappingType() const {
      return _mapping_type;
    }
    ///Sets mapping type

    ///Sets mapping type.
    ///
    ///The mapping type is set to \ref SUBGRAPH by default.
    ///
    ///\sa See \ref MappingType for the possible values.
    void mappingType(MappingType m_type) {
      _mapping_type = m_type;
    }

    ///Finds a mapping

    ///It finds a mapping from g1 into g2 according to the mapping
    ///type set by \ref mappingType(MappingType) "mappingType()".
    ///
    ///By subsequent calls, it returns all possible mappings one-by-one.
    ///
    ///\retval true if a mapping is found.
    ///\retval false if there is no (more) mapping.
    bool find() {
      switch(_mapping_type) {
      case SUBGRAPH:
        return extMatch<SUBGRAPH>();
      case INDUCED:
        return extMatch<INDUCED>();
      case ISOMORPH:
        return extMatch<ISOMORPH>();
      default:
        return false;
      }
    }
  };

  template<class G1, class G2>
  class Vf2WizardBase {
  protected:
    typedef G1 Graph1;
    typedef G2 Graph2;

    const G1 &_g1;
    const G2 &_g2;

    MappingType _mapping_type;

    typedef typename G1::template NodeMap<typename G2::Node> Mapping;
    bool _local_mapping;
    void *_mapping;
    void createMapping() {
      _mapping = new Mapping(_g1);
    }

    void *myVf2; //used in Vf2Wizard::find


    typedef bits::vf2::AlwaysEq NodeEq;
    NodeEq _node_eq;

    Vf2WizardBase(const G1 &g1,const G2 &g2)
      : _g1(g1), _g2(g2), _mapping_type(SUBGRAPH), _local_mapping(true) { }
  };


  /// Auxiliary class for the function-type interface of %VF2 algorithm.
  ///
  /// This auxiliary class implements the named parameters of
  /// \ref vf2() "function-type interface" of \ref Vf2 algorithm.
  ///
  /// \warning This class is not to be used directly.
  ///
  /// \tparam TR The traits class that defines various types used by the
  /// algorithm.
  template<class TR>
  class Vf2Wizard : public TR {
    typedef TR Base;
    typedef typename TR::Graph1 Graph1;
    typedef typename TR::Graph2 Graph2;

    typedef typename TR::Mapping Mapping;
    typedef typename TR::NodeEq NodeEq;

    using TR::_g1;
    using TR::_g2;
    using TR::_mapping_type;
    using TR::_mapping;
    using TR::_node_eq;

  public:
    ///Constructor
    Vf2Wizard(const Graph1 &g1,const Graph2 &g2) : Base(g1,g2) {}

    ///Copy constructor
    Vf2Wizard(const Base &b) : Base(b) {}

    ///Copy constructor
    Vf2Wizard(const Vf2Wizard &b) : Base(b) {}


    template<class T>
    struct SetMappingBase : public Base{
      typedef T Mapping;
      SetMappingBase(const Base &b) : Base(b) {}
    };

    ///\brief \ref named-templ-param "Named parameter" for setting
    ///the mapping.
    ///
    ///\ref named-templ-param "Named parameter" function for setting
    ///the map that stores the found embedding.
    template<class T>
    Vf2Wizard< SetMappingBase<T> > mapping(const T &t) {
      Base::_mapping=reinterpret_cast<void*>(const_cast<T*>(&t));
      Base::_local_mapping = false;
      return Vf2Wizard<SetMappingBase<T> >(*this);
    }

    template<class NE>
    struct SetNodeEqBase : public Base {
      typedef NE NodeEq;
      NodeEq _node_eq;
      SetNodeEqBase(const Base &b, const NE &node_eq)
        : Base(b), _node_eq(node_eq){
      }
    };

    ///\brief \ref named-templ-param "Named parameter" for setting
    ///the node equivalence relation.
    ///
    ///\ref named-templ-param "Named parameter" function for setting
    ///the equivalence relation between the nodes.
    ///
    ///\param node_eq A bool-valued binary functor determinining
    ///whether a node is mappable to another. By default it is an
    ///always true operator.
    template<class T>
    Vf2Wizard< SetNodeEqBase<T> > nodeEq(const T &node_eq) {
      return Vf2Wizard<SetNodeEqBase<T> >(SetNodeEqBase<T>(*this,node_eq));
    }

    ///\brief \ref named-templ-param "Named parameter" for setting
    ///the node labels.
    ///
    ///\ref named-templ-param "Named parameter" function for setting
    ///the node labels defining equivalence relation between them.
    ///
    ///\param m1 An arbitrary \ref concepts::ReadMap "readable node map"
    ///of g1.
    ///\param m2 An arbitrary \ref concepts::ReadMap "readable node map"
    ///of g2.
    ///
    ///The value type of these maps must be equal comparable.
    template<class M1, class M2>
    Vf2Wizard< SetNodeEqBase<bits::vf2::MapEq<M1,M2> > >
    nodeLabels(const M1 &m1,const M2 &m2){
      return nodeEq(bits::vf2::MapEq<M1,M2>(m1,m2));
    }

    ///\brief \ref named-templ-param "Named parameter" for setting
    ///the mapping type.
    ///
    ///\ref named-templ-param "Named parameter" for setting
    ///the mapping type.
    ///
    ///The mapping type is set to \ref SUBGRAPH by default.
    ///
    ///\sa See \ref MappingType for the possible values.
    Vf2Wizard<Base> &mappingType(MappingType m_type) {
      _mapping_type = m_type;
      return *this;
    }

    ///\brief \ref named-templ-param "Named parameter" for setting
    ///the mapping type to \ref INDUCED.
    ///
    ///\ref named-templ-param "Named parameter" for setting
    ///the mapping type to \ref INDUCED.
    Vf2Wizard<Base> &induced() {
      _mapping_type = INDUCED;
      return *this;
    }

    ///\brief \ref named-templ-param "Named parameter" for setting
    ///the mapping type to \ref ISOMORPH.
    ///
    ///\ref named-templ-param "Named parameter" for setting
    ///the mapping type to \ref ISOMORPH.
    Vf2Wizard<Base> &iso() {
      _mapping_type = ISOMORPH;
      return *this;
    }


    ///Runs VF2 algorithm.

    ///This method runs VF2 algorithm.
    ///
    ///\retval true if a mapping is found.
    ///\retval false if there is no mapping.
    bool run(){
      if(Base::_local_mapping)
        Base::createMapping();

      Vf2<Graph1, Graph2, Mapping, NodeEq >
        alg(_g1, _g2, *reinterpret_cast<Mapping*>(_mapping), _node_eq);

      alg.mappingType(_mapping_type);

      bool ret = alg.find();

      if(Base::_local_mapping)
        delete reinterpret_cast<Mapping*>(_mapping);

      return ret;
    }

    ///Get a pointer to the generated Vf2 object.

    ///Gives a pointer to the generated Vf2 object.
    ///
    ///\return Pointer to the generated Vf2 object.
    ///\warning Don't forget to delete the referred Vf2 object after use.
    Vf2<Graph1, Graph2, Mapping, NodeEq >* getPtrToVf2Object() {
      if(Base::_local_mapping)
        Base::createMapping();
      Vf2<Graph1, Graph2, Mapping, NodeEq >* ptr =
        new Vf2<Graph1, Graph2, Mapping, NodeEq>
        (_g1, _g2, *reinterpret_cast<Mapping*>(_mapping), _node_eq);
      ptr->mappingType(_mapping_type);
      if(Base::_local_mapping)
        ptr->_deallocMappingAfterUse = true;
      return ptr;
    }

    ///Counts the number of mappings.

    ///This method counts the number of mappings.
    ///
    /// \return The number of mappings.
    int count() {
      if(Base::_local_mapping)
        Base::createMapping();

      Vf2<Graph1, Graph2, Mapping, NodeEq>
        alg(_g1, _g2, *reinterpret_cast<Mapping*>(_mapping), _node_eq);
      if(Base::_local_mapping)
        alg._deallocMappingAfterUse = true;
      alg.mappingType(_mapping_type);

      int ret = 0;
      while(alg.find())
        ++ret;

      return ret;
    }
  };

  ///Function-type interface for VF2 algorithm.

  /// \ingroup graph_isomorphism
  ///Function-type interface for VF2 algorithm \cite cordella2004sub.
  ///
  ///This function has several \ref named-func-param "named parameters"
  ///declared as the members of class \ref Vf2Wizard.
  ///The following examples show how to use these parameters.
  ///\code
  ///  // Find an embedding of graph g1 into graph g2
  ///  ListGraph::NodeMap<ListGraph::Node> m(g);
  ///  vf2(g1,g2).mapping(m).run();
  ///
  ///  // Check whether graphs g1 and g2 are isomorphic
  ///  bool is_iso = vf2(g1,g2).iso().run();
  ///
  ///  // Count the number of isomorphisms
  ///  int num_isos = vf2(g1,g2).iso().count();
  ///
  ///  // Iterate through all the induced subgraph mappings of graph g1 into g2
  ///  auto* myVf2 = vf2(g1,g2).mapping(m).nodeLabels(c1,c2)
  ///  .induced().getPtrToVf2Object();
  ///  while(myVf2->find()){
  ///    //process the current mapping m
  ///  }
  ///  delete myVf22;
  ///\endcode
  ///\warning Don't forget to put the \ref Vf2Wizard::run() "run()",
  ///\ref Vf2Wizard::count() "count()" or
  ///the \ref Vf2Wizard::getPtrToVf2Object() "getPtrToVf2Object()"
  ///to the end of the expression.
  ///\sa Vf2Wizard
  ///\sa Vf2
  template<class G1, class G2>
  Vf2Wizard<Vf2WizardBase<G1,G2> > vf2(const G1 &g1, const G2 &g2) {
    return Vf2Wizard<Vf2WizardBase<G1,G2> >(g1,g2);
  }

}

#endif
