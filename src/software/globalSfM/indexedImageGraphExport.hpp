
// Copyright (c) 2012, 2013, 2014 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_INDEXED_IMAGE_GRAPH_EXPORT_H_
#define OPENMVG_INDEXED_IMAGE_GRAPH_EXPORT_H_

#include <lemon/core.h>

#include <iostream>
#include <fstream>

namespace openMVG  {
namespace imageGraph  {
  using namespace std;

// Export an Image connection graph
//   to graphviz file format.
// Example :
// graph 1 {
//   node [shape=circle]
//   n2
//   n1
//   n0
//   n1 --  n0
//   n2 --  n0
//   n2 --  n1
// }
template <typename GraphT>
bool exportToGraphvizFormat_Nodal(
  const GraphT & g,
  ostream & os)
{
  os << "graph 1 {" << endl;
  os << "node [shape=circle]" << endl;
  //Export node label
  for(typename GraphT::NodeIt n(g); n!= lemon::INVALID; ++n)
  {
    os << "  n" << g.id(n) << std::endl;
  }

  //-- Export arc (as the graph is bi-directional, export arc only one time)

  map< std::pair<size_t,size_t>, size_t > map_arcs;
  for(typename GraphT::ArcIt e(g); e!=lemon::INVALID; ++e) {
    if( map_arcs.end() == map_arcs.find(std::make_pair(size_t (g.id(g.source(e))), size_t (g.id(g.target(e)))))
      &&
      map_arcs.end() == map_arcs.find(std::make_pair(size_t (g.id(g.target(e))), size_t (g.id(g.source(e))))))
    {
      map_arcs[std::pair<size_t,size_t>(size_t (g.id(g.source(e))),
        size_t (g.id(g.target(e)))) ] = 1.0;
    }
  }
  //os << "edge [style=bold]" << endl;
  for ( map< std::pair<size_t,size_t>, size_t >::const_iterator iter = map_arcs.begin();
    iter != map_arcs.end();
    ++iter)
  {
    os << "  n" << iter->first.first << " -- " << " n" << iter->first.second << endl;
  }

  os << "}" << endl;
  return os.good();
}

// Export the Image connection graph
//   to the graphviz file format.
// Add the image name and preview in the graph by using HTML label.
template <typename GraphT, typename NodeMap, typename EdgeMap>
bool exportToGraphvizFormat_Image(
  const GraphT & g,
  const NodeMap & nodeMap,
  const EdgeMap & edgeMap,
  ostream & os, bool bWeightedEdge=false)
{
  os << "graph 1 {" << endl;
  os << "node [shape=none]" << endl;
  //Export node label
  for(typename GraphT::NodeIt n(g); n!=lemon::INVALID; ++n)
  {
    os << "  n" << g.id(n)
       << "[ label ="
      <<
      "< "<< endl
      <<"<table>"<< endl
      <<"<tr><td>" << "\"" << nodeMap[n] <<"\"" <<"</td></tr>"<< endl
      <<"<tr><td><img src=\"" << nodeMap[n] <<"\"/></td></tr>"<< endl
      <<"</table>"<< endl
      <<">, cluster=1];"<< endl;

    //os << "  n" << g.id(n)
    //  << " [ "
    //  << " image=\"" << nodeMap[n] << "\" cluster=1]; " << endl;
  }

  //Export arc value
  map< std::pair<size_t,size_t>, size_t > map_arcs;
  for(typename GraphT::ArcIt e(g); e!=lemon::INVALID; ++e) {
    if( map_arcs.end() == map_arcs.find(std::make_pair(size_t (g.id(g.source(e))), size_t (g.id(g.target(e)))))
      &&
      map_arcs.end() == map_arcs.find(std::make_pair(size_t (g.id(g.target(e))), size_t (g.id(g.source(e))))))
    {
      map_arcs[std::pair<size_t,size_t>(size_t (g.id(g.source(e))),
        size_t (g.id(g.target(e)))) ] = edgeMap[e];
    }
  }

  os << "edge [style=bold]" << endl;
  for ( map< std::pair<size_t,size_t>, size_t >::const_iterator iter = map_arcs.begin();
    iter != map_arcs.end();
    ++iter)
  {
    if (bWeightedEdge)
    {
      os << "  n" << iter->first.first << " -- " << " n" << iter->first.second
        << " [label=\"" << iter->second << "\"]" << endl;
    }
    else
    {
      os << "  n" << iter->first.first << " -- " << " n" << iter->first.second << endl;
    }
  }
  os << "}" << endl;
  return os.good();
}

template <typename GraphT>
void exportToGraphvizData(const std::string& sfile, const GraphT & graph){
  //Prepare Data

  std::ofstream file(sfile.c_str());
  openMVG::imageGraph::exportToGraphvizFormat_Nodal(graph, file);
  file.close();

  //Use Graphviz
  const std::string cmd = "neato -Tsvg -O -Goverlap=scale -Gsplines=false " + sfile;
  system(cmd.c_str());
}

} // namespace imageGraph
} // namespace openMVG
#endif // OPENMVG_INDEXED_IMAGE_GRAPH_EXPORT_H_
