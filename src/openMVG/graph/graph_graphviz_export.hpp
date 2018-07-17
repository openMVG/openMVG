// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2012, 2013, 2014 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_GRAPH_GRAPH_EXPORT_HPP
#define OPENMVG_GRAPH_GRAPH_EXPORT_HPP

#include <lemon/list_graph.h>

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <string>

#if __APPLE__
#include "TargetConditionals.h"
#endif

#include "openMVG/types.hpp"

namespace openMVG {
namespace graph {

    /**
    * @brief Export an Image connection graph
    *   to graphviz file format.
    * Example :
    * graph 1 {
    *   node [shape=circle]
    *   n1 --  n0
    *   n2 --  n0
    *   n2 --  n1
    * }
    * @param g Graph to export
    * @param os stream where graph is exported
    * @retval true If export is correct
    * @retval false If export is incorrect
    */
inline bool exportToGraphvizFormat_Nodal
(
  const indexedGraph & graph,
  std::ostream & os
)
{
  os << "graph 1 {" << std::endl;
  os << "node [shape=circle]" << std::endl;

  //-- Export graph edges (just a link between the nodes)
  for (indexedGraph::GraphT::EdgeIt e(graph.g); e != lemon::INVALID; ++e)
  {
    os
      << " n" << (*graph.node_map_id)[graph.g.u(e)]
      << " -- "
      << " n" << (*graph.node_map_id)[graph.g.v(e)] << '\n';
  }

  os << "}" << std::endl;
  return os.good();
}

  /**
  * @brief Export a graph and generate it using graphviz
  * @param sfile File in which graph is exported
  * @param graph Graph to export
  */
inline void exportToGraphvizData
(
  const std::string& sfile,
  const indexedGraph & graph
)
{
  // Export the graph as a DOT (graph description language) file
  std::ofstream file(sfile.c_str());
  openMVG::graph::exportToGraphvizFormat_Nodal(graph, file);
  file.close();

  //Use Graphviz
  const std::string cmd = "neato -Tsvg -O -Goverlap=scale -Gsplines=false " + sfile;
#ifndef TARGET_OS_IPHONE
  const int ret = std::system(cmd.c_str());
  (void)ret;
#endif
}

} // namespace graph
} // namespace openMVG

#endif // OPENMVG_GRAPH_GRAPH_EXPORT_HPP
