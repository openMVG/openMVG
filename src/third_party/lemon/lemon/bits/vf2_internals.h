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

#ifndef VF2_INTERNALS_H
#define VF2_INTERNALS_H


///\ingroup graph_properties
///\file
///\brief Mapping types for graph matching algorithms.

namespace lemon {
  ///\ingroup graph_isomorphism
  ///The \ref Vf2 "VF2" algorithm is capable of finding different kind of
  ///graph embeddings, this enum specifies their types.
  ///
  ///See \ref graph_isomorphism for a more detailed description.
  enum MappingType {
    /// Subgraph isomorphism
    SUBGRAPH = 0,
    /// Induced subgraph isomorphism
    INDUCED = 1,
    /// Graph isomorphism
    ///
    /// If the two graphs have the same number of nodes, than it is
    /// equivalent to \ref INDUCED, and if they also have the same
    /// number of edges, then it is also equivalent to \ref SUBGRAPH.
    ///
    /// However, using this setting is faster than the other two options.
    ISOMORPH = 2
  };
}
#endif
