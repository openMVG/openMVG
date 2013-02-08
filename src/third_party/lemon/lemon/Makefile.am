EXTRA_DIST += \
	lemon/lemon.pc.in \
	lemon/lemon.pc.cmake \
	lemon/CMakeLists.txt \
	lemon/config.h.cmake

pkgconfig_DATA += lemon/lemon.pc

lib_LTLIBRARIES += lemon/libemon.la

lemon_libemon_la_SOURCES = \
	lemon/arg_parser.cc \
	lemon/base.cc \
	lemon/color.cc \
	lemon/lp_base.cc \
	lemon/lp_skeleton.cc \
	lemon/random.cc \
	lemon/bits/windows.cc

nodist_lemon_HEADERS = lemon/config.h	

lemon_libemon_la_CXXFLAGS = \
	$(AM_CXXFLAGS) \
	$(GLPK_CFLAGS) \
	$(CPLEX_CFLAGS) \
	$(SOPLEX_CXXFLAGS) \
	$(CLP_CXXFLAGS) \
	$(CBC_CXXFLAGS)

lemon_libemon_la_LDFLAGS = \
	$(GLPK_LIBS) \
	$(CPLEX_LIBS) \
	$(SOPLEX_LIBS) \
	$(CLP_LIBS) \
	$(CBC_LIBS)

if HAVE_GLPK
lemon_libemon_la_SOURCES += lemon/glpk.cc
endif

if HAVE_CPLEX
lemon_libemon_la_SOURCES += lemon/cplex.cc
endif

if HAVE_SOPLEX
lemon_libemon_la_SOURCES += lemon/soplex.cc
endif

if HAVE_CLP
lemon_libemon_la_SOURCES += lemon/clp.cc
endif

if HAVE_CBC
lemon_libemon_la_SOURCES += lemon/cbc.cc
endif

lemon_HEADERS += \
	lemon/adaptors.h \
	lemon/arg_parser.h \
	lemon/assert.h \
	lemon/bellman_ford.h \
	lemon/bfs.h \
	lemon/bin_heap.h \
	lemon/binomial_heap.h \
	lemon/bucket_heap.h \
	lemon/capacity_scaling.h \
	lemon/cbc.h \
	lemon/circulation.h \
	lemon/clp.h \
	lemon/color.h \
	lemon/concept_check.h \
	lemon/connectivity.h \
	lemon/core.h \
	lemon/cost_scaling.h \
	lemon/counter.h \
	lemon/cplex.h \
	lemon/cycle_canceling.h \
	lemon/dfs.h \
	lemon/dheap.h \
	lemon/dijkstra.h \
	lemon/dim2.h \
	lemon/dimacs.h \
	lemon/edge_set.h \
	lemon/elevator.h \
	lemon/error.h \
	lemon/euler.h \
	lemon/fib_heap.h \
	lemon/fractional_matching.h \
	lemon/full_graph.h \
	lemon/glpk.h \
	lemon/gomory_hu.h \
	lemon/graph_to_eps.h \
	lemon/grid_graph.h \
	lemon/hartmann_orlin_mmc.h \
	lemon/howard_mmc.h \
	lemon/hypercube_graph.h \
	lemon/karp_mmc.h \
	lemon/kruskal.h \
	lemon/hao_orlin.h \
	lemon/lgf_reader.h \
	lemon/lgf_writer.h \
	lemon/list_graph.h \
	lemon/lp.h \
	lemon/lp_base.h \
	lemon/lp_skeleton.h \
	lemon/maps.h \
	lemon/matching.h \
	lemon/math.h \
	lemon/min_cost_arborescence.h \
	lemon/nauty_reader.h \
	lemon/network_simplex.h \
	lemon/pairing_heap.h \
	lemon/path.h \
	lemon/planarity.h \
	lemon/preflow.h \
	lemon/quad_heap.h \
	lemon/radix_heap.h \
	lemon/radix_sort.h \
	lemon/random.h \
	lemon/smart_graph.h \
	lemon/soplex.h \
	lemon/static_graph.h \
	lemon/suurballe.h \
	lemon/time_measure.h \
	lemon/tolerance.h \
	lemon/unionfind.h \
	lemon/bits/windows.h

bits_HEADERS += \
	lemon/bits/alteration_notifier.h \
	lemon/bits/array_map.h \
	lemon/bits/bezier.h \
	lemon/bits/default_map.h \
	lemon/bits/edge_set_extender.h \
	lemon/bits/enable_if.h \
	lemon/bits/graph_adaptor_extender.h \
	lemon/bits/graph_extender.h \
	lemon/bits/map_extender.h \
	lemon/bits/path_dump.h \
	lemon/bits/solver_bits.h \
	lemon/bits/traits.h \
	lemon/bits/variant.h \
	lemon/bits/vector_map.h

concept_HEADERS += \
	lemon/concepts/digraph.h \
	lemon/concepts/graph.h \
	lemon/concepts/graph_components.h \
	lemon/concepts/heap.h \
	lemon/concepts/maps.h \
	lemon/concepts/path.h
