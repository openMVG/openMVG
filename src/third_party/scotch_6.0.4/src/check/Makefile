## Copyright 2011,2012,2014,2015 IPB, Universite de Bordeaux, INRIA & CNRS
##
## This file is part of the Scotch software package for static mapping,
## graph partitioning and sparse matrix ordering.
##
## This software is governed by the CeCILL-C license under French law
## and abiding by the rules of distribution of free software. You can
## use, modify and/or redistribute the software under the terms of the
## CeCILL-C license as circulated by CEA, CNRS and INRIA at the following
## URL: "http://www.cecill.info".
## 
## As a counterpart to the access to the source code and rights to copy,
## modify and redistribute granted by the license, users are provided
## only with a limited warranty and the software's author, the holder of
## the economic rights, and the successive licensors have only limited
## liability.
## 
## In this respect, the user's attention is drawn to the risks associated
## with loading, using, modifying and/or developing or reproducing the
## software by the user in light of its specific status of free software,
## that may mean that it is complicated to manipulate, and that also
## therefore means that it is reserved for developers and experienced
## professionals having in-depth computer knowledge. Users are therefore
## encouraged to load and test the software's suitability as regards
## their requirements in conditions enabling the security of their
## systems and/or data to be ensured and, more generally, to use and
## operate it in the same conditions as regards security.
## 
## The fact that you are presently reading this means that you have had
## knowledge of the CeCILL-C license and that you accept its terms.
##

SCOTCHINCLUDEDIR	= ../../include
SCOTCHLIBDIR		= ../../lib
SCOTCHBINDIR		= ../../bin

LIBPTSCOTCH		= -lptscotch -lscotch -lptscotcherrexit
LIBSCOTCH		= -lscotch -lscotcherrexit

DEPPTSCOTCH		= $(SCOTCHLIBDIR)/libptscotch$(LIB)
DEPSCOTCH		= $(SCOTCHLIBDIR)/libscotch$(LIB)

EXECS			=
EXECP			= mpirun -n 4

##
##  General inference rules.
##

include ../Makefile.inc

%$(OBJ)	:	%.c
				$(CC) $(CFLAGS) $(CLIBFLAGS) -I$(SCOTCHINCLUDEDIR) -c $(<) -o $(@)

%$(EXE)	:	%.c
		 		$(CC) $(CFLAGS) -I$(SCOTCHINCLUDEDIR) -L$(SCOTCHLIBDIR) $(<) -o $(@) $(SCOTCHLIBS) $(LDFLAGS)

##
##  Project rules.
##

.PHONY				:	check	ptcheck	clean	realclean

check				:	realclean
					$(MAKE) CC="$(CCS)" CCD="$(CCS)" SCOTCHLIBS="$(LIBSCOTCH)" $(CHECKSCOTCH)

ptcheck				:	realclean
					$(MAKE) CFLAGS="$(CFLAGS) -DSCOTCH_PTSCOTCH" CC="$(CCP)" SCOTCHLIBS="$(LIBPTSCOTCH)" $(CHECKPTSCOTCH)

clean				:
					-$(RM) *~ *$(OBJ)

realclean			:	clean
					-$(RM) 				\
					test_common_random		\
					test_common_thread		\
					test_scotch_graph_coarsen	\
					test_scotch_graph_coarsen_build	\
					test_scotch_graph_color		\
					test_scotch_graph_map		\
					test_scotch_graph_map_copy	\
					test_scotch_graph_order		\
					test_scotch_graph_part_ovl	\
					test_scotch_dgraph_band		\
					test_scotch_dgraph_check	\
					test_scotch_dgraph_coarsen	\
					test_scotch_dgraph_grow		\
					test_scotch_dgraph_redist	\
					test_strat_par			\
					test_strat_seq

##
##  Test cases dependencies.
##

CHECKSCOTCH = 		check_common_random			\
			check_common_thread			\
			check_strat_seq				\
			check_scotch_graph_coarsen		\
			check_scotch_graph_coarsen_build	\
			check_scotch_graph_color		\
			check_scotch_graph_map			\
			check_scotch_graph_map_copy		\
			check_scotch_graph_order		\
			check_scotch_graph_part_ovl		\
			check_prog_gord				\
			check_prog_gpart_remap

CHECKPTSCOTCH = 	check_strat_par				\
			check_scotch_dgraph_check		\
			check_scotch_dgraph_band		\
			check_scotch_dgraph_coarsen		\
			check_scotch_dgraph_grow		\
			check_scotch_dgraph_redist		\

##
##  Todo list.
##

check_common_random		:	test_common_random
					$(EXECS) ./test_common_random /tmp/rand.dat 0
					$(EXECS) ./test_common_random /tmp/rand.dat 1

test_common_random		:	test_common_random.c		\
					$(SCOTCHLIBDIR)/libscotch$(LIB)

##

check_common_thread		:	test_common_thread
					$(EXECS) ./test_common_thread

test_common_thread		:	test_common_thread.c		\
					$(SCOTCHLIBDIR)/libscotch$(LIB)

##

check_dgraph_fold_comm		:	test_dgraph_fold_comm
					$(EXECS) ./test_dgraph_fold_comm data/dgraph_fold_comm_1.txt

test_dgraph_fold_comm		:	test_dgraph_fold_comm.c

##

check_scotch_dgraph_redist	:	test_scotch_dgraph_redist
					$(EXECP) ./test_scotch_dgraph_redist data/bump.grf

test_scotch_dgraph_redist	:	test_scotch_dgraph_redist.c

##

check_prog_gord			:
					$(EXECS) $(SCOTCHBINDIR)/gord data/bump.grf /dev/null -vt
					$(EXECS) $(SCOTCHBINDIR)/gord data/bump_b1.grf /dev/null -vt

##

check_prog_gpart_remap		:
					$(EXECS) $(SCOTCHBINDIR)/gpart 32 data/bump_imbal_32.grf /dev/null -rodata/bump_old.map -vmt

##

check_scotch_graph_coarsen	:	test_scotch_graph_coarsen
					$(EXECS) ./test_scotch_graph_coarsen data/bump_b1.grf

test_scotch_graph_coarsen	:	test_scotch_graph_coarsen.c	\
					$(SCOTCHLIBDIR)/libscotch$(LIB)

##

check_scotch_graph_coarsen_build :	test_scotch_graph_coarsen_build
					$(EXECS) ./test_scotch_graph_coarsen_build data/bump.grf
					$(EXECS) ./test_scotch_graph_coarsen_build data/bump_b1.grf

test_scotch_graph_coarsen_build	:	test_scotch_graph_coarsen_build.c \
					$(SCOTCHLIBDIR)/libscotch$(LIB)

##

check_scotch_graph_color	:	test_scotch_graph_color
					$(EXECS) ./test_scotch_graph_color data/bump.grf

test_scotch_graph_color		:	test_scotch_graph_color.c	\
					$(SCOTCHLIBDIR)/libscotch$(LIB)

##

check_scotch_graph_map		:	test_scotch_graph_map
					$(EXECS) ./test_scotch_graph_map data/m4x4.grf
					$(EXECS) ./test_scotch_graph_map data/m4x4_b1.grf
					$(EXECS) ./test_scotch_graph_map data/m16x16.grf
					$(EXECS) ./test_scotch_graph_map data/m16x16_b1.grf

test_scotch_graph_map		:	test_scotch_graph_map.c		\
					$(SCOTCHLIBDIR)/libscotch$(LIB)

##

check_scotch_graph_map_copy	:	test_scotch_graph_map_copy
					$(EXECS) ./test_scotch_graph_map_copy data/bump_b1.grf

test_scotch_graph_map_copy	:	test_scotch_graph_map_copy.c	\
					$(SCOTCHLIBDIR)/libscotch$(LIB)

##

check_scotch_graph_order	:	test_scotch_graph_order
					$(EXECS) ./test_scotch_graph_order data/bump.grf
					$(EXECS) ./test_scotch_graph_order data/bump_b1.grf

test_scotch_graph_order		:	test_scotch_graph_order.c	\
					$(SCOTCHLIBDIR)/libscotch$(LIB)

##

check_scotch_graph_part_ovl	:	test_scotch_graph_part_ovl
					$(EXECS) ./test_scotch_graph_part_ovl 4 data/m16x16.grf /dev/null
					$(EXECS) ./test_scotch_graph_part_ovl 4 data/m16x16_b1.grf /dev/null

test_scotch_graph_part_ovl	:	test_scotch_graph_part_ovl.c	\
					$(SCOTCHLIBDIR)/libscotch$(LIB)

##

check_strat_seq			:	test_strat_seq
					$(EXECS) ./test_strat_seq

test_strat_seq			:	test_strat_seq.c		\
					$(SCOTCHLIBDIR)/libscotch$(LIB)

##

check_strat_par			:	test_strat_par
					$(EXECS) ./test_strat_par

test_strat_par			:	test_strat_par.c		\
					$(SCOTCHLIBDIR)/libptscotch$(LIB)

##

check_scotch_dgraph_band	:	test_scotch_dgraph_band
					$(EXECP) ./test_scotch_dgraph_band data/bump_b1.grf

test_scotch_dgraph_band		:	test_scotch_dgraph_band.c	\
					$(SCOTCHLIBDIR)/libptscotch$(LIB)

##

check_scotch_dgraph_coarsen	:	test_scotch_dgraph_coarsen
					$(EXECP) ./test_scotch_dgraph_coarsen data/bump.grf
					$(EXECP) ./test_scotch_dgraph_coarsen data/bump_b1.grf

test_scotch_dgraph_coarsen	:	test_scotch_dgraph_coarsen.c	\
					$(SCOTCHLIBDIR)/libptscotch$(LIB)

##

check_scotch_dgraph_check	:	test_scotch_dgraph_check
					$(EXECP) ./test_scotch_dgraph_check data/bump.grf
					$(EXECP) ./test_scotch_dgraph_check data/bump_b1.grf

test_scotch_dgraph_check	:	test_scotch_dgraph_check.c	\
					$(SCOTCHLIBDIR)/libptscotch$(LIB)

##

check_scotch_dgraph_grow	:	test_scotch_dgraph_grow
					$(EXECP) ./test_scotch_dgraph_grow data/bump.grf

test_scotch_dgraph_grow		:	test_scotch_dgraph_grow.c	\
					$(SCOTCHLIBDIR)/libptscotch$(LIB)
