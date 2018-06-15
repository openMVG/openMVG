## Copyright 2004,2007-2010,2012,2015 IPB, Universite de Bordeaux, INRIA & CNRS
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

bindir		= ../../bin
includedir	= ../../include
libdir		= ../../lib

##
##  General inference rules.
##

include ../Makefile.inc

%$(OBJ)	:	%.c
				$(CC) $(CFLAGS) $(CLIBFLAGS) -I$(includedir) -c $(<) -o $(@)

%$(EXE)	:	%.c
		 		$(CC) $(CFLAGS) -I$(includedir) $(<) -o $(@) -L$(libdir) $(LDFLAGS) -L. -l$(ESMUMPSLIB) -l$(SCOTCHLIB) -lscotch -l$(SCOTCHLIB)errexit -lm

##
##  Project rules.
##

.PHONY				:	ptscotch	scotch	ptinstall	install	clean	realclean

scotch				:	clean
					$(MAKE) CC="$(CCS)" CCD="$(CCS)" SCOTCHLIB=scotch ESMUMPSLIB=esmumps	\
					libesmumps$(LIB)							\
					main_esmumps$(EXE)

ptscotch			:	clean
					$(MAKE) CFLAGS="$(CFLAGS) -DSCOTCH_PTSCOTCH" CC="$(CCP)" SCOTCHLIB=ptscotch ESMUMPSLIB=ptesmumps	\
					libesmumps$(LIB)											\
					main_esmumps$(EXE)

install				:
					-$(CP) esmumps.h $(includedir)
					-$(CP) libesmumps$(LIB) $(libdir)

ptinstall			:
					-$(CP) esmumps.h $(includedir)
					-$(CP) libptesmumps$(LIB) $(libdir)

clean				:
					-$(RM) *~ common.h *$(OBJ) lib*$(LIB) main_esmumps$(EXE)

realclean			:	clean

##
##  Todo list.
##

common.h			:	../libscotch/common.h			\
					module.h
					$(CAT) module.h ../libscotch/common.h > $(@)

graph_graph$(OBJ)		:	graph_graph.c				\
					common.h				\
					$(includedir)/$(SCOTCHLIB).h		\
					graph.h

order$(OBJ)			:	order.c					\
					order.h					\
					common.h				\
					graph.h

order_scotch_graph$(OBJ)	:	order_scotch_graph.c			\
					order.h					\
					common.h				\
					graph.h					\
					$(includedir)/$(SCOTCHLIB).h

dof$(OBJ)			:	dof.c					\
					dof.h					\
					common.h				\
					graph.h					\
					$(includedir)/$(SCOTCHLIB).h

symbol$(OBJ)			:	symbol.c				\
					symbol.h				\
					common.h

symbol_fax$(OBJ)		:	symbol_fax.c				\
					common.h				\
					graph.h					\
					symbol.h				\
					order.h					\
					fax.h

symbol_fax_graph$(OBJ)		:	symbol_fax_graph.c			\
					symbol_fax.c				\
					common.h				\
					graph.h					\
					symbol.h				\
					order.h					\
					fax.h					\
					$(includedir)/$(SCOTCHLIB).h

esmumps$(OBJ)			:	esmumps.c				\
					common.h				\
					graph.h					\
					symbol.h				\
					order.h					\
					fax.h					\
					esmumps.h				\
					$(includedir)/$(SCOTCHLIB).h

esmumps_f$(OBJ)			:	esmumps_f.c				\
					common.h				\
					esmumps.h

esmumps_strats$(OBJ)		:	esmumps_strats.c			\
					common.h				\
					esmumps.h

libesmumps$(LIB)		:	graph_graph$(OBJ)			\
					order$(OBJ)				\
					order_scotch_graph$(OBJ)		\
					dof$(OBJ)				\
					symbol$(OBJ)				\
					symbol_fax_graph$(OBJ)			\
					esmumps$(OBJ)				\
					esmumps_f$(OBJ)				\
					esmumps_strats$(OBJ)
					$(AR) $(ARFLAGS) lib$(ESMUMPSLIB)$(LIB) $(?)
					-$(RANLIB) lib$(ESMUMPSLIB)$(LIB)

main_esmumps$(EXE)		:	main_esmumps.c				\
					common.h				\
					graph.h					\
					order.h					\
					symbol.h				\
					esmumps.h				\
					lib$(ESMUMPSLIB)$(LIB)			\
					$(libdir)/lib$(SCOTCHLIB)$(LIB)		\
					$(libdir)/lib$(SCOTCHLIB)errexit$(LIB)
