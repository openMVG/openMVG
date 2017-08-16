## Copyright 2007-2010 ENSEIRB, INRIA & CNRS
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
				$(CC) $(CFLAGS) $(CLIBFLAGS) -I$(includedir) -I../libscotch -c $(<) -o $(@)

%$(EXE)	:	%.c
		 		$(CC) $(CFLAGS) -I$(includedir) -I../libscotch $(<) -o $(@) -L$(libdir) $(LDFLAGS)

##
##  Project rules.
##

.PHONY				:	ptscotch	scotch	ptinstall	install	clean	realclean

scotch				:
					$(MAKE) CC="$(CCS)" SCOTCHLIB=ptscotch						\
					libscotchmetis$(LIB)

ptscotch			:
					$(MAKE) CFLAGS="$(CFLAGS) -DSCOTCH_PTSCOTCH" CC="$(CCP)" SCOTCHLIB=ptscotch	\
					libptscotchparmetis$(LIB)

install				:	scotch
					-$(CP) metis.h $(includedir)
					-$(CP) libscotchmetis$(LIB) $(libdir)

ptinstall			:	ptscotch
					-$(CP) parmetis.h $(includedir)
					-$(CP) libptscotchparmetis$(LIB) $(libdir)

clean				:
					-$(RM) *~ *$(OBJ) lib*$(LIB)

realclean			:	clean

##
##  Todo list.
##

metis_graph_order$(OBJ)		:	metis_graph_order.c			\
					../libscotch/module.h			\
					../libscotch/common.h			\
					$(includedir)/scotch.h			\
					metis.h

metis_graph_order_f$(OBJ)	:	metis_graph_order_f.c			\
					../libscotch/module.h			\
					../libscotch/common.h			\
					$(includedir)/scotch.h			\
					metis.h

metis_graph_part$(OBJ)		:	metis_graph_part.c			\
					../libscotch/module.h			\
					../libscotch/common.h			\
					$(includedir)/scotch.h			\
					metis.h

metis_graph_part_f$(OBJ)	:	metis_graph_part_f.c			\
					../libscotch/module.h			\
					../libscotch/common.h			\
					$(includedir)/scotch.h			\
					metis.h

parmetis_dgraph_order$(OBJ)	:	parmetis_dgraph_order.c			\
					../libscotch/module.h			\
					../libscotch/common.h			\
					$(includedir)/ptscotch.h		\
					parmetis.h

parmetis_dgraph_order_f$(OBJ)	:	parmetis_dgraph_order_f.c		\
					../libscotch/module.h			\
					../libscotch/common.h			\
					$(includedir)/ptscotch.h		\
					parmetis.h

parmetis_dgraph_part$(OBJ)	:	parmetis_dgraph_part.c			\
					../libscotch/module.h			\
					../libscotch/common.h			\
					$(includedir)/ptscotch.h		\
					parmetis.h

parmetis_dgraph_part_f$(OBJ)	:	parmetis_dgraph_part_f.c		\
					../libscotch/module.h			\
					../libscotch/common.h			\
					$(includedir)/ptscotch.h		\
					parmetis.h

libptscotchparmetis$(LIB)	:	parmetis_dgraph_order$(OBJ)		\
					parmetis_dgraph_order_f$(OBJ)		\
					parmetis_dgraph_part$(OBJ)		\
					parmetis_dgraph_part_f$(OBJ)
					$(AR) $(ARFLAGS) $(@) $(^)
					-$(RANLIB) $(@)

libscotchmetis$(LIB)		:	metis_graph_order$(OBJ)			\
					metis_graph_order_f$(OBJ)		\
					metis_graph_part$(OBJ)			\
					metis_graph_part_f$(OBJ)
					$(AR) $(ARFLAGS) $(@) $(^)
					-$(RANLIB) $(@)
