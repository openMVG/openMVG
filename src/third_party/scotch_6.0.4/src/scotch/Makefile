## Copyright 2004,2007-2012 IPB, Universite de Bordeaux, INRIA & CNRS
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

%$(EXE)	:	%$(OBJ)
		 		$(CC) $(CFLAGS) -I$(includedir) -I../libscotch $(<) -o $(@) -L$(libdir) -l$(SCOTCHLIB) -lscotch -l$(SCOTCHLIB)errexit $(LDFLAGS)

%$(OBJ)	:	%.c
				$(CC) $(CFLAGS) -I$(includedir) -I../libscotch -c $(<) -o $(@)

%$(EXE)	:	%.c
		 		$(CC) $(CFLAGS) -I$(includedir) -I../libscotch $(<) -o $(@) -L$(libdir) -l$(SCOTCHLIB) -lscotch -l$(SCOTCHLIB)errexit $(LDFLAGS)

##
##  Project rules.
##

.PHONY				:	ptscotch	scotch	ptinstall	install	clean	realclean

scotch				:
					$(MAKE) CC="$(CCS)" SCOTCHLIB=scotch	\
					acpl$(EXE)				\
					amk_ccc$(EXE)				\
					amk_fft2$(EXE)				\
					amk_grf$(EXE)				\
					amk_hy$(EXE)				\
					amk_m2$(EXE)				\
					amk_p2$(EXE)				\
					atst$(EXE)				\
					gbase$(EXE)				\
					gcv$(EXE)				\
					gmap$(EXE)				\
					gmk_hy$(EXE)				\
					gmk_m2$(EXE)				\
					gmk_m3$(EXE)				\
					gmk_msh$(EXE)				\
					gmk_ub2$(EXE)				\
					gmtst$(EXE)				\
					gord$(EXE)				\
					gotst$(EXE)				\
					gout$(EXE)				\
					gpart$(EXE)				\
					gscat$(EXE)				\
					gtst$(EXE)				\
					mcv$(EXE)				\
					mmk_m2$(EXE)				\
					mmk_m3$(EXE)				\
					mord$(EXE)				\
					mtst$(EXE)

ptscotch			:
					$(MAKE) CC="$(CCP)" SCOTCHLIB=ptscotch	\
					dggath$(EXE)				\
					dgmap$(EXE)				\
					dgord$(EXE)				\
					dgpart$(EXE)				\
					dgscat$(EXE)				\
					dgtst$(EXE)

install				:	scotch
					-$(CP) acpl$(EXE) amk_ccc$(EXE) amk_fft2$(EXE) amk_grf$(EXE) amk_hy$(EXE) amk_m2$(EXE) amk_p2$(EXE) atst$(EXE) gbase$(EXE) gcv$(EXE) gmap$(EXE) gmk_hy$(EXE) gmk_m2$(EXE) gmk_m3$(EXE) gmk_msh$(EXE) gmk_ub2$(EXE) gmtst$(EXE) gord$(EXE) gotst$(EXE) gout$(EXE) gpart$(EXE) *gtst$(EXE) gscat$(EXE) mcv$(EXE) mmk_m2$(EXE) mmk_m3$(EXE) mord$(EXE) mtst$(EXE) $(bindir)

ptinstall			:	ptscotch
					-$(CP) dggath$(EXE) dgmap$(EXE) dgord$(EXE) dgpart$(EXE) dgscat$(EXE) dgtst$(EXE) $(bindir)

clean				:
					-$(RM) *~ *$(OBJ) acpl$(EXE) amk_ccc$(EXE) amk_fft2$(EXE) amk_grf$(EXE) amk_hy$(EXE) amk_m2$(EXE) amk_p2$(EXE) atst$(EXE) gbase$(EXE) gcv$(EXE) *ggath$(EXE) *gmap$(EXE) gmk_hy$(EXE) gmk_m2$(EXE) gmk_m3$(EXE) gmk_msh$(EXE) gmk_ub2$(EXE) gmtst$(EXE) *gord$(EXE) gotst$(EXE) gout$(EXE) *gpart$(EXE) *gscat$(EXE) *gtst$(EXE) mcv$(EXE) mmk_m2$(EXE) mmk_m3$(EXE) mord$(EXE) mtst$(EXE)

realclean			:	clean

##
##  Todo list.
##

acpl$(EXE)			:	acpl.c					\
					../libscotch/module.h			\
					../libscotch/common.h			\
					$(includedir)/scotch.h			\
					$(libdir)/libscotch$(LIB)		\
					$(libdir)/libscotcherrexit$(LIB)	\
					acpl.h

amk_ccc$(EXE)			:	amk_ccc.c				\
					../libscotch/module.h			\
					../libscotch/common.h			\
					$(includedir)/scotch.h			\
					$(libdir)/libscotch$(LIB)		\
					$(libdir)/libscotcherrexit$(LIB)	\
					amk_ccc.h

amk_fft2$(EXE)			:	amk_fft2.c				\
					../libscotch/module.h			\
					../libscotch/common.h			\
					$(includedir)/scotch.h			\
					$(libdir)/libscotch$(LIB)		\
					$(libdir)/libscotcherrexit$(LIB)	\
					amk_fft2.h

amk_grf$(EXE)			:	amk_grf.c				\
					../libscotch/module.h			\
					../libscotch/common.h			\
					$(includedir)/scotch.h			\
					$(libdir)/libscotch$(LIB)		\
					$(libdir)/libscotcherrexit$(LIB)	\
					amk_grf.h

amk_hy$(EXE)			:	amk_hy.c				\
					../libscotch/module.h			\
					../libscotch/common.h			\
					$(includedir)/scotch.h			\
					$(libdir)/libscotch$(LIB)		\
					$(libdir)/libscotcherrexit$(LIB)	\
					amk_hy.h

amk_m2$(EXE)			:	amk_m2.c				\
					../libscotch/module.h			\
					../libscotch/common.h			\
					$(includedir)/scotch.h			\
					$(libdir)/libscotch$(LIB)		\
					$(libdir)/libscotcherrexit$(LIB)	\
					amk_m2.h

amk_p2$(EXE)			:	amk_p2.c				\
					../libscotch/module.h			\
					../libscotch/common.h			\
					$(includedir)/scotch.h			\
					$(libdir)/libscotch$(LIB)		\
					$(libdir)/libscotcherrexit$(LIB)	\
					amk_p2.h

atst$(EXE)			:	atst.c					\
					../libscotch/module.h			\
					../libscotch/common.h			\
					$(includedir)/scotch.h			\
					$(libdir)/libscotch$(LIB)		\
					$(libdir)/libscotcherrexit$(LIB)	\
					atst.h

dggath$(EXE)			:	dggath.c				\
					../libscotch/module.h			\
					../libscotch/common.h			\
					$(includedir)/ptscotch.h		\
					$(libdir)/libscotch$(LIB)		\
					$(libdir)/libptscotch$(LIB)		\
					$(libdir)/libptscotcherrexit$(LIB)	\
					dggath.h

dgmap$(EXE)			:	dgmap.c					\
					../libscotch/module.h			\
					../libscotch/common.h			\
					$(includedir)/ptscotch.h		\
					$(libdir)/libscotch$(LIB)		\
					$(libdir)/libptscotch$(LIB)		\
					$(libdir)/libptscotcherrexit$(LIB)	\
					dgmap.h

dgord$(EXE)			:	dgord.c					\
					../libscotch/module.h			\
					../libscotch/common.h			\
					$(includedir)/ptscotch.h		\
					$(libdir)/libscotch$(LIB)		\
					$(libdir)/libptscotch$(LIB)		\
					$(libdir)/libptscotcherrexit$(LIB)	\
					dgord.h

dgpart$(EXE)			:	dgmap.c					\
					../libscotch/module.h			\
					../libscotch/common.h			\
					$(includedir)/ptscotch.h		\
					$(libdir)/libscotch$(LIB)		\
					$(libdir)/libptscotch$(LIB)		\
					$(libdir)/libptscotcherrexit$(LIB)	\
					dgmap.h
			 		$(CC) $(CFLAGS) -I$(includedir) -I../libscotch $(<) -DSCOTCH_COMPILE_PART -o $(@) -L$(libdir) -lptscotch -lscotch -lptscotcherrexit $(LDFLAGS)

dgscat$(EXE)			:	dgscat.c				\
					../libscotch/module.h			\
					../libscotch/common.h			\
					$(includedir)/ptscotch.h		\
					$(libdir)/libscotch$(LIB)		\
					$(libdir)/libptscotch$(LIB)		\
					$(libdir)/libptscotcherrexit$(LIB)	\
					dgscat.h

dgtst$(EXE)			:	dgtst.c					\
					../libscotch/module.h			\
					../libscotch/common.h			\
					$(includedir)/ptscotch.h		\
					$(libdir)/libscotch$(LIB)		\
					$(libdir)/libptscotch$(LIB)		\
					$(libdir)/libptscotcherrexit$(LIB)	\
					dgtst.h

gbase$(EXE)			:	gbase.c					\
					../libscotch/module.h			\
					../libscotch/common.h			\
					$(includedir)/scotch.h			\
					$(libdir)/libscotch$(LIB)		\
					$(libdir)/libscotcherrexit$(LIB)	\
					gbase.h

gcv$(EXE)			:	gcv.c					\
					../libscotch/module.h			\
					../libscotch/common.h			\
					$(includedir)/scotch.h			\
					$(libdir)/libscotch$(LIB)		\
					$(libdir)/libscotcherrexit$(LIB)	\
					gcv.h

gmap$(EXE)			:	gmap.c					\
					../libscotch/module.h			\
					../libscotch/common.h			\
					$(includedir)/scotch.h			\
					$(libdir)/libscotch$(LIB)		\
					$(libdir)/libscotcherrexit$(LIB)	\
					gmap.h

gmk_hy$(EXE)			:	gmk_hy.c				\
					../libscotch/module.h			\
					../libscotch/common.h			\
					$(includedir)/scotch.h			\
					$(libdir)/libscotch$(LIB)		\
					$(libdir)/libscotcherrexit$(LIB)	\
					gmk_hy.h

gmk_m2$(EXE)			:	gmk_m2.c				\
					../libscotch/module.h			\
					../libscotch/common.h			\
					$(includedir)/scotch.h			\
					$(libdir)/libscotch$(LIB)		\
					$(libdir)/libscotcherrexit$(LIB)	\
					gmk_m2.h

gmk_m3$(EXE)			:	gmk_m3.c				\
					../libscotch/module.h			\
					../libscotch/common.h			\
					$(includedir)/scotch.h			\
					$(libdir)/libscotch$(LIB)		\
					$(libdir)/libscotcherrexit$(LIB)	\
					gmk_m3.h

gmk_msh$(EXE)			:	gmk_msh.c				\
					../libscotch/module.h			\
					../libscotch/common.h			\
					$(includedir)/scotch.h			\
					$(libdir)/libscotch$(LIB)		\
					$(libdir)/libscotcherrexit$(LIB)	\
					gmk_msh.h

gmk_ub2$(EXE)			:	gmk_ub2.c				\
					../libscotch/module.h			\
					../libscotch/common.h			\
					$(includedir)/scotch.h			\
					$(libdir)/libscotch$(LIB)		\
					$(libdir)/libscotcherrexit$(LIB)	\
					gmk_ub2.h

gmtst$(EXE)			:	gmtst.c					\
					../libscotch/module.h			\
					../libscotch/common.h			\
					$(includedir)/scotch.h			\
					$(libdir)/libscotch$(LIB)		\
					$(libdir)/libscotcherrexit$(LIB)	\
					gmtst.h

gord$(EXE)			:	gord.c					\
					../libscotch/module.h			\
					../libscotch/common.h			\
					$(includedir)/scotch.h			\
					$(libdir)/libscotch$(LIB)		\
					$(libdir)/libscotcherrexit$(LIB)	\
					gord.h

gotst$(EXE)			:	gotst.c					\
					../libscotch/module.h			\
					../libscotch/common.h			\
					$(includedir)/scotch.h			\
					$(libdir)/libscotch$(LIB)		\
					$(libdir)/libscotcherrexit$(LIB)	\
					gotst.h

gout$(EXE)			:	gout_c.c				\
					gout_o.c				\
					$(includedir)/scotch.h			\
					$(libdir)/libscotch$(LIB)
					$(CC) $(CFLAGS) -I../libscotch -I$(includedir) gout_c.c gout_o.c -o $(@) -L$(libdir) -lscotch -lscotcherrexit $(LDFLAGS)

gpart$(EXE)			:	gmap.c					\
					../libscotch/module.h			\
					../libscotch/common.h			\
					$(includedir)/scotch.h			\
					$(libdir)/libscotch$(LIB)		\
					$(libdir)/libscotcherrexit$(LIB)	\
					gmap.h
			 		$(CC) $(CFLAGS) -I$(includedir) -I../libscotch $(<) -DSCOTCH_COMPILE_PART -o $(@) -L$(libdir) -lscotch -lscotcherrexit $(LDFLAGS)

gscat$(EXE)			:	gscat.c					\
					../libscotch/module.h			\
					../libscotch/common.h			\
					$(includedir)/scotch.h			\
					$(libdir)/libscotch$(LIB)		\
					$(libdir)/libscotcherrexit$(LIB)	\
					gscat.h

gtst$(EXE)			:	gtst.c					\
					../libscotch/module.h			\
					../libscotch/common.h			\
					$(includedir)/scotch.h			\
					$(libdir)/libscotch$(LIB)		\
					$(libdir)/libscotcherrexit$(LIB)	\
					gtst.h

mcv$(EXE)			:	mcv.c					\
					../libscotch/module.h			\
					../libscotch/common.h			\
					$(includedir)/scotch.h			\
					$(libdir)/libscotch$(LIB)		\
					$(libdir)/libscotcherrexit$(LIB)	\
					mcv.h

mmk_m2$(EXE)			:	mmk_m2.c				\
					../libscotch/module.h			\
					../libscotch/common.h			\
					$(includedir)/scotch.h			\
					$(libdir)/libscotch$(LIB)		\
					$(libdir)/libscotcherrexit$(LIB)	\
					mmk_m2.h

mmk_m3$(EXE)			:	mmk_m3.c				\
					../libscotch/module.h			\
					../libscotch/common.h			\
					$(includedir)/scotch.h			\
					$(libdir)/libscotch$(LIB)		\
					$(libdir)/libscotcherrexit$(LIB)	\
					mmk_m3.h

mord$(EXE)			:	mord.c					\
					../libscotch/module.h			\
					../libscotch/common.h			\
					$(includedir)/scotch.h			\
					$(libdir)/libscotch$(LIB)		\
					$(libdir)/libscotcherrexit$(LIB)	\
					mord.h

mtst$(EXE)			:	mtst.c					\
					../libscotch/module.h			\
					../libscotch/common.h			\
					$(includedir)/scotch.h			\
					$(libdir)/libscotch$(LIB)		\
					$(libdir)/libscotcherrexit$(LIB)	\
					mtst.h
