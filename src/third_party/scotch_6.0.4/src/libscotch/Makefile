## Copyright 2004,2007-2013 IPB, Universite de Bordeaux, INRIA & CNRS
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

includedir	= ../../include
libdir		= ../../lib

##
##  General inference rules.
##

include ../Makefile.inc

%$(OBJ)	:	%.c
				$(CC) $(CFLAGS) $(CLIBFLAGS) -c $(<) -o $(@)

%$(EXE)	:	%.c
		 		$(CC) $(CFLAGS) -DSCOTCH_VERSION=$(VERSION) -DSCOTCH_RELEASE=$(RELEASE) -DSCOTCH_PATCHLEVEL=$(PATCHLEVEL) $(<) -o $(@) $(LDFLAGS)

##
##  Project rules.
##

.PHONY				:	ptscotch	scotch	ptinstall	install	clean	realclean

scotch				:
					$(MAKE) CC="$(CCS)" CCD="$(CCS)"	\
					scotch.h				\
					scotchf.h				\
					libscotch$(LIB)				\
					libscotcherr$(LIB)			\
					libscotcherrexit$(LIB)

ptscotch			:	scotch
					$(MAKE) CFLAGS="$(CFLAGS) -DSCOTCH_PTSCOTCH" CC="$(CCP)"	\
					ptscotch.h							\
					ptscotchf.h							\
					libptscotch$(LIB)						\
					libptscotcherr$(LIB)						\
					libptscotcherrexit$(LIB)

install				:	$(includedir)/scotch.h			\
					$(includedir)/scotchf.h			\
					$(libdir)/libscotch$(LIB)

ptinstall			:	install					\
					$(includedir)/ptscotch.h		\
					$(includedir)/ptscotchf.h		\
					$(libdir)/libptscotch$(LIB)

clean				:
					-$(RM) *~ *$(OBJ) lib*$(LIB) parser_yy.c parser_ly.h parser_ll.c *scotch.h *scotchf.h y.output *dummysizes$(EXE)

realclean			:	clean


##
##  Installation dependencies
##

$(includedir)/scotch.h		:	scotch.h
					-$(CP) scotch.h $(includedir)

$(includedir)/scotchf.h		:	scotchf.h
					-$(CP) scotchf.h $(includedir)

$(includedir)/ptscotch.h	:	ptscotch.h
					-$(CP) ptscotch.h $(includedir)

$(includedir)/ptscotchf.h	:	ptscotchf.h
					-$(CP) ptscotchf.h $(includedir)

$(libdir)/libscotch$(LIB)	:	libscotch$(LIB)
					-$(CP) libscotch$(LIB) libscotcherr*$(LIB) $(libdir)

$(libdir)/libptscotch$(LIB)	:	libptscotch$(LIB)
					-$(CP) libptscotch$(LIB) libptscotcherr*$(LIB) $(libdir)

##
##  Library dependencies.
##

LIBPTSCOTCHDEPS =	bdgraph$(OBJ)				\
			bdgraph_bipart_bd$(OBJ)			\
			bdgraph_bipart_df$(OBJ)			\
			bdgraph_bipart_ex$(OBJ)			\
			bdgraph_bipart_ml$(OBJ)			\
			bdgraph_bipart_sq$(OBJ)			\
			bdgraph_bipart_st$(OBJ)			\
			bdgraph_bipart_zr$(OBJ)			\
			bdgraph_check$(OBJ)			\
			bdgraph_gather_all$(OBJ)		\
			bdgraph_store$(OBJ)			\
			comm$(OBJ)				\
			dgraph$(OBJ)				\
			dgraph_allreduce$(OBJ)			\
			dgraph_band$(OBJ)			\
			dgraph_build$(OBJ)			\
			dgraph_build_grid3d$(OBJ)		\
			dgraph_build_hcub$(OBJ)			\
			dgraph_check$(OBJ)			\
			dgraph_coarsen$(OBJ)			\
			dgraph_fold$(OBJ)			\
			dgraph_fold_comm$(OBJ)			\
			dgraph_fold_dup$(OBJ)			\
			dgraph_gather$(OBJ)			\
			dgraph_gather_all$(OBJ)			\
			dgraph_ghst$(OBJ)			\
			dgraph_halo$(OBJ)			\
			dgraph_induce$(OBJ)			\
			dgraph_io_load$(OBJ)			\
			dgraph_io_save$(OBJ)			\
			dgraph_match$(OBJ)			\
			dgraph_match_sync_coll$(OBJ)		\
			dgraph_match_sync_ptop$(OBJ)		\
			dgraph_match_check$(OBJ)		\
			dgraph_redist$(OBJ)			\
			dgraph_scatter$(OBJ)			\
			dgraph_view$(OBJ)			\
			dmapping$(OBJ)				\
			dmapping_io$(OBJ)			\
			dorder$(OBJ)				\
			dorder_gather$(OBJ)			\
			dorder_io$(OBJ)				\
			dorder_io_block$(OBJ)			\
			dorder_io_tree$(OBJ)			\
			dorder_perm$(OBJ)			\
			dorder_tree_dist$(OBJ)			\
			hdgraph$(OBJ)				\
			hdgraph_check$(OBJ)			\
			hdgraph_fold$(OBJ)			\
			hdgraph_gather$(OBJ)			\
			hdgraph_induce$(OBJ)			\
			hdgraph_order_nd$(OBJ)			\
			hdgraph_order_si$(OBJ)			\
			hdgraph_order_sq$(OBJ)			\
			hdgraph_order_st$(OBJ)			\
			kdgraph$(OBJ)				\
			kdgraph_gather$(OBJ)			\
			kdgraph_map_rb$(OBJ)			\
			kdgraph_map_rb_map$(OBJ)		\
			kdgraph_map_rb_part$(OBJ)		\
			kdgraph_map_st$(OBJ)			\
			library_dgraph$(OBJ)			\
			library_dgraph_f$(OBJ)			\
			library_dgraph_band$(OBJ)		\
			library_dgraph_band_f$(OBJ)		\
			library_dgraph_build$(OBJ)		\
			library_dgraph_build_f$(OBJ)		\
			library_dgraph_build_grid3d$(OBJ)	\
			library_dgraph_build_grid3d_f$(OBJ)	\
			library_dgraph_check$(OBJ)		\
			library_dgraph_check_f$(OBJ)		\
			library_dgraph_coarsen$(OBJ)		\
			library_dgraph_coarsen_f$(OBJ)		\
			library_dgraph_gather$(OBJ)		\
			library_dgraph_gather_f$(OBJ)		\
			library_dgraph_grow$(OBJ)		\
			library_dgraph_halo$(OBJ)		\
			library_dgraph_halo_f$(OBJ)		\
			library_dgraph_induce$(OBJ)		\
			library_dgraph_induce_f$(OBJ)		\
			library_dgraph_io_load$(OBJ)		\
			library_dgraph_io_load_f$(OBJ)		\
			library_dgraph_io_save$(OBJ)		\
			library_dgraph_io_save_f$(OBJ)		\
			library_dgraph_map$(OBJ)		\
			library_dgraph_map_f$(OBJ)		\
			library_dgraph_map_view$(OBJ)		\
			library_dgraph_map_view_f$(OBJ)		\
			library_dgraph_order$(OBJ)		\
			library_dgraph_order_f$(OBJ)		\
			library_dgraph_order_gather$(OBJ)	\
			library_dgraph_order_gather_f$(OBJ)	\
			library_dgraph_order_io$(OBJ)		\
			library_dgraph_order_io_f$(OBJ)		\
			library_dgraph_order_io_block$(OBJ)	\
			library_dgraph_order_io_block_f$(OBJ)	\
			library_dgraph_order_perm$(OBJ)		\
			library_dgraph_order_perm_f$(OBJ)	\
			library_dgraph_order_tree_dist$(OBJ)	\
			library_dgraph_order_tree_dist_f$(OBJ)	\
			library_dgraph_redist$(OBJ)		\
			library_dgraph_redist_f$(OBJ)		\
			library_dgraph_scatter$(OBJ)		\
			library_dgraph_scatter_f$(OBJ)		\
			library_dgraph_stat$(OBJ)		\
			library_dgraph_stat_f$(OBJ)		\
			library_dmapping$(OBJ)			\
			library_dorder$(OBJ)			\
			vdgraph$(OBJ)				\
			vdgraph_check$(OBJ)			\
			vdgraph_gather_all$(OBJ)		\
			vdgraph_separate_bd$(OBJ)		\
			vdgraph_separate_df$(OBJ)		\
			vdgraph_separate_ml$(OBJ)		\
			vdgraph_separate_sq$(OBJ)		\
			vdgraph_separate_st$(OBJ)		\
			vdgraph_separate_zr$(OBJ)		\
			vdgraph_store$(OBJ)

LIBSCOTCHDEPS	=	arch$(OBJ)				\
			arch_build$(OBJ)			\
			arch_cmplt$(OBJ)			\
			arch_cmpltw$(OBJ)			\
			arch_deco$(OBJ)				\
			arch_dist$(OBJ)				\
			arch_hcub$(OBJ)				\
			arch_mesh$(OBJ)				\
			arch_tleaf$(OBJ)			\
			arch_torus$(OBJ)			\
			arch_vcmplt$(OBJ)			\
			arch_vhcub$(OBJ)			\
			bgraph$(OBJ)				\
			bgraph_bipart_bd$(OBJ)			\
			bgraph_bipart_df$(OBJ)			\
			bgraph_bipart_ex$(OBJ)			\
			bgraph_bipart_fm$(OBJ)			\
			bgraph_bipart_gg$(OBJ)			\
			bgraph_bipart_gp$(OBJ)			\
			bgraph_bipart_ml$(OBJ)			\
			bgraph_bipart_st$(OBJ)			\
			bgraph_bipart_zr$(OBJ)			\
			bgraph_check$(OBJ)			\
			bgraph_store$(OBJ)			\
			common$(OBJ)				\
			common_file$(OBJ)			\
			common_file_compress$(OBJ)		\
			common_file_uncompress$(OBJ)		\
			common_integer$(OBJ)			\
			common_memory$(OBJ)			\
			common_string$(OBJ)			\
			common_stub$(OBJ)			\
			common_thread$(OBJ)			\
			fibo$(OBJ)				\
			gain$(OBJ)				\
			geom$(OBJ)				\
			graph$(OBJ)				\
			graph_base$(OBJ)			\
			graph_band$(OBJ)			\
			graph_check$(OBJ)			\
			graph_coarsen$(OBJ)			\
			graph_induce$(OBJ)			\
			graph_io$(OBJ)				\
			graph_io_chac$(OBJ)			\
			graph_io_habo$(OBJ)			\
			graph_io_mmkt$(OBJ)			\
			graph_io_scot$(OBJ)			\
			graph_list$(OBJ)			\
			graph_match$(OBJ)			\
			hall_order_hd$(OBJ)			\
			hall_order_hf$(OBJ)			\
			hall_order_hx$(OBJ)			\
			hgraph$(OBJ)				\
			hgraph_check$(OBJ)			\
			hgraph_induce$(OBJ)			\
			hgraph_order_bl$(OBJ)			\
			hgraph_order_cp$(OBJ)			\
			hgraph_order_gp$(OBJ)			\
			hgraph_order_hd$(OBJ)			\
			hgraph_order_hf$(OBJ)			\
			hgraph_order_hx$(OBJ)			\
			hgraph_order_kp$(OBJ)			\
			hgraph_order_nd$(OBJ)			\
			hgraph_order_si$(OBJ)			\
			hgraph_order_st$(OBJ)			\
			hmesh$(OBJ)				\
			hmesh_check$(OBJ)			\
			hmesh_hgraph$(OBJ)			\
			hmesh_induce$(OBJ)			\
			hmesh_mesh$(OBJ)			\
			hmesh_order_bl$(OBJ)			\
			hmesh_order_cp$(OBJ)			\
			hmesh_order_gr$(OBJ)			\
			hmesh_order_gp$(OBJ)			\
			hmesh_order_hd$(OBJ)			\
			hmesh_order_hf$(OBJ)			\
			hmesh_order_hx$(OBJ)			\
			hmesh_order_nd$(OBJ)			\
			hmesh_order_si$(OBJ)			\
			hmesh_order_st$(OBJ)			\
			kgraph$(OBJ)				\
			kgraph_band$(OBJ)			\
			kgraph_check$(OBJ)			\
			kgraph_map_bd$(OBJ)			\
			kgraph_map_cp$(OBJ)			\
			kgraph_map_df$(OBJ)			\
			kgraph_map_ex$(OBJ)			\
			kgraph_map_fm$(OBJ)			\
			kgraph_map_ml$(OBJ)			\
			kgraph_map_rb$(OBJ)			\
			kgraph_map_rb_map$(OBJ)			\
			kgraph_map_rb_part$(OBJ)		\
			kgraph_map_st$(OBJ)			\
			kgraph_store$(OBJ)			\
			library_arch$(OBJ)			\
			library_arch_f$(OBJ)			\
			library_arch_build$(OBJ)		\
			library_arch_build_f$(OBJ)		\
			library_common_f$(OBJ)			\
			library_geom$(OBJ)			\
			library_geom_f$(OBJ)			\
			library_graph$(OBJ)			\
			library_graph_f$(OBJ)			\
			library_graph_base$(OBJ)		\
			library_graph_base_f$(OBJ)		\
			library_graph_check$(OBJ)		\
			library_graph_check_f$(OBJ)		\
			library_graph_coarsen$(OBJ)		\
			library_graph_coarsen_f$(OBJ)		\
			library_graph_color$(OBJ)		\
			library_graph_color_f$(OBJ)		\
			library_graph_io_chac$(OBJ)		\
			library_graph_io_chac_f$(OBJ)		\
			library_graph_io_habo$(OBJ)		\
			library_graph_io_habo_f$(OBJ)		\
			library_graph_io_mmkt$(OBJ)		\
			library_graph_io_mmkt_f$(OBJ)		\
			library_graph_io_scot$(OBJ)		\
			library_graph_io_scot_f$(OBJ)		\
			library_graph_map$(OBJ)			\
			library_graph_map_f$(OBJ)		\
			library_graph_map_io$(OBJ)		\
			library_graph_map_io_f$(OBJ)		\
			library_graph_map_view$(OBJ)		\
			library_graph_map_view_f$(OBJ)		\
			library_graph_order$(OBJ)		\
			library_graph_order_f$(OBJ)		\
			library_graph_part_ovl$(OBJ)		\
			library_graph_part_ovl_f$(OBJ)		\
			library_mapping$(OBJ)			\
			library_memory$(OBJ)			\
			library_memory_f$(OBJ)			\
			library_mesh$(OBJ)			\
			library_mesh_f$(OBJ)			\
			library_mesh_graph$(OBJ)		\
			library_mesh_graph_f$(OBJ)		\
			library_mesh_io_habo$(OBJ)		\
			library_mesh_io_habo_f$(OBJ)		\
			library_mesh_io_scot$(OBJ)		\
			library_mesh_io_scot_f$(OBJ)		\
			library_mesh_order$(OBJ)		\
			library_mesh_order_f$(OBJ)		\
			library_order$(OBJ)			\
			library_parser$(OBJ)			\
			library_parser_f$(OBJ)			\
			library_random$(OBJ)			\
			library_random_f$(OBJ)			\
			library_strat$(OBJ)			\
			library_version$(OBJ)			\
			library_version_f$(OBJ)			\
			mapping$(OBJ)				\
			mapping_io$(OBJ)			\
			mesh$(OBJ)				\
			mesh_check$(OBJ)			\
			mesh_coarsen$(OBJ)			\
			mesh_graph$(OBJ)			\
			mesh_induce_sepa$(OBJ)			\
			mesh_io$(OBJ)				\
			mesh_io_habo$(OBJ)			\
			mesh_io_scot$(OBJ)			\
			order$(OBJ)				\
			order_check$(OBJ)			\
			order_io$(OBJ)				\
			parser$(OBJ)				\
			parser_ll$(OBJ)				\
			parser_yy$(OBJ)				\
			vgraph$(OBJ)				\
			vgraph_check$(OBJ)			\
			vgraph_separate_bd$(OBJ)		\
			vgraph_separate_df$(OBJ)		\
			vgraph_separate_es$(OBJ)		\
			vgraph_separate_fm$(OBJ)		\
			vgraph_separate_gg$(OBJ)		\
			vgraph_separate_gp$(OBJ)		\
			vgraph_separate_ml$(OBJ)		\
			vgraph_separate_st$(OBJ)		\
			vgraph_separate_th$(OBJ)		\
			vgraph_separate_vw$(OBJ)		\
			vgraph_separate_zr$(OBJ)		\
			vgraph_store$(OBJ)			\
			vmesh$(OBJ)				\
			vmesh_check$(OBJ)			\
			vmesh_separate_fm$(OBJ)			\
			vmesh_separate_gg$(OBJ)			\
			vmesh_separate_gr$(OBJ)			\
			vmesh_separate_ml$(OBJ)			\
			vmesh_separate_zr$(OBJ)			\
			vmesh_separate_st$(OBJ)			\
			vmesh_store$(OBJ)			\
			wgraph$(OBJ)				\
			wgraph_check$(OBJ)			\
			wgraph_part_fm$(OBJ)  			\
			wgraph_part_gg$(OBJ)          		\
			wgraph_part_gp$(OBJ)	        	\
			wgraph_part_ml$(OBJ)   			\
			wgraph_part_rb$(OBJ)  			\
			wgraph_part_st$(OBJ)  			\
			wgraph_part_zr$(OBJ)  			\
			wgraph_store$(OBJ)

##
##  Todo list.
##

common$(OBJ)			:	common.c				\
					module.h				\
					common.h
					$(CC) $(CFLAGS) $(CLIBFLAGS) -c $(<) -DSCOTCH_COMMON_RENAME -o $(@)

common_file$(OBJ)		:	common_file.c				\
					module.h				\
					common.h				\
					common_file.h
					$(CC) $(CFLAGS) $(CLIBFLAGS) -c $(<) -DSCOTCH_COMMON_RENAME -o $(@)

common_file_compress$(OBJ)	:	common_file_compress.c			\
					module.h				\
					common.h				\
					common_file.h				\
					common_file_compress.h
					$(CC) $(CFLAGS) $(CLIBFLAGS) -c $(<) -DSCOTCH_COMMON_RENAME -o $(@)

common_file_uncompress$(OBJ)	:	common_file_uncompress.c		\
					module.h				\
					common.h				\
					common_file.h				\
					common_file_compress.h
					$(CC) $(CFLAGS) $(CLIBFLAGS) -c $(<) -DSCOTCH_COMMON_RENAME -o $(@)

common_integer$(OBJ)		:	common_integer.c			\
					module.h				\
					common.h				\
					common_sort.c
					$(CC) $(CFLAGS) $(CLIBFLAGS) -c $(<) -DSCOTCH_COMMON_RENAME -o $(@)

common_memory$(OBJ)		:	common_memory.c				\
					module.h				\
					common.h
					$(CC) $(CFLAGS) $(CLIBFLAGS) -c $(<) -DSCOTCH_COMMON_RENAME -o $(@)

common_string$(OBJ)		:	common_string.c				\
					module.h				\
					common.h
					$(CC) $(CFLAGS) $(CLIBFLAGS) -c $(<) -DSCOTCH_COMMON_RENAME -o $(@)

common_stub$(OBJ)		:	common_stub.c				\
					module.h				\
					common.h
					$(CC) $(CFLAGS) $(CLIBFLAGS) -c $(<) -DSCOTCH_COMMON_RENAME -o $(@)

common_thread$(OBJ)		:	common_thread.c				\
					module.h				\
					common.h
					$(CC) $(CFLAGS) $(CLIBFLAGS) -c $(<) -DSCOTCH_COMMON_RENAME -o $(@)

arch$(OBJ)			:	arch.c					\
					module.h				\
					common.h				\
					arch.h					\
					arch_cmplt.h				\
					arch_cmpltw.h				\
					arch_deco.h				\
					arch_dist.h				\
					arch_hcub.h				\
					arch_mesh.h				\
					arch_tleaf.h				\
					arch_torus.h				\
					arch_vcmplt.h				\
					arch_vhcub.h

arch_build$(OBJ)		:	arch_build.c				\
					module.h				\
					common.h				\
					parser.h				\
					graph.h					\
					arch$(OBJ)				\
					arch_deco.h				\
					arch_vcmplt.h				\
					mapping.h				\
					bgraph.h				\
					bgraph_bipart_st.h			\
					arch_build.h

arch_cmplt$(OBJ)		:	arch_cmplt.c				\
					module.h				\
					common.h				\
					arch.h					\
					arch_cmplt.h

arch_cmpltw$(OBJ)		:	arch_cmpltw.c				\
					module.h				\
					common.h				\
					arch.h					\
					arch_cmpltw.h

arch_deco$(OBJ)			:	arch_deco.c				\
					module.h				\
					common.h				\
					arch.h					\
					arch_deco.h

arch_dist$(OBJ)			:	arch_dist.c				\
					module.h				\
					common.h				\
					arch.h					\
					arch_dist.h

arch_hcub$(OBJ)			:	arch_hcub.c				\
					module.h				\
					common.h				\
					arch.h					\
					arch_hcub.h

arch_mesh$(OBJ)			:	arch_mesh.c				\
					module.h				\
					common.h				\
					arch.h					\
					arch_mesh.h

arch_tleaf$(OBJ)		:	arch_tleaf.c				\
					module.h				\
					common.h				\
					arch.h					\
					arch_tleaf.h

arch_torus$(OBJ)		:	arch_torus.c				\
					module.h				\
					common.h				\
					arch.h					\
					arch_torus.h

arch_vcpmlt$(OBJ)		:	arch_vcmplt.c				\
					module.h				\
					common.h				\
					arch.h					\
					arch_vcmplt.h

arch_vhcub$(OBJ)		:	arch_vhcub.c				\
					module.h				\
					common.h				\
					arch.h					\
					arch_vhcub.h

bdgraph$(OBJ)			:	bdgraph.c				\
					module.h				\
					common.h				\
					arch$(OBJ)				\
					dgraph.h				\
					dmapping.h				\
					bdgraph.h

bdgraph_bipart_bd$(OBJ)		:	bdgraph_bipart_bd.c			\
					module.h				\
					common.h				\
					parser.h				\
					arch$(OBJ)				\
					dgraph.h				\
					dgraph_halo.h				\
					bdgraph.h				\
					bdgraph_bipart_bd.h			\
					bdgraph_bipart_st.h

bdgraph_bipart_df$(OBJ)		:	bdgraph_bipart_df.c			\
					module.h				\
					common.h				\
					arch$(OBJ)				\
					dgraph.h				\
					bdgraph.h				\
					bdgraph_bipart_df.h

bdgraph_bipart_ex$(OBJ)		:	bdgraph_bipart_ex.c			\
					module.h				\
					common.h				\
					arch$(OBJ)				\
					dgraph.h				\
					bdgraph.h				\
					bdgraph_bipart_ex.h

bdgraph_bipart_ml$(OBJ)		:	bdgraph_bipart_ml.c			\
					module.h				\
					common.h				\
					parser.h				\
					arch$(OBJ)				\
					dgraph.h				\
					dgraph_coarsen.h			\
					bdgraph.h				\
					bdgraph_bipart_ml.h			\
					bdgraph_bipart_st.h

bdgraph_bipart_sq$(OBJ)		:	bdgraph_bipart_sq.c			\
					module.h				\
					common.h				\
					comm.h					\
					parser.h				\
					arch$(OBJ)				\
					graph.h					\
					bgraph.h				\
					bgraph_bipart_st.h			\
					dgraph.h				\
					bdgraph.h				\
					bdgraph_bipart_sq.h

bdgraph_bipart_st$(OBJ)		:	bdgraph_bipart_st.c			\
					module.h				\
					common.h				\
					gain.h					\
					parser.h				\
					graph.h					\
					arch$(OBJ)				\
					bgraph.h				\
					bgraph_bipart_st.h			\
					dgraph.h				\
					dgraph_coarsen.h			\
					bdgraph.h				\
					bdgraph_bipart_bd.h			\
					bdgraph_bipart_df.h			\
					bdgraph_bipart_ex.h			\
					bdgraph_bipart_ml.h			\
					bdgraph_bipart_sq.h			\
					bdgraph_bipart_st.h			\
					bdgraph_bipart_zr.h

bdgraph_bipart_zr$(OBJ)		:	bdgraph_bipart_zr.c			\
					module.h				\
					common.h				\
					arch$(OBJ)				\
					dgraph.h				\
					bdgraph.h				\
					bdgraph_bipart_zr.h

bdgraph_check$(OBJ)		:	bdgraph_check.c				\
					module.h				\
					common.h				\
					arch$(OBJ)				\
					dgraph.h				\
					bdgraph.h

bdgraph_gather_all$(OBJ)	:	bdgraph_gather_all.c			\
					module.h				\
					common.h				\
					comm.h					\
					arch$(OBJ)				\
					graph.h					\
					bgraph.h				\
					dgraph.h				\
					bdgraph.h

bdgraph_store$(OBJ)		:	bdgraph_store.c				\
					module.h				\
					common.h				\
					arch$(OBJ)				\
					dgraph.h				\
					bdgraph.h

bgraph$(OBJ)			:	bgraph.c				\
					module.h				\
					common.h				\
					graph.h					\
					arch$(OBJ)				\
					mapping.h				\
					bgraph.h

bgraph_bipart_bd$(OBJ)		:	bgraph_bipart_bd.c			\
					module.h				\
					common.h				\
					parser.h				\
					graph.h					\
					arch$(OBJ)				\
					bgraph.h				\
					bgraph_bipart_bd.h			\
					bgraph_bipart_st.h

bgraph_bipart_df$(OBJ)		:	bgraph_bipart_df.c			\
					bgraph_bipart_df_loop.c			\
					module.h				\
					common.h				\
					graph.h					\
					arch$(OBJ)				\
					bgraph.h				\
					bgraph_bipart_df.h

bgraph_bipart_ex$(OBJ)		:	bgraph_bipart_ex.c			\
					module.h				\
					common.h				\
					gain.h					\
					graph.h					\
					arch$(OBJ)				\
					bgraph.h				\
					bgraph_bipart_ex.h			\
					bgraph_bipart_fm.h			\
					bgraph_bipart_gg.h

bgraph_bipart_fm$(OBJ)		:	bgraph_bipart_fm.c			\
					gain.h					\
					fibo.h					\
					module.h				\
					common.h				\
					graph.h					\
					arch$(OBJ)				\
					bgraph.h				\
					bgraph_bipart_fm.h

bgraph_bipart_gg$(OBJ)		:	bgraph_bipart_gg.c			\
					gain.h					\
					fibo.h					\
					module.h				\
					common.h				\
					graph.h					\
					arch$(OBJ)				\
					bgraph.h				\
					bgraph_bipart_gg.h

bgraph_bipart_gp$(OBJ)		:	bgraph_bipart_gp.c			\
					module.h				\
					common.h				\
					graph.h					\
					arch$(OBJ)				\
					bgraph.h				\
					bgraph_bipart_gp.h

bgraph_bipart_ml$(OBJ)		:	bgraph_bipart_ml.c			\
					module.h				\
					common.h				\
					parser.h				\
					graph.h					\
					graph_coarsen.h				\
					arch$(OBJ)				\
					bgraph.h				\
					bgraph_bipart_ml.h			\
					bgraph_bipart_st.h

bgraph_bipart_st$(OBJ)		:	bgraph_bipart_st.c			\
					gain.h					\
					fibo.h					\
					module.h				\
					common.h				\
					parser.h				\
					graph.h					\
					arch$(OBJ)				\
					bgraph.h				\
					bgraph_bipart_bd.h			\
					bgraph_bipart_df.h			\
					bgraph_bipart_ex.h			\
					bgraph_bipart_fm.h			\
					bgraph_bipart_gg.h			\
					bgraph_bipart_gp.h			\
					bgraph_bipart_ml.h			\
					bgraph_bipart_st.h			\
					bgraph_bipart_zr.h

bgraph_bipart_zr$(OBJ)		:	bgraph_bipart_zr.c			\
					module.h				\
					common.h				\
					graph.h					\
					arch$(OBJ)				\
					bgraph.h				\
					bgraph_bipart_zr.h

bgraph_check$(OBJ)		:	bgraph_check.c				\
					module.h				\
					common.h				\
					graph.h					\
					arch$(OBJ)				\
					bgraph.h

bgraph_store$(OBJ)		:	bgraph_store.c				\
					module.h				\
					common.h				\
					graph.h					\
					arch$(OBJ)				\
					bgraph.h

comm$(OBJ)			:	comm.c					\
					module.h				\
					common.h				\
					comm.h

dgraph$(OBJ)			:	dgraph.c				\
					module.h				\
					common.h				\
					dgraph.h

dgraph_allreduce$(OBJ)		:	dgraph_allreduce.c			\
					module.h				\
					common.h				\
					dgraph.h				\
					dgraph_allreduce.h

dgraph_band$(OBJ)		:	dgraph_band.c				\
					module.h				\
					common.h				\
					dgraph.h

dgraph_build$(OBJ)		:	dgraph_build.c				\
					module.h				\
					common.h				\
					dgraph.h				\
					dgraph_allreduce.h			\
					dgraph_build.h

dgraph_build_grid3d$(OBJ)	:	dgraph_build_grid3d.c			\
					module.h				\
					common.h				\
					dgraph.h				\
					dgraph_build_grid3d.h

dgraph_build_hcub$(OBJ)		:	dgraph_build_hcub.c			\
					module.h				\
					common.h				\
					dgraph.h

dgraph_check$(OBJ)		:	dgraph_check.c				\
					module.h				\
					common.h				\
					dgraph.h

dgraph_coarsen$(OBJ)		:	dgraph_coarsen.c			\
					module.h				\
					common.h				\
					dgraph.h				\
					dgraph_allreduce.h			\
					dgraph_coarsen.h			\
					dgraph_match.h

dgraph_fold$(OBJ)		:	dgraph_fold.c				\
					module.h				\
					common.h				\
					comm.h					\
					dgraph.h				\
					dgraph_fold.h				\
					dgraph_fold_comm.h

dgraph_fold_comm$(OBJ)		:	dgraph_fold_comm.c			\
					module.h				\
					common.h				\
					comm.h					\
					dgraph.h				\
					dgraph_fold_comm.h

dgraph_fold_dup$(OBJ)		:	dgraph_fold_dup.c			\
					module.h				\
					common.h				\
					dgraph.h				\
					dgraph_fold_dup.h

dgraph_gather$(OBJ)		:	dgraph_gather.c				\
					module.h				\
					common.h				\
					graph.h					\
					dgraph.h

dgraph_gather_all$(OBJ)		:	dgraph_gather_all.c			\
					module.h				\
					common.h				\
					comm.h					\
					graph.h					\
					dgraph.h

dgraph_ghst$(OBJ)		:	dgraph_ghst.c				\
					module.h				\
					common.h				\
					dgraph.h				\
					dgraph_ghst.h

dgraph_halo$(OBJ)		:	dgraph_halo.c				\
					dgraph_halo_fill.c			\
					module.h				\
					common.h				\
					dgraph.h				\
					dgraph_halo.h

dgraph_induce$(OBJ)		:	dgraph_induce.c				\
					module.h				\
					common.h				\
					dgraph.h

dgraph_io_load$(OBJ)		:	dgraph_io_load.c			\
					module.h				\
					common.h				\
					dgraph.h				\
					dgraph_io_load.h

dgraph_io_save$(OBJ)		:	dgraph_io_save.c			\
					module.h				\
					common.h				\
					dgraph.h

dgraph_match$(OBJ)		:	dgraph_match.c				\
					dgraph_match_scan.c			\
					module.h				\
					common.h				\
					dgraph.h				\
					dgraph_coarsen.h			\
					dgraph_match.h

dgraph_match_check$(OBJ)	:	dgraph_match_check.c			\
					module.h				\
					common.h				\
					dgraph.h				\
					dgraph_coarsen.h			\
					dgraph_match.h

dgraph_match_sync_coll$(OBJ)	:	dgraph_match_sync_coll.c		\
					module.h				\
					common.h				\
					dgraph.h				\
					dgraph_coarsen.h                        \
					dgraph_match.h				\
					graph.h

dgraph_match_sync_ptop$(OBJ)	:	dgraph_match_sync_ptop.c		\
					module.h				\
					common.h				\
					dgraph.h				\
					dgraph_coarsen.h                        \
					dgraph_match.h				\
					graph.h

dgraph_redist$(OBJ)		:	dgraph_redist.c				\
					module.h				\
					common.h				\
					graph.h					\
					dgraph.h

dgraph_scatter$(OBJ)		:	dgraph_scatter.c			\
					module.h				\
					common.h				\
					graph.h					\
					dgraph.h

dgraph_view$(OBJ)		:	dgraph_view.c				\
					module.h				\
					common.h				\
					dgraph.h

dmapping$(OBJ)			:	dmapping.c				\
					module.h				\
					common.h				\
					arch$(OBJ)				\
					dgraph.h				\
					dmapping.h

dmapping_io$(OBJ)		:	dmapping_io.c				\
					module.h				\
					common.h				\
					comm.h					\
					arch$(OBJ)				\
					dgraph.h				\
					dgraph_allreduce.h			\
					dmapping.h

dorder$(OBJ)			:	dorder.c				\
					module.h				\
					common.h				\
					dgraph.h				\
					dorder.h

dorder_gather$(OBJ)		:	dorder_gather.c				\
					module.h				\
					common.h				\
					dgraph.h				\
					dgraph_allreduce.h			\
					dorder.h				\
					dorder_gather.h				\
					order.h

dorder_io$(OBJ)			:	dorder_io.c				\
					module.h				\
					common.h				\
					comm.h					\
					dgraph.h				\
					dorder.h				\
					order.h

dorder_io_block$(OBJ)		:	dorder_io_block.c			\
					module.h				\
					common.h				\
					dgraph.h				\
					dorder.h				\
					order.h

dorder_io_tree$(OBJ)		:	dorder_io_tree.c			\
					module.h				\
					common.h				\
					comm.h					\
					dgraph.h				\
					dorder.h				\
					order.h

dorder_perm$(OBJ)		:	dorder_perm.c				\
					module.h				\
					common.h				\
					dgraph.h				\
					dorder.h				\
					dorder_perm.h

dorder_tree_dist$(OBJ)		:	dorder_tree_dist.c			\
					module.h				\
					common.h				\
					dgraph.h				\
					dorder.h

fibo$(OBJ)			:	fibo.c					\
					module.h				\
					common.h				\
					fibo.h                                  

gain$(OBJ)			:	gain.c					\
					module.h				\
					common.h				\
					gain.h

geom$(OBJ)			:	geom.c					\
					module.h				\
					common.h				\
					geom.h

graph$(OBJ)			:	graph.c					\
					module.h				\
					common.h				\
					graph.h

graph_band$(OBJ)		:	graph_band.c				\
					module.h				\
					common.h				\
					graph.h

graph_base$(OBJ)		:	graph_base.c				\
					module.h				\
					common.h				\
					graph.h

graph_check$(OBJ)		:	graph_check.c				\
					module.h				\
					common.h				\
					graph.h

graph_coarsen$(OBJ)		:	graph_coarsen.c				\
					graph_coarsen_edge.c			\
					module.h				\
					common.h				\
					graph.h					\
					graph_coarsen.h

graph_induce$(OBJ)		:	graph_induce.c				\
					module.h				\
					common.h				\
					graph.h					\
					graph_induce.h

graph_io$(OBJ)			:	graph_io.c				\
					module.h				\
					common.h				\
					graph.h					\
					graph_io.h

graph_io_chac$(OBJ)		:	graph_io_chac.c				\
					module.h				\
					common.h				\
					geom.h					\
					graph.h

graph_io_habo$(OBJ)		:	graph_io_habo.c				\
					module.h				\
					common.h				\
					geom.h					\
					graph.h					\
					graph_io_habo.h

graph_io_mmkt$(OBJ)		:	graph_io_mmkt.c				\
					module.h				\
					common.h				\
					geom.h					\
					graph.h					\
					graph_io_mmkt.h

graph_io_scot$(OBJ)		:	graph_io_scot.c				\
					module.h				\
					common.h				\
					geom.h					\
					graph.h					\
					graph_io_scot.h

graph_list$(OBJ)		:	graph_list.c				\
					module.h				\
					common.h				\
					graph.h

graph_match$(OBJ)		:	graph_match.c				\
					graph_match_scan.c			\
					module.h				\
					common.h				\
					graph.h					\
					graph_coarsen.h				\
					graph_match.h

hall_order_hd$(OBJ)		:	hall_order_hd.c				\
					module.h				\
					common.h				\
					graph.h					\
					hall_order_hd.h

hall_order_hf$(OBJ)		:	hall_order_hf.c				\
					module.h				\
					common.h				\
					graph.h					\
					hall_order_hf.h

hall_order_hx$(OBJ)		:	hall_order_hx.c				\
					module.h				\
					common.h				\
					graph.h					\
					order.h					\
					hall_order_hx.h

hdgraph$(OBJ)			:	hdgraph.c				\
					module.h				\
					common.h				\
					dgraph.h				\
					hdgraph.h

hdgraph_check$(OBJ)		:	hdgraph_check.c				\
					module.h				\
					common.h				\
					dgraph.h				\
					hdgraph.h

hdgraph_fold$(OBJ)		:	hdgraph_fold.c				\
					module.h				\
					common.h				\
					comm.h					\
					dgraph.h				\
					dgraph_fold_comm.h			\
					hdgraph.h				\
					hdgraph_fold.h

hdgraph_gather$(OBJ)		:	hdgraph_gather.c			\
					module.h				\
					common.h				\
					comm.h					\
					graph.h					\
					hgraph.h				\
					dgraph.h				\
					hdgraph.h

hdgraph_induce$(OBJ)		:	hdgraph_induce.c			\
					module.h				\
					common.h				\
					dgraph.h				\
					hdgraph.h

hdgraph_order_nd$(OBJ)		:	hdgraph_order_nd.c			\
					module.h				\
					common.h				\
					parser.h				\
					graph.h					\
					order.h					\
					hgraph.h				\
					hgraph_order_st.h			\
					dgraph.h				\
					dorder.h				\
					hdgraph.h				\
					hdgraph_order_nd.h			\
					hdgraph_order_sq.h			\
					hdgraph_order_st.h			\
					vdgraph.h				\
					vdgraph_separate_st.h

hdgraph_order_si$(OBJ)		:	hdgraph_order_si.c			\
					module.h				\
					common.h				\
					dgraph.h				\
					dorder.h				\
					hdgraph.h				\
					hdgraph_order_si.h

hdgraph_order_sq$(OBJ)		:	hdgraph_order_sq.c			\
					module.h				\
					common.h				\
					parser.h				\
					graph.h					\
					order.h					\
					hgraph.h				\
					hgraph_order_st.h			\
					dgraph.h				\
					dorder.h				\
					hdgraph.h				\
					hdgraph_order_sq.h

hdgraph_order_st$(OBJ)		:	hdgraph_order_st.c			\
					module.h				\
					common.h				\
					parser.h				\
					dgraph.h				\
					dgraph_coarsen.h			\
					dorder.h				\
					hdgraph.h				\
					vdgraph.h				\
					vdgraph_separate_st.h			\
					hdgraph_order_nd.h			\
					hdgraph_order_si.h			\
					hdgraph_order_sq.h			\
					hdgraph_order_st.h

hgraph$(OBJ)			:	hgraph.c				\
					module.h				\
					common.h				\
					graph.h					\
					hgraph.h

hgraph_check$(OBJ)		:	hgraph_check.c				\
					module.h				\
					common.h				\
					graph.h					\
					hgraph.h

hgraph_induce$(OBJ)		:	hgraph_induce.c				\
					hgraph_induce_edge.c			\
					module.h				\
					common.h				\
					graph.h					\
					hgraph.h				\
					hgraph_induce.h

hgraph_order_bl$(OBJ)		:	hgraph_order_bl.c			\
					module.h				\
					common.h				\
					parser.h				\
					graph.h					\
					order.h					\
					hgraph.h				\
					hgraph_order_bl.h

hgraph_order_cp$(OBJ)		:	hgraph_order_cp.c			\
					module.h				\
					common.h				\
					parser.h				\
					graph.h					\
					order.h					\
					hgraph.h				\
					hgraph_order_cp.h			\
					hgraph_order_st.h

hgraph_order_gp$(OBJ)		:	hgraph_order_gp.c			\
					module.h				\
					common.h				\
					graph.h					\
					order.h					\
					hgraph.h				\
					hgraph_order_gp.h

hgraph_order_hd$(OBJ)		:	hgraph_order_hd.c			\
					module.h				\
					common.h				\
					graph.h					\
					order.h					\
					hgraph.h				\
					hall_order_hd.h				\
					hall_order_hx.h				\
					hgraph_order_hd.h			\
					hgraph_order_hx.h			\
					hgraph_order_si.h

hgraph_order_hf$(OBJ)		:	hgraph_order_hf.c			\
					module.h				\
					common.h				\
					graph.h					\
					order.h					\
					hgraph.h				\
					hall_order_hf.h				\
					hall_order_hx.h				\
					hgraph_order_hf.h			\
					hgraph_order_hx.h			\
					hgraph_order_si.h

hgraph_order_hx$(OBJ)		:	hgraph_order_hx.c			\
					module.h				\
					common.h				\
					graph.h					\
					hgraph.h				\
					hgraph_order_hx.h

hgraph_order_kp$(OBJ)		:	hgraph_order_kp.c			\
					module.h				\
					common.h				\
					parser.h				\
					graph.h					\
					arch.h					\
					mapping.h				\
					order.h					\
					hgraph.h				\
					hgraph_order_kp.h			\
					hgraph_order_si.h			\
					hgraph_order_st.h			\
					kgraph.h				\
					kgraph_map_st.h				\
					scotch.h

hgraph_order_nd$(OBJ)		:	hgraph_order_nd.c			\
					module.h				\
					common.h				\
					parser.h				\
					graph.h					\
					order.h					\
					hgraph.h				\
					hgraph_order_nd.h			\
					hgraph_order_st.h			\
					vgraph.h				\
					vgraph_separate_st.h

hgraph_order_si$(OBJ)		:	hgraph_order_si.c			\
					module.h				\
					common.h				\
					graph.h					\
					order.h					\
					hgraph.h				\
					hgraph_order_si.h

hgraph_order_st$(OBJ)		:	hgraph_order_st.c			\
					module.h				\
					common.h				\
					parser.h				\
					graph.h					\
					arch$(OBJ)				\
					mapping.h				\
					order.h					\
					hgraph.h				\
					hgraph_order_bl.h			\
					hgraph_order_cp.h			\
					hgraph_order_gp.h			\
					hgraph_order_hd.h			\
					hgraph_order_hf.h			\
					hgraph_order_kp.h			\
					hgraph_order_nd.h			\
					hgraph_order_si.h			\
					hgraph_order_st.h			\
					kgraph.h				\
					kgraph_map_st.h				\
					vgraph.h				\
					vgraph_separate_st.h

hmesh$(OBJ)			:	hmesh.c					\
					module.h				\
					common.h				\
					graph.h					\
					mesh.h					\
					hmesh.h

hmesh_check$(OBJ)		:	hmesh_check.c				\
					module.h				\
					common.h				\
					graph.h					\
					mesh.h					\
					hmesh.h

hmesh_hgraph$(OBJ)		:	hmesh_hgraph.c				\
					module.h				\
					common.h				\
					graph.h					\
					hgraph.h				\
					mesh.h					\
					hmesh.h					\
					hmesh_hgraph.h

hmesh_induce$(OBJ)		:	hmesh_induce.c				\
					module.h				\
					common.h				\
					graph.h					\
					mesh.h					\
					hmesh.h

hmesh_mesh$(OBJ)		:	hmesh_mesh.c				\
					module.h				\
					common.h				\
					graph.h					\
					mesh.h					\
					hmesh.h

hmesh_order_bl$(OBJ)		:	hmesh_order_bl.c			\
					module.h				\
					common.h				\
					parser.h				\
					graph.h					\
					order.h					\
					mesh.h					\
					hmesh.h					\
					hmesh_order_bl.h			\
					hmesh_order_st.h

hmesh_order_cp$(OBJ)		:	hmesh_order_cp.c			\
					module.h				\
					common.h				\
					parser.h				\
					graph.h					\
					order.h					\
					mesh.h					\
					hmesh.h					\
					hmesh_order_cp.h			\
					hmesh_order_st.h

hmesh_order_gp$(OBJ)		:	hmesh_order_gp.c			\
					module.h				\
					common.h				\
					graph.h					\
					order.h					\
					mesh.h					\
					hmesh.h					\
					hmesh_order_gp.h

hmesh_order_gp$(OBJ)		:	hmesh_order_gp.c			\
					module.h				\
					common.h				\
					graph.h					\
					order.h					\
					mesh.h					\
					hmesh.h					\
					hmesh_order_gp.h

hmesh_order_hd$(OBJ)		:	hmesh_order_hd.c			\
					module.h				\
					common.h				\
					graph.h					\
					order.h					\
					mesh.h					\
					hmesh.h					\
					hall_order_hd.h				\
					hall_order_hx.h				\
					hmesh_order_hd.h			\
					hmesh_order_si.h

hmesh_order_hf$(OBJ)		:	hmesh_order_hf.c			\
					module.h				\
					common.h				\
					graph.h					\
					order.h					\
					mesh.h					\
					hmesh.h					\
					hall_order_hf.h				\
					hall_order_hx.h				\
					hmesh_order_hf.h			\
					hmesh_order_si.h

hmesh_order_hx$(OBJ)		:	hmesh_order_hx.c			\
					module.h				\
					common.h				\
					graph.h					\
					mesh.h					\
					hmesh.h					\
					hmesh_order_hx.h

hmesh_order_nd$(OBJ)		:	hmesh_order_nd.c			\
					module.h				\
					common.h				\
					parser.h				\
					graph.h					\
					order.h					\
					mesh.h					\
					hmesh.h					\
					hmesh_order_nd.h			\
					hmesh_order_st.h			\
					vmesh.h					\
					vmesh_separate_st.h

hmesh_order_si$(OBJ)		:	hmesh_order_si.c			\
					module.h				\
					common.h				\
					graph.h					\
					order.h					\
					mesh.h					\
					hmesh.h					\
					hmesh_order_si.h

hmesh_order_st$(OBJ)		:	hmesh_order_st.c			\
					module.h				\
					common.h				\
					parser.h				\
					graph.h					\
					hgraph.h				\
					hgraph_order_st.h			\
					mesh.h					\
					hmesh.h					\
					order.h					\
					vmesh_separate_st.h			\
					hmesh_order_bl.h			\
					hmesh_order_gp.h			\
					hmesh_order_gr.h			\
					hmesh_order_hd.h			\
					hmesh_order_hf.h			\
					hmesh_order_nd.h			\
					hmesh_order_si.h			\
					hmesh_order_st.h

kdgraph$(OBJ)			:	kdgraph.c				\
					module.h				\
					common.h				\
					arch$(OBJ)				\
					dgraph.h				\
					dmapping.h				\
					kdgraph.h

kdgraph_gather$(OBJ)		:	kdgraph_gather.c			\
					module.h				\
					common.h				\
					arch$(OBJ)				\
					graph.h					\
					mapping.h				\
					kgraph.h				\
					dgraph.h				\
					kdgraph.h

kdgraph_map_rb$(OBJ)		:	kdgraph_map_rb.c			\
					module.h				\
					common.h				\
					parser.h				\
					arch$(OBJ)				\
					graph.h					\
					dgraph.h				\
					dmapping.h				\
					kdgraph.h				\
					kdgraph_map_rb.h			\
					kdgraph_map_rb_map.h			\
					kdgraph_map_rb_part.h

kdgraph_map_rb_map$(OBJ)	:	kdgraph_map_rb_map.c			\
					module.h				\
					common.h				\
					parser.h				\
					arch$(OBJ)				\
					graph.h					\
					bgraph.h				\
					bgraph_bipart_st.h			\
					mapping.h				\
					kgraph.h				\
					kgraph_map_st.h				\
					dgraph.h				\
					dmapping.h				\
					bdgraph.h				\
					bdgraph_bipart_st.h			\
					kdgraph.h				\
					kdgraph_map_rb.h			\
					kdgraph_map_rb_map.h			\
					kdgraph_map_st.h

kdgraph_map_rb_part$(OBJ)	:	kdgraph_map_rb_part.c			\
					module.h				\
					common.h				\
					parser.h				\
					arch$(OBJ)				\
					graph.h					\
					bgraph.h				\
					bgraph_bipart_st.h			\
					mapping.h				\
					kgraph.h				\
					kgraph_map_st.h				\
					dgraph.h				\
					dmapping.h				\
					bdgraph.h				\
					bdgraph_bipart_st.h			\
					kdgraph.h				\
					kdgraph_map_rb.h			\
					kdgraph_map_rb_part.h			\
					kdgraph_map_st.h

kdgraph_map_st$(OBJ)		:	kdgraph_map_st.c			\
					module.h				\
					common.h				\
					parser.h				\
					arch$(OBJ)				\
					graph.h					\
					dgraph.h				\
					dgraph_coarsen.h			\
					mapping.h				\
					dmapping.h				\
					bdgraph.h				\
					bdgraph_bipart_st.h			\
					kgraph.h				\
					kgraph_map_st.h				\
					kdgraph.h				\
					kdgraph_map_rb.h			\
					kdgraph_map_st.h

kgraph$(OBJ)			:	kgraph.c				\
					module.h				\
					common.h				\
					graph.h					\
					arch$(OBJ)				\
					mapping.h				\
					kgraph.h

kgraph_band$(OBJ)		:	kgraph_band.c		   		\
					module.h				\
					common.h				\
					arch$(OBJ)				\
					graph.h			        	\
					mapping.h				\
					kgraph.h

kgraph_check$(OBJ)		:	kgraph_check.c				\
					module.h				\
					common.h				\
					graph.h					\
					arch$(OBJ)				\
					mapping.h				\
					kgraph.h

kgraph_map_bd$(OBJ)		:	kgraph_map_bd.c				\
					module.h				\
					common.h				\
					parser.h				\
					graph.h					\
					arch$(OBJ)				\
					mapping.h				\
					kgraph.h				\
					kgraph_map_bd.h

kgraph_map_cp$(OBJ)		:	kgraph_map_cp.c				\
					module.h				\
					common.h				\
					parser.h				\
					graph.h					\
					arch$(OBJ)				\
					mapping.h				\
					kgraph.h				\
					kgraph_map_cp.h

kgraph_map_df$(OBJ)		:	kgraph_map_df.c				\
					kgraph_map_df_loop.c			\
					module.h				\
					common.h				\
					parser.h				\
					graph.h					\
					arch$(OBJ)				\
					mapping.h				\
					kgraph.h				\
					kgraph_map_df.h

kgraph_map_ex$(OBJ)		:	kgraph_map_ex.c				\
					module.h				\
					common.h				\
					graph.h					\
					arch$(OBJ)				\
					mapping.h				\
					kgraph.h				\
					kgraph_map_ex.h

kgraph_map_fm$(OBJ)		:	kgraph_map_fm.c				\
					module.h				\
					common.h				\
					parser.h				\
					graph.h					\
					arch$(OBJ)				\
					mapping.h				\
					kgraph.h				\
					kgraph_map_fm.h                         \
					kgraph_map_rb.h

kgraph_map_ml$(OBJ)		:	kgraph_map_ml.c				\
					module.h				\
					common.h				\
					parser.h				\
					graph.h					\
					graph_coarsen.h				\
					arch$(OBJ)				\
					mapping.h				\
					kgraph.h				\
					kgraph_map_ml.h

kgraph_map_rb$(OBJ)		:	kgraph_map_rb.c				\
					module.h				\
					common.h				\
					parser.h				\
					graph.h					\
					arch$(OBJ)				\
					mapping.h				\
					bgraph.h				\
					bgraph_bipart_st.h			\
					kgraph.h				\
					kgraph_map_rb.h				\
					kgraph_map_rb_map.h			\
					kgraph_map_rb_part.h

kgraph_map_rb_map$(OBJ)		:	kgraph_map_rb_map.c			\
					module.h				\
					common.h				\
					parser.h				\
					graph.h					\
					arch$(OBJ)				\
					mapping.h				\
					bgraph.h				\
					bgraph_bipart_st.h			\
					kgraph.h				\
					kgraph_map_rb.h				\
					kgraph_map_rb_map.h

kgraph_map_rb_part$(OBJ)	:	kgraph_map_rb_part.c			\
					module.h				\
					common.h				\
					parser.h				\
					graph.h					\
					arch$(OBJ)				\
					mapping.h				\
					bgraph.h				\
					bgraph_bipart_st.h			\
					kgraph.h				\
					kgraph_map_rb.h				\
					kgraph_map_rb_part.h

kgraph_map_st$(OBJ)		:	kgraph_map_st.c				\
					module.h				\
					common.h				\
					parser.h				\
					graph.h					\
					arch$(OBJ)				\
					mapping.h				\
					bgraph.h				\
					bgraph_bipart_st.h			\
					kgraph.h				\
					kgraph_map_bd.h				\
					kgraph_map_cp.h				\
					kgraph_map_df.h				\
					kgraph_map_fm.h				\
					kgraph_map_ml.h				\
					kgraph_map_rb.h

kgraph_store$(OBJ)		:	kgraph_store.c				\
					module.h				\
					common.h				\
					graph.h					\
					arch$(OBJ)				\
					mapping.h				\
					kgraph.h

library_arch$(OBJ)		:	library_arch.c				\
					module.h				\
					common.h				\
					graph.h					\
					arch$(OBJ)				\
					arch_cmplt.h				\
					scotch.h

library_arch_f$(OBJ)		:	library_arch_f.c			\
					module.h				\
					common.h				\
					scotch.h

library_arch_build$(OBJ)	:	library_arch_build.c			\
					module.h				\
					common.h				\
					parser.h				\
					graph.h					\
					arch$(OBJ)				\
					arch_build.h				\
					mapping.h				\
					bgraph.h				\
					bgraph_bipart_st.h			\
					scotch.h

library_arch_build_f$(OBJ)	:	library_arch_build_f.c			\
					module.h				\
					common.h				\
					scotch.h

library_common_f$(OBJ)		:	library_common_f.c			\
					module.h				\
					common.h				\
					scotch.h

library_dgraph$(OBJ)		:	library_dgraph.c			\
					module.h				\
					common.h				\
					graph.h					\
					dgraph.h				\
					ptscotch.h

library_dgraph_f$(OBJ)		:	library_dgraph_f.c			\
					module.h				\
					common.h				\
					ptscotch.h

library_dgraph_band$(OBJ)	:	library_dgraph_band.c			\
					module.h				\
					common.h				\
					dgraph.h				\
					ptscotch.h

library_dgraph_band_f$(OBJ)	:	library_dgraph_band_f.c			\
					module.h				\
					common.h				\
					ptscotch.h

library_dgraph_build$(OBJ)	:	library_dgraph_build.c			\
					module.h				\
					common.h				\
					dgraph.h				\
					ptscotch.h

library_dgraph_build_f$(OBJ)	:	library_dgraph_build_f.c		\
					module.h				\
					common.h				\
					ptscotch.h

library_dgraph_build_grid3d$(OBJ):	library_dgraph_build_grid3d.c		\
					module.h				\
					common.h				\
					dgraph.h				\
					ptscotch.h

library_dgraph_build_grid3d_f$(OBJ):	library_dgraph_build_grid3d_f.c		\
					module.h				\
					common.h				\
					ptscotch.h

library_dgraph_check$(OBJ)	:	library_dgraph_check.c			\
					module.h				\
					common.h				\
					dgraph.h				\
					ptscotch.h

library_dgraph_check_f$(OBJ)	:	library_dgraph_check_f.c		\
					module.h				\
					common.h				\
					ptscotch.h

library_dgraph_coarsen$(OBJ)	:	library_dgraph_coarsen.c		\
					module.h				\
					common.h				\
					dgraph.h				\
					dgraph_coarsen.h			\
					ptscotch.h

library_dgraph_coarsen_f$(OBJ)	:	library_dgraph_coarsen_f.c		\
					module.h				\
					common.h				\
					ptscotch.h

library_dgraph_gather$(OBJ)	:	library_dgraph_gather.c			\
					module.h				\
					common.h				\
					dgraph.h				\
					ptscotch.h

library_dgraph_gather_f$(OBJ)	:	library_dgraph_gather_f.c		\
					module.h				\
					common.h				\
					ptscotch.h

library_dgraph_grow$(OBJ)	:	library_dgraph_grow.c			\
					module.h				\
					common.h				\
					dgraph.h				\
					ptscotch.h

library_dgraph_halo$(OBJ)	:	library_dgraph_halo.c			\
					module.h				\
					common.h				\
					dgraph.h				\
					dgraph_halo.h				\
					ptscotch.h

library_dgraph_halo_f$(OBJ)	:	library_dgraph_halo_f.c			\
					module.h				\
					common.h				\
					ptscotch.h

library_dgraph_induce$(OBJ)	:	library_dgraph_induce.c			\
					module.h				\
					common.h				\
					dgraph.h				\
					ptscotch.h

library_dgraph_induce_f$(OBJ)	:	library_dgraph_induce_f.c		\
					module.h				\
					common.h				\
					ptscotch.h

library_dgraph_io_load$(OBJ)	:	library_dgraph_io_load.c		\
					module.h				\
					common.h				\
					dgraph.h				\
					ptscotch.h

library_dgraph_io_load_f$(OBJ)	:	library_dgraph_io_load_f.c		\
					module.h				\
					common.h				\
					ptscotch.h

library_dgraph_io_save$(OBJ)	:	library_dgraph_io_save.c		\
					module.h				\
					common.h				\
					dgraph.h				\
					ptscotch.h

library_dgraph_io_save_f$(OBJ)	:	library_dgraph_io_save_f.c		\
					module.h				\
					common.h				\
					ptscotch.h

library_dgraph_map$(OBJ)	:	library_dgraph_map.c			\
					module.h				\
					common.h				\
					parser.h				\
					arch$(OBJ)				\
					dgraph.h				\
					dmapping.h				\
					kdgraph.h				\
					kdgraph_map_st.h			\
					library_dmapping.h			\
					ptscotch.h

library_dgraph_map_f$(OBJ)	:	library_dgraph_map_f.c			\
					module.h				\
					common.h				\
					ptscotch.h

library_dgraph_map_view$(OBJ)	:	library_dgraph_map_view.c		\
					module.h				\
					common.h				\
					parser.h				\
					dgraph.h				\
					dgraph_halo.h				\
					arch$(OBJ)				\
					dmapping.h				\
					kdgraph.h				\
					library_dmapping.h			\
					ptscotch.h

library_dgraph_map_view_f$(OBJ)	:	library_dgraph_map_view_f.c		\
					module.h				\
					common.h				\
					ptscotch.h

library_dgraph_order$(OBJ)	:	library_dgraph_order.c			\
					module.h				\
					common.h				\
					parser.h				\
					dgraph.h				\
					dorder.h				\
					hdgraph.h				\
					hdgraph_order_st.h			\
					ptscotch.h

library_dgraph_order_f$(OBJ)	:	library_dgraph_order_f.c		\
					module.h				\
					common.h				\
					ptscotch.h

library_dgraph_order_gather$(OBJ):	library_dgraph_order_gather.c		\
					module.h				\
					common.h				\
					dgraph.h				\
					order.h					\
					dorder.h				\
					library_order.h				\
					ptscotch.h

library_dgraph_order_gather_f$(OBJ):	library_dgraph_order_gather_f.c		\
					module.h				\
					common.h				\
					ptscotch.h

library_dgraph_order_io$(OBJ)	:	library_dgraph_order_io.c		\
					module.h				\
					common.h				\
					dgraph.h				\
					dorder.h				\
					ptscotch.h

library_dgraph_order_io_f$(OBJ)	:	library_dgraph_order_io_f.c		\
					module.h				\
					common.h				\
					ptscotch.h

library_dgraph_order_io_block$(OBJ):	library_dgraph_order_io_block.c		\
					module.h				\
					common.h				\
					dgraph.h				\
					dorder.h				\
					ptscotch.h

library_dgraph_order_io_block_f$(OBJ):	library_dgraph_order_io_block_f.c	\
					module.h				\
					common.h				\
					ptscotch.h

library_dgraph_order_perm$(OBJ)	:	library_dgraph_order_perm.c		\
					module.h				\
					common.h				\
					dgraph.h				\
					dorder.h				\
					ptscotch.h

library_dgraph_order_perm_f$(OBJ):	library_dgraph_order_perm_f.c		\
					module.h				\
					common.h				\
					ptscotch.h

library_dgraph_order_tree_dist$(OBJ):	library_dgraph_order_tree_dist.c	\
					module.h				\
					common.h				\
					dgraph.h				\
					dorder.h				\
					ptscotch.h

library_dgraph_order_tree_dist_f$(OBJ):	library_dgraph_order_tree_dist_f.c	\
					module.h				\
					common.h				\
					ptscotch.h

library_dgraph_redist$(OBJ)	:	library_dgraph_redist.c			\
					module.h				\
					common.h				\
					dgraph.h				\
					ptscotch.h

library_dgraph_redist_f$(OBJ)	:	library_dgraph_redist_f.c		\
					module.h				\
					common.h				\
					ptscotch.h

library_dgraph_scatter$(OBJ)	:	library_dgraph_scatter.c		\
					module.h				\
					common.h				\
					dgraph.h				\
					ptscotch.h

library_dgraph_scatter_f$(OBJ)	:	library_dgraph_scatter_f.c		\
					module.h				\
					common.h				\
					ptscotch.h

library_dgraph_stat$(OBJ)	:	library_dgraph_stat.c			\
					module.h				\
					common.h				\
					dgraph.h				\
					ptscotch.h				\
					library_dgraph_stat.h

library_dgraph_stat_f$(OBJ)	:	library_dgraph_stat_f.c			\
					module.h				\
					common.h				\
					ptscotch.h

library_dmapping$(OBJ)		:	library_dmapping.c			\
					module.h				\
					common.h				\
					ptscotch.h

library_dorder$(OBJ)		:	library_dorder.c			\
					module.h				\
					common.h				\
					ptscotch.h

library_error$(OBJ)		:	library_error.c				\
					module.h				\
					common.h				\
					scotch.h

library_error_exit$(OBJ)	:	library_error_exit.c			\
					module.h				\
					common.h				\
					scotch.h

library_error_pt$(OBJ)		:	library_error.c				\
					module.h				\
					common.h				\
					scotch.h
					$(CC) $(CFLAGS) $(CLIBFLAGS) -c $(<) -o $(@)

library_error_exit_pt$(OBJ)	:	library_error_exit.c			\
					module.h				\
					common.h				\
					scotch.h
					$(CC) $(CFLAGS) $(CLIBFLAGS) -c $(<) -o $(@)

library_geom$(OBJ)		:	library_geom.c				\
					module.h				\
					common.h				\
					geom.h					\
					graph.h					\
					scotch.h

library_geom_f$(OBJ)		:	library_geom_f.c			\
					module.h				\
					common.h				\
					scotch.h

library_graph$(OBJ)		:	library_graph.c				\
					module.h				\
					common.h				\
					graph.h					\
					graph_io.h				\
					scotch.h

library_graph_f$(OBJ)		:	library_graph_f.c			\
					module.h				\
					common.h				\
					scotch.h

library_graph_base$(OBJ)	:	library_graph_base.c			\
					module.h				\
					common.h				\
					graph.h					\
					scotch.h

library_graph_base_f$(OBJ)	:	library_graph_base_f.c			\
					module.h				\
					common.h				\
					scotch.h

library_graph_check$(OBJ)	:	library_graph_check.c			\
					module.h				\
					common.h				\
					graph.h					\
					scotch.h

library_graph_check_f$(OBJ)	:	library_graph_check_f.c			\
					module.h				\
					common.h				\
					scotch.h

library_graph_coarsen$(OBJ)	:	library_graph_coarsen.c			\
					module.h				\
					common.h				\
					graph.h					\
					graph_coarsen.h				\
					scotch.h

library_graph_coarsen_f$(OBJ)	:	library_graph_coarsen_f.c		\
					module.h				\
					common.h				\
					scotch.h

library_graph_color$(OBJ)	:	library_graph_color.c			\
					module.h				\
					common.h				\
					graph.h					\
					scotch.h

library_graph_color_f$(OBJ)	:	library_graph_color_f.c			\
					module.h				\
					common.h				\
					scotch.h

library_graph_io_chac$(OBJ)	:	library_graph_io_chac.c			\
					module.h				\
					common.h				\
					geom.h					\
					graph.h					\
					scotch.h

library_graph_io_chac_f$(OBJ)	:	library_graph_io_chac_f.c		\
					module.h				\
					common.h				\
					scotch.h

library_graph_io_habo$(OBJ)	:	library_graph_io_habo.c			\
					module.h				\
					common.h				\
					geom.h					\
					graph.h					\
					scotch.h

library_graph_io_habo_f$(OBJ)	:	library_graph_io_habo_f.c		\
					module.h				\
					common.h				\
					scotch.h

library_graph_io_mmkt$(OBJ)	:	library_graph_io_mmkt.c			\
					module.h				\
					common.h				\
					geom.h					\
					graph.h					\
					scotch.h

library_graph_io_mmkt_f$(OBJ)	:	library_graph_io_mmkt_f.c		\
					module.h				\
					common.h				\
					scotch.h

library_graph_io_scot$(OBJ)	:	library_graph_io_scot.c			\
					module.h				\
					common.h				\
					geom.h					\
					graph.h					\
					scotch.h

library_graph_io_scot_f$(OBJ)	:	library_graph_io_scot_f.c		\
					module.h				\
					common.h				\
					scotch.h

library_graph_map$(OBJ)		:	library_graph_map.c			\
					module.h				\
					common.h				\
					parser.h				\
					graph.h					\
					arch$(OBJ)				\
					mapping.h				\
					kgraph.h				\
					kgraph_map_st.h				\
					library_mapping.h			\
					scotch.h

library_graph_map_f$(OBJ)	:	library_graph_map_f.c			\
					module.h				\
					common.h				\
					scotch.h

library_graph_map_io$(OBJ)	:	library_graph_map_io.c			\
					module.h				\
					common.h				\
					graph.h					\
					arch$(OBJ)				\
					library_mapping.h			\
					library_graph_map_io.h			\
					scotch.h

library_graph_map_io_f$(OBJ)	:	library_graph_map_io_f.c		\
					module.h				\
					common.h				\
					scotch.h

library_graph_map_view$(OBJ)	:	library_graph_map_view.c		\
					module.h				\
					common.h				\
					parser.h				\
					graph.h					\
					arch$(OBJ)				\
					mapping.h				\
					kgraph.h				\
					library_mapping.h			\
					library_graph_map_view.h		\
					scotch.h

library_graph_map_view_f$(OBJ)	:	library_graph_map_view_f.c		\
					module.h				\
					common.h				\
					scotch.h

library_graph_order$(OBJ)	:	library_graph_order.c			\
					module.h				\
					common.h				\
					parser.h				\
					graph.h					\
					order.h					\
					hgraph.h				\
					hgraph_order_st.h			\
					library_order.h				\
					scotch.h

library_graph_order_f$(OBJ)	:	library_graph_order_f.c			\
					module.h				\
					common.h				\
					scotch.h

library_graph_part_ovl$(OBJ)	:	library_graph_part_ovl.c		\
					module.h				\
					common.h				\
					parser.h				\
					graph.h					\
					wgraph.h				\
					wgraph_part_st.h			\
					scotch.h

library_graph_part_ovl_f$(OBJ)	:	library_graph_part_ovl_f.c		\
					module.h				\
					common.h				\
					scotch.h

library_mapping$(OBJ)		:	library_mapping.c			\
					module.h				\
					common.h				\
					scotch.h

library_memory$(OBJ)		:	library_memory.c			\
					module.h				\
					common.h				\
					scotch.h

library_memory_f$(OBJ)		:	library_memory_f.c			\
					module.h				\
					common.h				\
					scotch.h

library_mesh$(OBJ)		:	library_mesh.c				\
					module.h				\
					common.h				\
					graph.h					\
					mesh.h					\
					scotch.h

library_mesh_f$(OBJ)		:	library_mesh_f.c			\
					module.h				\
					common.h				\
					scotch.h

library_mesh_graph$(OBJ)	:	library_mesh_graph.c			\
					module.h				\
					common.h				\
					graph.h					\
					mesh.h					\
					scotch.h

library_mesh_graph_f$(OBJ)	:	library_mesh_graph_f.c			\
					module.h				\
					common.h				\
					scotch.h

library_mesh_io_habo$(OBJ)	:	library_mesh_io_habo.c			\
					module.h				\
					common.h				\
					geom.h					\
					graph.h					\
					mesh.h					\
					scotch.h

library_mesh_io_habo_f$(OBJ)	:	library_mesh_io_habo_f.c		\
					module.h				\
					common.h				\
					scotch.h

library_mesh_io_scot$(OBJ)	:	library_mesh_io_scot.c			\
					module.h				\
					common.h				\
					geom.h					\
					graph.h					\
					mesh.h					\
					scotch.h

library_mesh_io_scot_f$(OBJ)	:	library_mesh_io_scot_f.c		\
					module.h				\
					common.h				\
					scotch.h

library_mesh_order$(OBJ)	:	library_mesh_order.c			\
					module.h				\
					common.h				\
					parser.h				\
					graph.h					\
					mesh.h					\
					hmesh.h					\
					order.h					\
					hmesh_order_st.h			\
					scotch.h				\
					library_order.h

library_mesh_order_f$(OBJ)	:	library_mesh_order_f.c			\
					module.h				\
					common.h				\
					scotch.h

library_order$(OBJ)		:	library_order.c				\
					module.h				\
					common.h				\
					scotch.h

library_parser$(OBJ)		:	library_parser.c			\
					module.h				\
					common.h				\
					parser.h				\
					scotch.h

library_parser_f$(OBJ)		:	library_parser_f.c			\
					module.h				\
					common.h				\
					scotch.h

library_random$(OBJ)		:	library_random.c			\
					module.h				\
					common.h				\
					scotch.h

library_random_f$(OBJ)		:	library_random_f.c			\
					module.h				\
					common.h				\
					scotch.h

library_strat$(OBJ)		:	library_strat.c				\
					module.h				\
					common.h				\
					scotch.h

library_version$(OBJ)		:	library_version.c			\
					module.h				\
					common.h				\
					scotch.h

library_version_f$(OBJ)		:	library_version_f.c			\
					module.h				\
					common.h				\
					scotch.h

mapping$(OBJ)			:	mapping.c				\
					module.h				\
					common.h				\
					graph.h					\
					arch$(OBJ)				\
					mapping.h

mapping_io$(OBJ)		:	mapping_io.c				\
					module.h				\
					common.h				\
					graph.h					\
					arch$(OBJ)				\
					mapping.h				\
					mapping_io.h

mesh$(OBJ)			:	mesh.c					\
					module.h				\
					common.h				\
					graph.h					\
					mesh.h

mesh_check$(OBJ)		:	mesh_check.c				\
					module.h				\
					common.h				\
					graph.h					\
					mesh.h

mesh_coarsen$(OBJ)		:	mesh_coarsen.c				\
					module.h				\
					common.h				\
					graph.h					\
					mesh.h

mesh_graph$(OBJ)		:	mesh_graph.c				\
					module.h				\
					common.h				\
					graph.h					\
					mesh.h					\
					mesh_graph.h

mesh_induce_sepa$(OBJ)		:	mesh_induce_sepa.c			\
					module.h				\
					common.h				\
					graph.h					\
					mesh.h					\
					mesh_induce_sepa.h

mesh_io$(OBJ)			:	mesh_io.c				\
					module.h				\
					common.h				\
					graph.h					\
					graph_io.h				\
					mesh.h					\
					mesh_io.h

mesh_io_habo$(OBJ)		:	mesh_io_habo.c				\
					module.h				\
					common.h				\
					geom.h					\
					graph.h					\
					mesh.h

mesh_io_scot$(OBJ)		:	mesh_io_scot.c				\
					module.h				\
					common.h				\
					geom.h					\
					graph.h					\
					mesh.h

order$(OBJ)			:	order.c					\
					module.h				\
					common.h				\
					graph.h					\
					order.h

order_check$(OBJ)		:	order_check.c				\
					module.h				\
					common.h				\
					graph.h					\
					order.h

order_io$(OBJ)			:	order_io.c				\
					module.h				\
					common.h				\
					graph.h					\
					order.h

order_tree$(OBJ)		:	order_tree.c				\
					module.h				\
					common.h				\
					graph.h					\
					order.h

parser$(OBJ)			:	parser.c				\
					module.h				\
					common.h				\
					parser.h				\
					parser_yy.h

parser_ll.c			:	parser_ll.l				\
					module.h				\
					common.h				\
					parser.h				\
					parser_ll.h				\
					parser_ly.h
					($(LEX) parser_ll.l &&			\
					 $(MV) lex.yy.c parser_ll.c) ||		\
					$(CP) last_resort/parser_ll.c . 

parser_ll$(OBJ)			:	parser_ll.c				\
					parser_ly.h

parser_ly.h			:	parser_yy$(OBJ)

parser_yy.c			:	parser_yy.y				\
					module.h				\
					common.h				\
					parser.h				\
					parser_yy.h
					($(YACC) -d -v parser_yy.y &&		\
					 $(MV) y.tab.c parser_yy.c &&		\
					 $(MV) y.tab.h parser_ly.h) ||		\
					$(CP) last_resort/parser_ly.h last_resort/parser_yy.c .

parser_yy$(OBJ)			:	parser_yy.c

vdgraph$(OBJ)			:	vdgraph.c				\
					module.h				\
					common.h				\
					dgraph.h				\
					vdgraph.h

vdgraph_check$(OBJ)		:	vdgraph_check.c				\
					module.h				\
					common.h				\
					dgraph.h				\
					vdgraph.h

vdgraph_gather_all$(OBJ)	:	vdgraph_gather_all.c			\
					module.h				\
					common.h				\
					comm.h					\
					graph.h					\
					vgraph.h				\
					dgraph.h				\
					vdgraph.h

vdgraph_separate_bd$(OBJ)	:	vdgraph_separate_bd.c			\
					module.h				\
					common.h				\
					parser.h				\
					dgraph.h				\
					vdgraph.h				\
					vdgraph_separate_bd.h			\
					vdgraph_separate_st.h

vdgraph_separate_df$(OBJ)	:	vdgraph_separate_df.c			\
					module.h				\
					common.h				\
					dgraph.h				\
					vdgraph.h				\
					vdgraph_separate_df.h

vdgraph_separate_ml$(OBJ)	:	vdgraph_separate_ml.c			\
					module.h				\
					common.h				\
					parser.h				\
					dgraph.h				\
					dgraph_coarsen.h			\
					vdgraph.h				\
					vdgraph_separate_ml.h			\
					vdgraph_separate_st.h

vdgraph_separate_sq$(OBJ)	:	vdgraph_separate_sq.c			\
					module.h				\
					common.h				\
					comm.h					\
					parser.h				\
					graph.h					\
					vgraph.h				\
					vgraph_separate_st.h			\
					dgraph.h				\
					vdgraph.h				\
					vdgraph_separate_sq.h

vdgraph_separate_st$(OBJ)	:	vdgraph_separate_st.c			\
					module.h				\
					common.h				\
					parser.h				\
					graph.h					\
					vgraph.h				\
					vgraph_separate_st.h			\
					dgraph.h				\
					dgraph_coarsen.h			\
					vdgraph.h				\
					vdgraph_separate_bd.h			\
					vdgraph_separate_df.h			\
					vdgraph_separate_ml.h			\
					vdgraph_separate_sq.h			\
					vdgraph_separate_st.h			\
					vdgraph_separate_zr.h

vdgraph_separate_zr$(OBJ)	:	vdgraph_separate_zr.c			\
					module.h				\
					common.h				\
					dgraph.h				\
					vdgraph.h				\
					vdgraph_separate_zr.h

vdgraph_store$(OBJ)		:	vdgraph_store.c				\
					module.h				\
					common.h				\
					dgraph.h				\
					vdgraph.h

vgraph$(OBJ)			:	vgraph.c				\
					module.h				\
					common.h				\
					graph.h					\
					vgraph.h

vgraph_check$(OBJ)		:	vgraph_check.c				\
					module.h				\
					common.h				\
					graph.h					\
					vgraph.h

vgraph_separate_bd$(OBJ)	:	vgraph_separate_bd.c			\
					module.h				\
					common.h				\
					parser.h				\
					graph.h					\
					vgraph.h				\
					vgraph_separate_bd.h			\
					vgraph_separate_st.h

vgraph_separate_df$(OBJ)	:	vgraph_separate_df.c			\
					module.h				\
					common.h				\
					graph.h					\
					vgraph.h				\
					vgraph_separate_df.h

vgraph_separate_es$(OBJ)	:	vgraph_separate_es.c			\
					module.h				\
					common.h				\
					parser.h				\
					graph.h					\
					arch$(OBJ)				\
					mapping.h				\
					bgraph.h				\
					bgraph_bipart_st.h			\
					vgraph.h				\
					vgraph_separate_es.h

vgraph_separate_fm$(OBJ)	:	vgraph_separate_fm.c			\
					module.h				\
					common.h				\
					gain.h					\
					graph.h					\
					vgraph.h				\
					vgraph_separate_gg.h			\
					vgraph_separate_fm.h

vgraph_separate_gg$(OBJ)	:	vgraph_separate_gg.c			\
					module.h				\
					common.h				\
					graph.h					\
					vgraph.h				\
					vgraph_separate_gg.h

vgraph_separate_gp$(OBJ)	:	vgraph_separate_gp.c			\
					module.h				\
					common.h				\
					graph.h					\
					vgraph.h				\
					vgraph_separate_gp.h

vgraph_separate_ml$(OBJ)	:	vgraph_separate_ml.c			\
					module.h				\
					common.h				\
					parser.h				\
					graph.h					\
					graph_coarsen.h				\
					vgraph.h				\
					vgraph_separate_ml.h			\
					vgraph_separate_st.h

vgraph_separate_th$(OBJ)	:	vgraph_separate_th.c			\
					module.h				\
					common.h				\
					graph.h					\
					vgraph.h				\
					vgraph_separate_th.h

vgraph_separate_vw$(OBJ)	:	vgraph_separate_vw.c			\
					module.h				\
					common.h				\
					graph.h					\
					vgraph.h				\
					vgraph_separate_vw.h

vgraph_separate_zr$(OBJ)	:	vgraph_separate_zr.c			\
					module.h				\
					common.h				\
					graph.h					\
					vgraph.h				\
					vgraph_separate_zr.h

vgraph_separate_st$(OBJ)	:	vgraph_separate_st.c			\
					gain.h					\
					module.h				\
					common.h				\
					parser.h				\
					graph.h					\
					arch$(OBJ)				\
					mapping.h				\
					bgraph.h				\
					bgraph_bipart_st.h			\
					vgraph.h				\
					vgraph_separate_bd.h			\
					vgraph_separate_df.h			\
					vgraph_separate_fm.h			\
					vgraph_separate_gg.h			\
					vgraph_separate_gp.h			\
					vgraph_separate_ml.h			\
					vgraph_separate_th.h			\
					vgraph_separate_vw.h			\
					vgraph_separate_zr.h			\
					vgraph_separate_st.h

vgraph_store$(OBJ)		:	vgraph_store.c				\
					module.h				\
					common.h				\
					graph.h					\
					vgraph.h

vmesh$(OBJ)			:	vmesh.c					\
					module.h				\
					common.h				\
					graph.h					\
					mesh.h					\
					vmesh.h

vmesh_check$(OBJ)		:	vmesh_check.c				\
					module.h				\
					common.h				\
					graph.h					\
					mesh.h					\
					vmesh.h

vmesh_separate_fm$(OBJ)		:	vmesh_separate_fm.c			\
					module.h				\
					common.h				\
					graph.h					\
					mesh.h					\
					vmesh.h					\
					vmesh_separate_fm.h

vmesh_separate_gg$(OBJ)		:	vmesh_separate_gg.c			\
					module.h				\
					common.h				\
					graph.h					\
					mesh.h					\
					vmesh.h					\
					vmesh_separate_gg.h

vmesh_separate_gr$(OBJ)		:	vmesh_separate_gr.c			\
					module.h				\
					common.h				\
					parser.h				\
					graph.h					\
					vgraph.h				\
					vgraph_separate_st.h			\
					mesh.h					\
					vmesh.h					\
					vmesh_separate_gr.h

vmesh_separate_ml$(OBJ)		:	vmesh_separate_ml.c			\
					module.h				\
					common.h				\
					parser.h				\
					graph.h					\
					mesh.h					\
					mesh_coarsen.h				\
					vmesh.h					\
					vmesh_separate_ml.h			\
					vmesh_separate_st.h

vmesh_separate_zr$(OBJ)		:	vmesh_separate_zr.c			\
					module.h				\
					common.h				\
					graph.h					\
					mesh.h					\
					vmesh.h					\
					vmesh_separate_zr.h

vmesh_separate_st$(OBJ)		:	vmesh_separate_st.c			\
					gain.h					\
					module.h				\
					common.h				\
					parser.h				\
					graph.h					\
					mesh.h					\
					mesh_coarsen.h				\
					vmesh.h					\
					vmesh_separate_fm.h			\
					vmesh_separate_gg.h			\
					vmesh_separate_gr.h			\
					vmesh_separate_ml.h			\
					vmesh_separate_zr.h			\
					vmesh_separate_st.h

vmesh_store$(OBJ)		:	vmesh_store.c				\
					module.h				\
					common.h				\
					graph.h					\
					mesh.h					\
					vmesh.h

wgraph$(OBJ)			:	wgraph.c				\
					module.h				\
					common.h				\
					graph.h					\
					wgraph.h

wgraph_check$(OBJ)		:	wgraph_check.c				\
					module.h				\
					common.h				\
					graph.h					\
					wgraph.h

wgraph_part_fm$(OBJ)		:	wgraph_part_fm.c			\
					module.h				\
					common.h				\
					graph.h					\
					wgraph.h				\
					wgraph_part_fm.h

wgraph_part_gg$(OBJ)		:	wgraph_part_gg.c			\
					module.h				\
					common.h				\
					graph.h					\
					wgraph.h				\
					wgraph_part_gg.h

wgraph_part_gp$(OBJ)		:	wgraph_part_gp.c			\
					module.h				\
					common.h				\
					graph.h					\
					wgraph.h				\
					wgraph_part_gp.h

wgraph_part_ml$(OBJ)		: 	wgraph_part_ml.c			\
					module.h				\
					common.h				\
					parser.h				\
					graph.h					\
					graph_coarsen.h				\
					wgraph.h				\
					wgraph_part_ml.h			\
					wgraph_part_st.h

wgraph_part_rb$(OBJ)		:	wgraph_part_rb.c			\
					module.h				\
					common.h                                \
					parser.h                                \
					graph.h                                 \
					arch.h                                  \
					mapping.h                               \
					vgraph.h                                \
					vgraph_separate_st.h                    \
					vgraph_separate_zr.h                    \
					wgraph.h                                \
					wgraph_part_rb.h			\
					scotch.h

wgraph_part_st$(OBJ)		:	wgraph_part_st.c			\
					gain.h					\
					module.h				\
					common.h				\
					parser.h				\
					graph.h					\
					graph_coarsen.h				\
					vgraph.h				\
					vgraph_separate_st.h			\
					wgraph.h				\
					wgraph_part_fm.h			\
					wgraph_part_gg.h			\
					wgraph_part_gp.h			\
					wgraph_part_ml.h			\
					wgraph_part_rb.h			\
					wgraph_part_zr.h			\
					wgraph_part_st.h

wgraph_part_zr$(OBJ)		: 	wgraph_part_zr.c			\
					module.h				\
					common.h				\
					graph.h					\
					wgraph.h				\
					wgraph_part_zr.h

wgraph_store$(OBJ)		:	wgraph_store.c				\
					module.h				\
					common.h				\
					graph.h					\
					wgraph.h

dummysizes$(EXE)		:	dummysizes.c				\
					module.h				\
					common.h				\
					arch$(OBJ)				\
					graph.h					\
					geom.h					\
					mesh.h					\
					mapping.h				\
					order.h					\
					parser.h
			 		$(CCD) $(CFLAGS) -DSCOTCH_VERSION=$(VERSION) -DSCOTCH_RELEASE=$(RELEASE) -DSCOTCH_PATCHLEVEL=$(PATCHLEVEL) $(<) -o $(@) $(LDFLAGS)

ptdummysizes$(EXE)		:	dummysizes.c				\
					module.h				\
					common.h				\
					dgraph.h				\
					dorder.h
			 		$(CCD) $(CFLAGS) -DSCOTCH_VERSION=$(VERSION) -DSCOTCH_RELEASE=$(RELEASE) -DSCOTCH_PATCHLEVEL=$(PATCHLEVEL) $(<) -o $(@) $(LDFLAGS)

scotch.h			:	dummysizes$(EXE)			\
					library.h
					./dummysizes$(EXE) library.h scotch.h

scotchf.h			:	dummysizes$(EXE) 			\
					library_f.h
					./dummysizes$(EXE) library_f.h scotchf.h

ptscotch.h			:	ptdummysizes$(EXE)			\
					library_pt.h
					./ptdummysizes$(EXE) library_pt.h ptscotch.h

ptscotchf.h			:	ptdummysizes$(EXE)			\
					library_pt_f.h
					./ptdummysizes$(EXE) library_pt_f.h ptscotchf.h

libscotch$(LIB)			:	$(LIBSCOTCHDEPS)
					$(AR) $(ARFLAGS) $(@) $(?)
					-$(RANLIB) $(@)

libscotcherr$(LIB)		:	library_error$(OBJ)
					$(AR) $(ARFLAGS) $(@) $(?)
					-$(RANLIB) $(@)

libscotcherrexit$(LIB)		:	library_error_exit$(OBJ)
					$(AR) $(ARFLAGS) $(@) $(?)
					-$(RANLIB) $(@)

libptscotch$(LIB)		:	$(LIBPTSCOTCHDEPS)
					$(AR) $(ARFLAGS) $(@) $(?)
					-$(RANLIB) $(@)

libptscotcherr$(LIB)		:	library_error_pt$(OBJ)
					$(AR) $(ARFLAGS) $(@) $(?)
					-$(RANLIB) $(@)

libptscotcherrexit$(LIB)	:	library_error_exit_pt$(OBJ)
					$(AR) $(ARFLAGS) $(@) $(?)
					-$(RANLIB) $(@)
