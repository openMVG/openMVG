# RPM spec file for Scotch, tested on RHEL 5 & 6.
# Licence:  CeCILL-C

# Build PT-Scotch stuff with openmpi.  Probably there should be an mpich
# version too, but we don't use mpich.
# The mpiwrappers-openmpi-devel required on RHEL5 is in the EPEL repo.
%bcond_with openmpi

# fixme:  Is there a better way?
%define rhel5 %(fgrep -q ' 5.' /etc/redhat-release 2>/dev/null && echo 1)

Name:		scotch
# The original tar file was called 5.1.12b, but unpacks into 5.1.12, so it
# was re-named.
Version:	5.1.12
Release:	1%{?dist}
Summary:	programs and libraries for graph, mesh and hypergraph partitioning
Group:		Applications/Engineering
License:	CeCILL-C
URL:		https://gforge.inria.fr/projects/scotch
Source0:	https://gforge.inria.fr/frs/download.php/28925/scotch_%{version}.tar.gz
BuildRoot:	%{_tmppath}/%{name}-%{version}-%{release}-root-%(%{__id_u} -n)
BuildRequires:	gcc, make, binutils, bison, flex
%if %{with openmpi}
%if 0%rhel5
BuildRequires:	mpiwrappers-openmpi-devel
%else
BuildRequires:	openmpi-devel
%endif
%endif

# Cribbed from the Debian package
%description
Scotch is a software package for graph and mesh/hypergraph partitioning,
graph clustering, and sparse matrix ordering.

Its purpose is to apply graph theory, with a divide and conquer
approach, to scientific computing problems such as graph and mesh
partitioning, static mapping, and sparse matrix ordering, in
application domains ranging from structural mechanics to operating
systems or bio-chemistry.

The SCOTCH distribution is a set of programs and libraries which
implement the static mapping and sparse matrix reordering algorithms
developed within the SCOTCH project.

SCOTCH has many interesting features:

o Its capabilities can be used through a set of stand-alone programs
as well as through the libSCOTCH library, which offers both C and
Fortran interfaces.

o It provides algorithms to partition graph structures, as well as
mesh structures defined as node-element bipartite graphs and which
can also represent hypergraphs.

o It can map any weighted source graph onto any weighted target
graph. The source and target graphs may have any topology, and their
vertices and edges may be weighted. Moreover, both source and target
graphs may be disconnected. This feature allows for the mapping of
programs onto disconnected subparts of a parallel architecture made
up of heterogeneous processors and communication links.

o It computes amalgamated block orderings of sparse matrices, for
efficient solving using BLAS routines.

o Its running time is linear in the number of edges of the source
graph, and logarithmic in the number of vertices of the target graph
for mapping computations.

o It can handle indifferently graph and mesh data structures created
within C or Fortran programs, with array indices starting from 0 or
1.

o It offers extended support for adaptive graphs and meshes through
the handling of disjoint edge arrays.

o It is dynamically parametrizable thanks to strategy strings that
are interpreted at run-time.

o It uses system memory efficiently, to process large graphs and
meshes without incurring out-of-memory faults;

o It can be easily interfaced to other programs. The programs
comprising the SCOTCH project have been designed to run in
command-line mode without any interactive prompting, so that they can
be called easily from other programs by means of system() or popen()
calls, or piped together on a single command line. Moreover, vertex
labeling capabilities allow for easy renumbering of vertices.

o It provides many tools to build, check, and display graphs, meshes
and matrix patterns.

%package libs
Summary: Shared library files for scotch.
Group:		Applications/Engineering

%description libs
Scotch is a software package for graph and mesh/hypergraph partitioning,
graph clustering, and sparse matrix ordering.
This package contains its library files.

%package devel
Summary: Development files for the scotch library.
Group:		Applications/Engineering
Requires:	%{name}-libs = %{version}-%{release}	

%description devel
Scotch is a software package for graph and mesh/hypergraph partitioning,
graph clustering, and sparse matrix ordering.
This package contains development files.

%if %{with openmpi}
%package pt
Summary: PT-Scotch (parallel programs).
Group:		Applications/Engineering

%description pt
Scotch is a software package for graph and mesh/hypergraph partitioning,
graph clustering, and sparse matrix ordering.
This package contains the MPI-based programs.

%package ptlibs
Summary: Shared library files for ptscotch.
Group:		Applications/Engineering

%description ptlibs
Scotch is a software package for graph and mesh/hypergraph partitioning,
graph clustering, and sparse matrix ordering.
This package contains the MPI-based libraries.
%endif

%prep
%setup -q -n %{name}_%{version}
cat >src/Makefile.inc <<EOF
# Taken from the OpenFOAM Thirdparty distribution.
# Differs from the distributed Makefile.inc.i686_pc_linux2.nothreads in
# the lines noted.
EXE		=
LIB		= .so
OBJ		= .o
MAKE		= make
CAT		= cat
CCS		= gcc
LDFLAGS		= -lz -lm -lrt
CP		= cp
LN		= ln
MKDIR		= mkdir
MV		= mv
# modified:
AR		= gcc
ARFLAGS		= -shared -o
CFLAGS		= -O3 -DCOMMON_FILE_COMPRESS_GZ -DCOMMON_RANDOM_FIXED_SEED -DSCOTCH_RENAME -DSCOTCH_RENAME_PARSER -Drestrict=__restrict
CLIBFLAGS	= -fPIC
LEX		= flex -Pscotchyy -olex.yy.c
RANLIB		= echo
YACC		= bison -pscotchyy -y -b y
%if %{with openmpi}
CCP		= mpicc
CCD		= mpicc
%endif
EOF

# We do need to define all the ...dir
%define MAKEOPTS %{?_smp_mflags} prefix=$RPM_BUILD_ROOT bindir=$RPM_BUILD_ROOT%{_bindir} libdir=$RPM_BUILD_ROOT%{_libdir} includedir=$RPM_BUILD_ROOT%{_includedir} datarootdir=$RPM_BUILD_ROOT%{_datadir}

%build
rm -rf $RPM_BUILD_ROOT
mkdir -p $RPM_BUILD_ROOT
# These aren't getting made by make
(cd $RPM_BUILD_ROOT
  mkdir -p .%{_bindir} .%{_libdir} .%{_includedir} .%{_datadir} .%{_mandir};)
%if %{with openmpi}
%{_openmpi_load}
%endif
cd src
make %MAKEOPTS scotch
# ptscoth needs the headers installed
make %MAKEOPTS install
%if %{with openmpi}
make %MAKEOPTS ptscotch
%endif

# installation is done by the build stanza

%clean
rm -rf $RPM_BUILD_ROOT

%files
%defattr(-,root,root,-)
%{_bindir}/[agm]*
%exclude %{_mandir}/man1/d*
%{_mandir}/man1

%files libs
%defattr(-,root,root,-)
%{_libdir}/libscotch*
%doc doc/CeCILL-C_V1-en.txt

%if %{with openmpi}
%files pt
%defattr(-,root,root,-)
%{_bindir}/d*
%{_mandir}/man1/d*

%files ptlibs
%defattr(-,root,root,-)
%{_libdir}/libpt*
%doc doc/CeCILL-C_V1-en.txt
%doc doc/ptscotch_user5.1.pdf
%endif

%files devel
%defattr(-,root,root,-)
%{_includedir}
%doc doc/ptscotch_user5.1.pdf

%changelog
* Tue Oct 23 2012 Dave Love <d.love@liverpool.ac.uk> - 5.1.12b-1
- Initial packaging
