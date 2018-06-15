#! /bin/sh
# (C) 2008 Yves Secretan (yves.secretan@ete.inrs.ca)
# This software is governed by the CeCILL-C license under French law
# and abiding by the rules of distribution of free software. You can
# use, modify and/or redistribute the software under the terms of the
# CeCILL-C license as circulated by CEA, CNRS and INRIA at the following
# URL: "http://www.cecill.info".
#
# To be executed in a MSYS window.
#
# This file creates the Exports definition file from the PT-Scotch DLL.
# It must be adapted to reflect your environment, in particular library
# path and library name.

echo EXPORTS > ../../../lib/libptscotch.def
nm ../../../lib/libptscotch.a | grep ' T _SCOTCH_' | sed 's/.* T _//' >> ../../../lib/libptscotch.def
