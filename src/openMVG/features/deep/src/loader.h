#ifndef CVPR15_MATCHER_LOADER_H
#define CVPR15_MATCHER_LOADER_H
// Copyright 2015 Sergey Zagoruyko, Nikos Komodakis
// sergey.zagoruyko@imagine.enpc.fr, nikos.komodakis@enpc.fr
// Ecole des Ponts ParisTech, Universite Paris-Est, IMAGINE
//
// The software is free to use only for non-commercial purposes.
// IF YOU WOULD LIKE TO USE IT FOR COMMERCIAL PURPOSES, PLEASE CONTACT
// Prof. Nikos Komodakis (nikos.komodakis@enpc.fr)
#include "../cunnproduction/cunn.h"

cunn::Sequential::Ptr loadNetwork(THCState* state, const char* filename);

#endif
