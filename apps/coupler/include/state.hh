/*----------------------------------------------------------------------------*
  Copyright (c) 2020 Triad National Security, LLC
  All rights reserved
 *----------------------------------------------------------------------------*/

#pragma once

#include "mpasoflecsi/specialization/cartmesh.hh"

using namespace mpas; 
using namespace flecsi; 

namespace coupler {

inline const field<double>::definition<cartmesh, cartmesh::cells> ud;
inline const field<double>::definition<cartmesh, cartmesh::cells> vd;

inline cartmesh::slot mesh_src;
inline cartmesh::cslot coloring_src;
inline cartmesh::slot mesh_trg;
inline cartmesh::cslot coloring_trg;
} // namespace coupler
