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

inline cartmesh::slot m;
inline cartmesh::cslot coloring;

} // namespace coupler
