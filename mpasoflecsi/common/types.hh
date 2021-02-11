/*----------------------------------------------------------------------------*
  Copyright (c) 2020 Triad National Security, LLC
  All rights reserved
 *----------------------------------------------------------------------------*/

#pragma once

#include <array>

#include "constants.hh"

namespace mpas {

template<class T>
using vdtensor = std::array<T, vertexDegree>;

template<class T>
using vltensor = std::array<T, nVertLevels>;

template<class T>
using metensor = std::array<T, maxEdges>;

}
