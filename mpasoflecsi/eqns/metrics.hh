/*----------------------------------------------------------------------------*
  Copyright (c) 2020 Triad National Security, LLC
  All rights reserved
 *----------------------------------------------------------------------------*/

#pragma once

#include <cmath>

namespace mpas {
namespace eqns {

template<class REAL>
REAL
sphere_distance(REAL lat1, REAL lon1, REAL lat2, REAL lon2, REAL radius) {
  using namespace std;
  auto arg1 = sqrt(pow(sin(0.5 * (lat2 - lat1)), 2) +
                   cos(lat1) * cos(lat2) * pow(sin(0.5 * (lon2 - lon1)), 2));
  return 2 * radius * asin(arg1);
}

}}
