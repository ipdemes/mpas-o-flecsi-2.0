/*----------------------------------------------------------------------------*
  Copyright (c) 2020 Triad National Security, LLC
  All rights reserved
 *----------------------------------------------------------------------------*/

#pragma once

#include <cstddef>

namespace mpas {
static constexpr std::size_t maxEdges = 8;
static constexpr std::size_t maxEdges2 = maxEdges * 2;
static constexpr std::size_t nAdvectionCells = maxEdges2;
static constexpr std::size_t vertexDegree = 3;
//static constexpr std::size_t nVertLevels = 20;
static constexpr std::size_t nVertLevels = 5;
//static constexpr std::size_t nVertLevels = 1;
static constexpr std::size_t maxTracers = 5;
//static constexpr std::size_t maxTracers = 1;
static constexpr std::size_t deuce = 2;
}

namespace mpas_constants {

// The Sphere_240QU.nc mesh has a radius of 6371000.0 m.
// static constexpr double a       = 6371229.0;          // Spherical Earth
// radius [m]
static constexpr double a = 6371000.0; // Spherical Earth radius [m]
static constexpr double pi = 3.141592653589793;
static constexpr double gravity =
  9.80616; // Acceleration due to gravity [m s-2]
static constexpr double omega =
  7.29212e-5; // Angular rotation rate of the Earth [s-1]
static constexpr double rho_sw = 1026.0; // density of salt water (kg/m^3)
}
