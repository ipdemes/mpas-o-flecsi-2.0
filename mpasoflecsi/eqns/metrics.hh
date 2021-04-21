/*----------------------------------------------------------------------------*
  Copyright (c) 2020 Triad National Security, LLC
  All rights reserved
 *----------------------------------------------------------------------------*/

#pragma once

#include <cmath>
#include <flecsi/flog.hh>

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

template<class REAL>
REAL
sphere_angle(
  REAL ax, REAL ay, REAL az,
  REAL bx, REAL by, REAL bz,
  REAL cx, REAL cy, REAL cz)
{
  using namespace std;

  // angular side lengths of spherical triangle ABC
  REAL a, b, c;

  a = acos(max(min(bx * cx + by * cy + bz * cz, 1.0), -1.0));
  b = acos(max(min(ax * cx + ay * cy + az * cz, 1.0), -1.0));
  c = acos(max(min(ax * bx + ay * by + az * bz, 1.0), -1.0));

  // components of vectors AB and AC
  REAL ABx, ABy, ABz, ACx, ACy, ACz;
  ABx = bx - ax;
  ABy = by - ay;
  ABz = bz - az;

  ACx = cx - ax;
  ACy = cy - ay;
  ACz = cz - az;

  // cross-product AB x AC
  REAL Dx, Dy, Dz;
  Dx = (ABy * ACz) - (ABz * ACy);
  Dy = (ABz * ACx) - (ABx * ACz);
  Dz = (ABx * ACy) - (ABy * ACx);

  // semiperimeter of triangle ABC
  REAL s;
  s = 0.5 * (a + b + c);

  // sine of desired angle
  REAL sin_angle;
  sin_angle = sqrt(
    min(1.0, max(0.0, (sin(s - b) * sin(s - c)) / (sin(b) * sin(c)))));

  REAL sphere_angle;

  if(Dx * ax + Dy * ay + Dz * az >= 0.0)
    sphere_angle = 2.0 * asin(max(min(sin_angle, 1.0), -1.0));
  else
    sphere_angle = -2.0 * asin(max(min(sin_angle, 1.0), -1.0));

  return sphere_angle;

}

template<class REAL>
REAL
sphere_arc_length(
  REAL ax, REAL ay, REAL az,
  REAL bx, REAL by, REAL bz)
{
  
  REAL cx, cy, cz;
  cx = bx - ax;
  cy = by - ay;
  cz = bz - az;
    
  REAL r = sqrt(ax * ax + ay * ay + az * az);
  REAL c = sqrt(cx * cx + cy * cy + cz * cz);
  
  return r * 2. * asin(c / (2. * r));
}

inline void arc_bisect(
  double ax, double ay, double az,
  double bx, double by, double bz,
  double & cx, double & cy, double & cz)
{
  double r = sqrt(ax*ax + ay*ay + az*az);

  cx = 0.5 * (ax + bx);
  cy = 0.5 * (ay + by);
  cz = 0.5 * (az + bz);

  flog_assert(cx != 0. && cy != 0. && cz != 0.,
    "arc_bisect: A and B are diametrically opposite");

  double d = sqrt(cx * cx + cy * cy + cz * cz);
  cx = r * cx / d;
  cy = r * cy / d;
  cz = r * cz / d;
}
}}
