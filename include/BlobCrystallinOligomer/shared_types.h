// shared_types.h

#ifndef SHARED_TYPES_H
#define SHARED_TYPES_H

#include <limits>

#include "Eigen/Dense"

namespace shared_types {

typedef Eigen::Vector3d vecT;
typedef Eigen::Matrix3d rotMatT;
typedef double distT;
typedef double eneT;
typedef unsigned long long int stepT;
typedef double timeT;
const double inf {std::numeric_limits<double>::infinity()};

/** For specifying which particle coordinate set to use: current or trial */
enum class CoorSet { current, trial };

// Exception structs
struct InputError {};

} // namespace shared_types

#endif // SHARED_TYPES_H
