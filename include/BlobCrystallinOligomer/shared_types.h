// shared_types.h

#ifndef SHARED_TYPES_H
#define SHARED_TYPES_H

#include <Eigen/Dense>

namespace shared_types {

    typedef Eigen::Vector3d vecT;
    typedef Eigen::Matrix3d rotMatT;
    typedef double distT;
    typedef double eneT;
    typedef unsigned long long int stepT;

    enum class CoorSet {
        current,
        trial
    };

    // Exception structs
    struct InputError {};

}


#endif // SHARED_TYPES_H
