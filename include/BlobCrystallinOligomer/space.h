// space.h 

#ifndef SPACE_H
#define SPACE_H

#include "BlobCrystallinOligomer/shared_types.h"

namespace space {

    using shared_types::distT;
    using shared_types::vecT;

    class CuboidPBC {
        public:
            CuboidPBC(); // Need to take box size as args

            distT calc_dist(vecT pos1, vecT pos2);
            vecT calc_diff(vecT pos1, vecT pos2);
    };
}

#endif // SPACE_H
