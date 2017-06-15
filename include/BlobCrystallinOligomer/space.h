// space.h 

#ifndef SPACE_H
#define SPACE_H

#include "BlobCrystallinOligomer/shared_types.h"

namespace space {

    using shared_types::distT;
    using shared_types::vecT;

    class CuboidPBC {
        public:
            CuboidPBC();
            CuboidPBC(distT len);

            void set_len(distT len);
            distT calc_dist(vecT pos1, vecT pos2);
            vecT calc_diff(vecT pos1, vecT pos2);
            vecT wrap(vecT pos);

        private:
            distT m_r;
    };
}

#endif // SPACE_H
