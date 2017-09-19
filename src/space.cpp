// space.cpp

#include "BlobCrystallinOligomer/space.h"
#include "BlobCrystallinOligomer/shared_types.h"

namespace space {

    using shared_types::distT;
    using shared_types::vecT;

    CuboidPBC::CuboidPBC() {}

    CuboidPBC::CuboidPBC(distT len): m_r {len/2}  {}

    void CuboidPBC::set_len(distT len) {
        m_r = len/2;
    }

    distT CuboidPBC::calc_dist(vecT pos1, vecT pos2) {
        vecT diff {calc_diff(pos1, pos2)};

        return diff.norm();
    }

    vecT CuboidPBC::calc_diff(vecT pos1, vecT pos2) {
        vecT diff;
        for (int i {0}; i != 3; i++) {
            distT comp_diff {pos1[i] - pos2[i]};
            if (comp_diff > m_r) {
                comp_diff = -2*m_r + comp_diff;
            }
            else if (comp_diff < -m_r) {
                comp_diff = 2*m_r + comp_diff;
            }
            diff[i] = comp_diff;
        }

        return diff;
    }

    vecT CuboidPBC::wrap(vecT pos) {
        for (int i {0}; i != 3; i++) {
            if (pos[i] > m_r) {
                pos[i] = -2*m_r + pos[i];
            }
            else if (pos[i] < -m_r) {
                pos[i] = 2*m_r + pos[i];
            }
        }

        return pos;
    }
}
