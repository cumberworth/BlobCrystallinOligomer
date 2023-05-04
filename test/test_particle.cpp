// test_particle.cpp

#include <memory>

#include "catch2/catch.hpp"

#include "BlobCrystallinOligomer/particle.h"
#include "BlobCrystallinOligomer/shared_types.h"
#include "BlobCrystallinOligomer/space.h"

#define private public
#define protected public

SCENARIO("Individual particles are moved in a box with PBC") {
    using particle::Orientation;
    using particle::Particle;
    using particle::PatchyParticle;
    using particle::OrientedPatchyParticle;
    using shared_types::vecT;
    using shared_types::CoorSet;
    using shared_types::rotMatT;
    using space::CuboidPBC;
    using std::unique_ptr;

    auto pbc_space {std::make_unique<CuboidPBC>(10)};

    int index {0};
    int type {0};
    vecT s_pos {0, 0, 0};
    vecT s_patch_norm {1, 0, 0};
    vecT s_patch_orient {0, 1, 0};
    Orientation s_ore {s_patch_norm, s_patch_orient};
    
    // Test for oriented particle will also check patchy and simple particle
    GIVEN("One oriented patchy particle") {
        OrientedPatchyParticle part {index, type, s_pos, s_ore, *pbc_space};

        WHEN("Translated by a given displacement vector") {
            vecT d_pos {5, 4, 3};
            part.translate(d_pos);
            THEN("Trial position is updated but not current") {
                REQUIRE(part.get_pos(CoorSet::current) == s_pos);
                REQUIRE(part.get_pos(CoorSet::trial) == d_pos);
            }

            WHEN("Current coordinates are updated") {
                part.trial_to_current();
                THEN("Trial position is updated but not current") {
                    REQUIRE(part.get_pos(CoorSet::current) == d_pos);
                }
            }
        }

        WHEN("Translated outside the periodic box") {
            vecT d_pos {6, -6, 10};
            part.translate(d_pos);
            THEN("Trial position is wrapped back to central box") {
                vecT e_pos {-4, 4, 0};
                REQUIRE(part.get_pos(CoorSet::trial) == e_pos);
            }
        }

        WHEN("Rotated about a given point") {
            vecT crot {1, 1, 1};
            rotMatT rmat;

            // This is positive pi/2 rotation in the z axis
            rmat << 0, -1, 0,
                    1,  0, 0,
                    0,  0, 1;
            part.rotate(crot, rmat);
            THEN("Trial position and orientation updated") {
                vecT e_pos {2, 0, 0};
                REQUIRE(part.get_pos(CoorSet::trial) == e_pos);
                REQUIRE(part.get_pos(CoorSet::current) == s_pos);
                vecT r_patch_norm {part.get_ore(CoorSet::trial).patch_norm};
                vecT c_patch_norm {part.get_ore(CoorSet::current).patch_norm};
                vecT r_patch_orient {part.get_ore(CoorSet::trial).patch_orient};
                vecT c_patch_orient {part.get_ore(CoorSet::current).patch_orient};
                vecT e_patch_norm {0, 1, 0};
                REQUIRE(r_patch_norm == e_patch_norm);
                REQUIRE(c_patch_norm == s_patch_norm);
                vecT e_patch_orient {-1, 0, 0};
                REQUIRE(r_patch_orient == e_patch_orient);
                REQUIRE(c_patch_orient == s_patch_orient);
            }
        }
        WHEN("Rotated to a point outside") {
            s_pos = {-4, 0, 0};
            vecT crot {4, 0, 0};
            rotMatT rmat;

            // This is positive pi/2 rotation in the z axis
            rmat << 0, -1, 0,
                    1,  0, 0,
                    0,  0, 1;
            part.set_pos(s_pos);
            part.current_to_trial();
            part.rotate(crot, rmat);
            THEN("Trial position and orientation updated") {
                vecT e_pos {4, 2, 0};
                REQUIRE(part.get_pos(CoorSet::trial) == e_pos);
                REQUIRE(part.get_pos(CoorSet::current) == s_pos);
                vecT r_patch_norm {part.get_ore(CoorSet::trial).patch_norm};
                vecT c_patch_norm {part.get_ore(CoorSet::current).patch_norm};
                vecT r_patch_orient {part.get_ore(CoorSet::trial).patch_orient};
                vecT c_patch_orient {part.get_ore(CoorSet::current).patch_orient};
                vecT e_patch_norm {0, 1, 0};
                REQUIRE(r_patch_norm == e_patch_norm);
                REQUIRE(c_patch_norm == s_patch_norm);
                vecT e_patch_orient {-1, 0, 0};
                REQUIRE(r_patch_orient == e_patch_orient);
                REQUIRE(c_patch_orient == s_patch_orient);
            }
        }
    }
}
