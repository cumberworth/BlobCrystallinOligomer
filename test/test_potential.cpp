// test_potential.cpp

#include "catch/catch.hpp"

#include "BlobCrystallinOligomer/particle.h"
#include "BlobCrystallinOligomer/potential.h"
#include "BlobCrystallinOligomer/shared_types.h"

SCENARIO("Compare potential values to hand calculated values") {
    using particle::Orientation;
    using potential::HardSpherePotential;
    using potential::ShiftedLJPotential;
    using potential::PatchyPotential;
    using potential::OrientedPatchyPotential;
    using shared_types::distT;
    using shared_types::eneT;
    using shared_types::inf;
    using shared_types::vecT;

    GIVEN("Hard sphere potential") {
        distT sigh {1};
        HardSpherePotential pot {sigh};
        vecT pl {0, 0, 0}; // Placeholder for vectors
        Orientation ore {pl, pl}; // Placeholder for orientations
        WHEN("Distance is greater than sphere overlap") {
            distT rdist {sigh*2};
            THEN("Energy is 0") {
                eneT c_ene {pot.calc_energy(rdist, pl, ore, ore)};
                REQUIRE(c_ene == 0);
            }
        }
        WHEN("Distance is less than sphere overlap") {
            distT rdist {sigh/2};
            THEN("Energy is infinity") {
                eneT c_ene {pot.calc_energy(rdist, pl, ore, ore)};
                REQUIRE(c_ene == inf);
            }
        }
    }

    GIVEN("Shifted Lennard Jones potential") {
        distT eps {1};
        distT sigl {1};
        distT rcut {4};
        ShiftedLJPotential pot {eps, sigl, rcut};
        vecT pl {0, 0, 0}; // Placeholder for vectors
        Orientation ore {pl, pl}; // Placeholder for orientations
        WHEN("Distance is in attractive regime") {
            distT rdist {2};
            THEN("Energy is negative (and agrees with calculated value)") {
                eneT c_ene {pot.calc_energy(rdist, pl, ore, ore)};
                /* 4*eps*(((sigl/rdist)**12 - (sigl/rdist)**6) -
                   ((sigl/rcut)**12 - (sigl/rcut)**6))
                */
                eneT e_ene {-0.0605471134185791};
                REQUIRE(c_ene == e_ene);
            }
        }
        WHEN("Distance is in repulsive regime") {
            distT rdist {0.5};
            THEN("Energy is positive (and agrees with calculated value)") {
                eneT c_ene {pot.calc_energy(rdist, pl, ore, ore)};
                /* 4*eps*(((sigl/rdist)**12 - (sigl/rdist)**6) -
                   ((sigl/rcut)**12 - (sigl/rcut)**6))
                */
                eneT e_ene {16128.000976324081};
                REQUIRE(c_ene == e_ene);
            }
        }
    }
}
