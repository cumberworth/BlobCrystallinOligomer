// test_potential.cpp

#include "catch2/catch.hpp"

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
        eneT eps {1};
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
        WHEN("Distance is past cutoff") {
            distT rdist {5};
            THEN("Energy is zero") {
                REQUIRE(pot.calc_energy(rdist, pl, ore, ore) == 0);
            }
        }
    }
    GIVEN("Patchy potential") {
        eneT eps {1};
        distT sigl {1};
        distT rcut {4};
        distT siga1 {0.9};
        distT siga2 {1.1};
        PatchyPotential pot {eps, sigl, rcut, siga1, siga2};
        vecT pl {0, 0, 0}; // Placeholder for vectors
        vecT diff {pl};
        distT rdist {0};
        Orientation ore1 {pl, pl};
        Orientation ore2 {pl, pl};
        WHEN("Patches are pointing directly at each other and in attractive regime") {
            diff = {2, 0, 0};
            rdist = 2;
            ore1.patch_norm = {1, 0, 0};
            ore2.patch_norm = {-1, 0, 0};
            eneT c_ene {pot.calc_energy(rdist, diff, ore1, ore2)};
            THEN("Energy reduces to shifted LJ") {
                /* 4*eps*(((sigl/rdist)**12 - (sigl/rdist)**6) -
                   ((sigl/rcut)**12 - (sigl/rcut)**6))
                */
                eneT e_ene {-0.0605471134185791};
                REQUIRE(c_ene == e_ene);
            }
        }
        WHEN("Particles are repulsive regime") {
            distT rdist {0.5};
            eneT c_ene {pot.calc_energy(rdist, diff, ore1, ore2)};
            THEN("Energy reduces to shifted LJ") {
                eneT e_ene {16128.000976324081};
                REQUIRE(c_ene == e_ene);
            }
        }
        WHEN("Distance is past cutoff") {
            distT rdist {5};
            eneT c_ene {pot.calc_energy(rdist, diff, ore1, ore2)};
            THEN("Energy is zero") {
                REQUIRE(c_ene == 0);
            }
        }
        WHEN("Patches are pointing fully away from each other and in attractive regime") {
            diff = {2, 0, 0};
            rdist = 2;
            ore1.patch_norm = {-1, 0, 0};
            ore2.patch_norm = {1, 0, 0};
            eneT c_ene {pot.calc_energy(rdist, diff, ore1, ore2)};
            THEN("Potential is reduced from shifted LJ value") {
                // np.exp(-np.arccos(-1)**2/(2*siga1**2) + -np.arccos(-1)**2/(2*siga2**2)) * -0.0605471134185791
                eneT e_ene {-2.3174785282120287e-06};
                REQUIRE(c_ene == Approx(e_ene));
            }
        }
        WHEN("Patches in random directions and in attractive regime") {
            diff = {2, 0, 0};
            rdist = 2;
            ore1.patch_norm = {-0.22434194136077418, 0.29142133850121443, -0.9299162848410816};
            ore2.patch_norm = {0.1781393764553412, 0.9554100741642645, 0.23549512254298774};
            eneT c_ene {pot.calc_energy(rdist, diff, ore1, ore2)};
            THEN("Potential is reduced from shifted LJ value") {
                // np.exp(-np.arccos((diff/rdist).dot(p1))**2/(2*siga1**2) + -np.arccos((-diff/rdist).dot(p2))**2/(2*siga2**2)) * -0.0605471134185791
                eneT e_ene {-0.0023270508981328374};
                REQUIRE(c_ene == Approx(e_ene));
            }
        }
    }
    GIVEN("Oriented patchy potential") {
        eneT eps {1};
        distT sigl {1};
        distT rcut {4};
        distT siga1 {0.9};
        distT siga2 {1.1};
        distT sigt {1.2};
        OrientedPatchyPotential pot {eps, sigl, rcut, siga1, siga2, sigt};
        vecT pl {0, 0, 0}; // Placeholder for vectors
        vecT diff {pl};
        distT rdist {0};
        Orientation ore1 {pl, pl};
        Orientation ore2 {pl, pl};
        WHEN("Directional are point at each other and in attractive regime") {
            diff = {2, 0, 0};
            rdist = 2;
            ore1.patch_norm = {1, 0, 0};
            ore2.patch_norm = {-1, 0, 0};
            WHEN("Orientational patches are parallel") {
                ore1.patch_orient = {0, 1, 0};
                ore2.patch_orient = {0, 1, 0};
                eneT c_ene {pot.calc_energy(rdist, diff, ore1, ore2)};
                THEN("Energy reduces to shifted LJ") {
                    eneT e_ene {-0.0605471134185791};
                    REQUIRE(c_ene == e_ene);
                }
            }
            WHEN("Are aligned but not parallel") {
                ore1.patch_orient = {-0.9969659883982727, 0.07783840939443323, 0};
                ore2.patch_orient = {-0.3833404063630678, 0.9236071312248503, 0};
                eneT c_ene {pot.calc_energy(rdist, diff, ore1, ore2)};
                THEN("Energy reduces to shifted LJ") {
                    eneT e_ene {-0.0605471134185791};
                    REQUIRE(c_ene == e_ene);
                }
            }
            WHEN("Have antiparallel alignments when projected") {
                ore1.patch_orient = {-0.9969659883982727, 0.07783840939443323, 0};
                ore2.patch_orient = {-0.3833404063630678, -0.9236071312248503, 0};
                eneT c_ene {pot.calc_energy(rdist, diff, ore1, ore2)};
                THEN("Energy reduces to shifted LJ") {
                // np.exp(-np.arccos(-1)**2/(2*sigt**2) * -0.0605471134185791
                    eneT e_ene {-0.0019669336820698456};
                    REQUIRE(c_ene == e_ene);
                }
            }
        }
        WHEN("Particles are repulsive regime") {
            distT rdist {0.5};
            eneT c_ene {pot.calc_energy(rdist, diff, ore1, ore2)};
            THEN("Energy reduces to shifted LJ") {
                eneT e_ene {16128.000976324081};
                REQUIRE(c_ene == e_ene);
            }
        }
        WHEN("Distance is past cutoff") {
            distT rdist {5};
            eneT c_ene {pot.calc_energy(rdist, diff, ore1, ore2)};
            THEN("Energy is zero") {
                REQUIRE(c_ene == 0);
            }
        }
        WHEN("Patches in random directions and in attractive regime") {
            diff = {2, 0, 0};
            rdist = 2;
            ore1.patch_norm = {-0.22434194136077418, 0.29142133850121443, -0.9299162848410816};
            ore1.patch_orient = {-0.25111817, 0.96737965, 0.03341075};
            ore2.patch_norm = {0.1781393764553412, 0.9554100741642645, 0.23549512254298774};
            ore2.patch_orient = {0.92184494, -0.09767422, -0.37504886};
            eneT c_ene {pot.calc_energy(rdist, diff, ore1, ore2)};
            THEN("Potential is reduced from shifted LJ value") {
                // proj1 = op1 - op1.dot(diff/rdist)*diff/rdist
                // proj2 = op2 - op2.dot(-diff/rdist)*-diff/rdist
                // phi = np.arccos(proj1.dot(proj2)/(np.linalg.norm(proj1)*np.linalg.norm(proj2)))
                // np.exp(-np.arccos((diff/rdist).dot(p1))**2/(2*siga1**2) -np.arccos((-diff/rdist).dot(p2))**2/(2*siga2**2) -phi**2/(2*sigt**2)) * -0.0605471134185791
                eneT e_ene {-0.00069993607228337282};
                REQUIRE(c_ene == Approx(e_ene));
            }
        }
    }
}
