// test_config.cpp

#include <memory>
#include <vector>

#include "catch2/catch.hpp"

#include "BlobCrystallinOligomer/config.h"
#include "BlobCrystallinOligomer/ifile.h"
#include "BlobCrystallinOligomer/monomer.h"
#include "BlobCrystallinOligomer/particle.h"
#include "BlobCrystallinOligomer/random_gens.h"
#include "BlobCrystallinOligomer/shared_types.h"

SCENARIO("Basic tests of configuration property calculation and access") {
    using config::Config;
    using ifile::ParticleData;
    using ifile::MonomerData;
    using monomer::Monomer;
    using particle::Particle;
    using random_gens::RandomGens;
    using shared_types::vecT;
    using shared_types::CoorSet;
    using shared_types::distT;
    using std::vector;
    using std::unique_ptr;

    //test
    GIVEN("System with two monomers of two simple particles in a box with PBC") {
        RandomGens random_num {};
        distT box_len {10};
        distT radius {1};
        vector<MonomerData> mds;
        int conformer {0};
        CoorSet coorset_current {CoorSet::current};
        for (int i {0}; i != 2; i++) {
            vector<ParticleData> pds;
            for (int j {0}; j != 2; j++) {
                vecT pos {i*2 + j, 0, 0};
                vecT ore {0, 0, 0};
                ParticleData pd {j, "", "SimpleParticle", 0, pos, ore, ore};
                pds.push_back(pd);
            }
            MonomerData md {i, conformer, pds};
            mds.push_back(md);
        }
        Config conf {mds, random_num, box_len, radius};

        // Quick reference
        Monomer& m1 {conf.get_monomer(0)};
        Monomer& m2 {conf.get_monomer(1)};

        // Get particles to test distance calculations on
        Particle& p1 {m1.get_particles()[0].get()};
        Particle& p2 {m2.get_particles()[1].get()};

        WHEN("Monomers are adjacent in middle of box") {
            THEN("Calculated distance is between those copies of the particles") {
                distT e_dist {3};
                distT c_dist {conf.calc_dist(p1, coorset_current, p2, coorset_current)};
                REQUIRE(e_dist == c_dist);
            }
        }
        WHEN("Monomers are moved to opposite sides of box") {
            m1.translate({-4, 0, 0});
            m1.trial_to_current();
            m2.translate({1, 0, 0});
            m2.trial_to_current();
            THEN("Calculated distance follows the minimum image convention") {
                distT e_dist {2};
                distT c_dist {conf.calc_dist(p1, coorset_current, p2, coorset_current)};
                REQUIRE(e_dist == c_dist);
            }
        }
    }
}
