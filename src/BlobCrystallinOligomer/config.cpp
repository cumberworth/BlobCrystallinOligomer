// config.cpp

#include <iostream>
#include <memory>

#include "BlobCrystallinOligomer/config.h"
#include "BlobCrystallinOligomer/file.h"
#include "BlobCrystallinOligomer/monomer.h"
#include "BlobCrystallinOligomer/shared_types.h"

namespace config {

    using file::InputConfigFile;
    using file::MonomerData;
    using file::ParticleData;
    using monomer::Monomer;
    using particle::Particle;
    using shared_types::vecT;
    using std::cout;
    using std::unique_ptr;

    Config::Config(InputParams params):
            m_space {space::CuboidPBC()} {

        InputConfigFile config_file {params.m_config_filename};
        vector<MonomerData> monomers {config_file.get_monomers()};
        create_monomers(monomers);

        // Create monomer reference array
        for (auto &m: m_monomers) {
            m_monomer_refs.emplace_back(m);
        }
    }

    Monomer& Config::get_monomer(int monomer_index) {
        return *m_monomers[monomer_index];
    }
    
    monomerArrayT Config::get_monomers() {
        return m_monomer_refs;
    }

    vecT Config::calc_interparticle_vector(Particle& particle1, CoorSet coorset1,
            Particle& particle2, CoorSet coorset2) {

        vecT pos1 {particle1.get_pos(coorset1)};
        vecT pos2 {particle2.get_pos(coorset2)};

        return m_space.calc_diff(pos1, pos2);
    }

    distT Config::calc_dist(Particle& particle1, CoorSet coorset1,
            Particle& particle2, CoorSet coorset2) {

        vecT pos1 {particle1.get_pos(coorset1)};
        vecT pos2 {particle2.get_pos(coorset2)};

        return m_space.calc_dist(pos1, pos2);
    }

    void Config::create_monomers(vector<MonomerData> monomers) {
        for (auto m_data: monomers) {
            m_monomers.emplace_back(m_data);
        }
    }
}
