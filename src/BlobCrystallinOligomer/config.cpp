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
    using monomer::particleArrayT;
    using particle::OrientedPatchyParticle;
    using particle::Particle;
    using particle::PatchyParticle;
    using shared_types::vecT;
    using std::cout;
    using std::shared_ptr;

    Config::Config(InputParams params):
            m_space {space::CuboidPBC()} {

        InputConfigFile config_file {params.m_config_filename};
        vector<MonomerData> monomers {config_file.get_monomers()};
        extract_config_from_monomers(monomers);
    }

    void Config::extract_config_from_monomers(vector<MonomerData> monomers) {
        for (auto m_data: monomers) {
            particleArrayT particles {};
            for (auto p_data: m_data.particles) {
                int type; // TODO: type string to type int table
                // TODO; also need to extract and set monomer state

                // Particle types specified as string
                shared_ptr<Particle> part;
                if (p_data.type == "SimpleParticle") {
                    part = std::make_shared<Particle>(
                            new Particle {p_data.index, type, p_data.pos});
                }
                else if (p_data.type == "PatchyParticle") {
                    part = std::make_shared<Particle>(
                            new PatchyParticle {p_data.index, type,
                            p_data.pos, p_data.patch_norm});
                }
                else if (p_data.type == "OrientedPatchyParticle") {
                    part = std::make_shared<Particle>(
                            new OrientedPatchyParticle {p_data.index, type,
                            p_data.pos, p_data.patch_norm, p_data.patch_orient});
                }
                else {
                    cout << "Particle type unknown\n";
                    throw shared_types::InputError {};
                }
                particles.push_back(*part);
            }
            shared_ptr<Monomer> monomer;
            int m_index {m_data.index};
            monomer = std::make_shared<Monomer>(new Monomer {m_index,
                    particles});
        }
    }

    Monomer& Config::get_monomer(int monomer_index) {
        return m_monomers[monomer_index];
    }

    distT Config::calc_dist(Particle& particle1, CoorSet coorset1,
            Particle& particle2, CoorSet coorset2) {

        vecT pos1 {particle1.get_pos(coorset1)};
        vecT pos2 {particle2.get_pos(coorset2)};
        m_space.calc_dist(pos1, pos2);
    }
}
