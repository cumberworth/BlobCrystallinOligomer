// monomer.cpp

#include <memory>
#include <vector>
#include <iostream>

#include "BlobCrystallinOligomer/file.h"
#include "BlobCrystallinOligomer/monomer.h"
#include "BlobCrystallinOligomer/particle.h"

namespace monomer {

    using file::MonomerData;
    using file::ParticleData;
    using particle::OrientedPatchyParticle;
    using particle::Particle;
    using particle::PatchyParticle;
    using std::cout;
    using std::unique_ptr;
    using std::vector;

    Monomer::Monomer(MonomerData m_data):
            m_index {m_data.index} {

        create_particles(m_data.particles);

        // Create reference array
        for (auto &p: m_particles) {
            m_particle_refs.emplace_back(p);
        }
    }

    particleArrayT Monomer::get_particles() {

        return m_particle_refs;
    }

    void Monomer::create_particles(vector<ParticleData> p_datas) {
        for (auto p_data: p_datas) {
            int type {p_data.type};
            Particle* part;
            if (p_data.form == "SimpleParticle") {
                part = new Particle {p_data.index, type,
                        p_data.pos};
            }
            else if (p_data.form == "PatchyParticle") {
                part = new PatchyParticle {p_data.index, type,
                        p_data.pos, p_data.patch_norm};
            }
            else if (p_data.form == "OrientedPatchyParticle") {
                part = new OrientedPatchyParticle {p_data.index, type,
                        p_data.pos, p_data.patch_norm, p_data.patch_orient};
            }
            else {
                cout << "Particle type unknown\n";
                throw shared_types::InputError {};
            }
            m_particles.emplace_back(part);
        }
    }
}
