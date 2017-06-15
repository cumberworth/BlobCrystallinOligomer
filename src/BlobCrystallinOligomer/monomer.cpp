// monomer.cpp

#include <memory>
#include <vector>
#include <iostream>

#include "BlobCrystallinOligomer/file.h"
#include "BlobCrystallinOligomer/monomer.h"
#include "BlobCrystallinOligomer/particle.h"
#include "BlobCrystallinOligomer/space.h"

namespace monomer {

    using file::MonomerData;
    using file::ParticleData;
    using particle::Orientation;
    using particle::OrientedPatchyParticle;
    using particle::Particle;
    using particle::PatchyParticle;
    using space::CuboidPBC;
    using std::cout;
    using std::unique_ptr;
    using std::vector;

    Monomer::Monomer(MonomerData m_data, CuboidPBC& pbc_space):
            m_index {m_data.index} {

        create_particles(m_data.particles, pbc_space);

        // Create reference array
        for (auto &p: m_particles) {
            m_particle_refs.emplace_back(p);
        }
    }

    int Monomer::get_index() {
        return m_index;
    }

    particleArrayT Monomer::get_particles() {

        return m_particle_refs;
    }

    void Monomer::create_particles(vector<ParticleData> p_datas,
            CuboidPBC& pbc_space) {
        for (auto p_data: p_datas) {
            int type {p_data.type};
            Particle* part;
            Orientation ore {};
            if (p_data.form == "SimpleParticle") {
                part = new Particle {p_data.index, type,
                        p_data.pos, ore, pbc_space};
            }
            else if (p_data.form == "PatchyParticle") {
                ore.patch_norm = p_data.patch_norm;
                part = new PatchyParticle {p_data.index, type,
                        p_data.pos, ore, pbc_space};
            }
            else if (p_data.form == "OrientedPatchyParticle") {
                ore.patch_norm = p_data.patch_norm;
                ore.patch_orient = p_data.patch_orient;
                part = new OrientedPatchyParticle {p_data.index, type,
                        p_data.pos, ore, pbc_space};
            }
            else {
                cout << "Particle type unknown\n";
                throw shared_types::InputError {};
            }
            m_particles.emplace_back(part);
        }
    }

    void Monomer::translate(vecT disv) {
        for (size_t i {0}; i != m_particles.size(); i++) {
            Particle& particle {m_particle_refs[i].get()};
            particle.translate(disv);
        }
    }

    void Monomer::rotate(?) {
        for (size_t i {0}; i != m_particles.size(); i++) {
            Particle& particle {m_particle_refs[i].get()};
            particle.rotate(?);
        }
    }

    void Monomer::trial_to_current() {
        for (size_t i {0}; i != m_particles.size(); i++) {
            Particle& particle {m_particle_refs[i].get()};
            particle.trial_to_current();
        }
    }
}
