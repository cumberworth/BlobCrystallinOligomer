// monomer.cpp

#include <memory>
#include <vector>
#include <iostream>

#include "BlobCrystallinOligomer/monomer.h"

namespace monomer {

    using particle::Orientation;
    using particle::OrientedPatchyParticle;
    using particle::PatchyParticle;
    using shared_types::CoorSet;
    using shared_types::distT;
    using std::cout;

    Monomer::Monomer(MonomerData m_data, CuboidPBC& pbc_space):
            m_index {m_data.index},
            m_space {pbc_space} {

        create_particles(m_data.particles, pbc_space);

        // Create reference array
        for (auto &p: m_particles) {
            m_particle_refs.emplace_back(*p);
        }

        m_num_particles = m_particles.size();
        calc_monomer_radius();
    }

    int Monomer::get_index() {
        return m_index;
    }

    particleArrayT Monomer::get_particles() {
        return m_particle_refs;
    }

    int Monomer::get_num_particles() {
        return m_num_particles;
    }

    vecT Monomer::get_center(CoorSet coorset) {
        vecT center {0, 0, 0};
        vecT& prev_pos {m_particle_refs[0].get().get_pos(coorset)};
        center += prev_pos;
        for (size_t i {1}; i != m_particle_refs.size(); i++) {
            vecT& pos {m_particle_refs[i].get().get_pos(coorset)};
            center += m_space.unwrap(prev_pos, pos);
        }
        center /= m_particles.size();

        return m_space.wrap(center);
    }

    distT Monomer::get_radius() {
        return m_r;
    }

    void Monomer::translate(vecT disv) {
        for (size_t i {0}; i != m_particles.size(); i++) {
            Particle& particle {m_particle_refs[i].get()};
            particle.translate(disv);
        }
    }

    void Monomer::rotate(vecT rot_c, rotMatT rot_mat) {
        for (size_t i {0}; i != m_particles.size(); i++) {
            Particle& particle {m_particle_refs[i].get()};
            particle.rotate(rot_c, rot_mat);
        }
    }

    void Monomer::current_to_trial() {
        for (size_t i {0}; i != m_particles.size(); i++) {
            Particle& particle {m_particle_refs[i].get()};
            particle.current_to_trial();
        }
    }

    void Monomer::trial_to_current() {
        for (size_t i {0}; i != m_particles.size(); i++) {
            Particle& particle {m_particle_refs[i].get()};
            particle.trial_to_current();
        }
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

    void Monomer::calc_monomer_radius() {
        vecT c {get_center(CoorSet::current)};
        distT max_d {};
        for (Particle& p: m_particle_refs) {
            distT d {m_space.calc_dist(p.get_pos(CoorSet::current), c)};
            if (d > max_d) {
                max_d = d;
            }
        }

        m_r = max_d;
    }
}
