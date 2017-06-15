// particle.cpp

#include "BlobCrystallinOligomer/particle.h"
#include "BlobCrystallinOligomer/shared_types.h"
#include "BlobCrystallinOligomer/space.h"

namespace particle {

    using shared_types::CoorSet;
    using space::CuboidPBC;

    Particle::Particle(int index, int type, vecT pos, Orientation ore,
            CuboidPBC& pbc_space):
            m_ore {ore}, m_index {index}, m_type {type}, m_pos {pos},
            m_space {pbc_space} {
    }

    int Particle::get_index() {
        return m_index;
    }

    int Particle::get_type() {
        return m_type;
    }

    vecT Particle::get_pos(CoorSet coorset) {
        if (coorset == CoorSet::current) {
            return m_pos;
        }
        else {
            return m_trial_pos;
        }
    }

    Orientation Particle::get_ore(CoorSet coorset) {
        if (coorset == CoorSet::current) {
            return m_ore;
        }
        else {
            return m_trial_ore;
        }
    }

    void Particle::translate(vecT disv) {
        vecT new_pos {m_pos + disv};
        new_pos = m_space.wrap(new_pos);
        m_trial_pos = new_pos;
    }

    void Particle::trial_to_current() {
        m_pos = m_trial_pos;
        m_ore = m_trial_ore;
    }

    PatchyParticle::PatchyParticle(int index, int type, vecT pos,
            Orientation ore, CuboidPBC& pbc_space):
            Particle {index, type, pos, ore, pbc_space} {
    }

    OrientedPatchyParticle::OrientedPatchyParticle(int index, int type,
            vecT pos, Orientation ore, CuboidPBC& pbc_space):
            PatchyParticle {index, type, pos, ore, pbc_space} {
    }

}
