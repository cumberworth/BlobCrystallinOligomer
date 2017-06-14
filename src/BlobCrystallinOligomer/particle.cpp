// particle.cpp

#include "BlobCrystallinOligomer/particle.h"
#include "BlobCrystallinOligomer/shared_types.h"

namespace particle {

    using shared_types::CoorSet;

    Particle::Particle(int index, int type, vecT pos):
            m_index {index}, m_type {type}, m_pos {pos} {
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

    void Particle::trial_to_current() {
        m_pos = m_trial_pos;
        m_ore = m_trial_ore;
    }
}
