// particle.cpp

#include <iostream>

#include "BlobCrystallinOligomer/particle.h"

namespace particle {

    using std::cout;

    Particle::Particle(int index, int type, vecT pos, Orientation ore,
            CuboidPBC& pbc_space):
            m_ore {ore},
            m_trial_ore {ore},
            m_index {index},
            m_type {type},
            m_pos {pos},
            m_trial_pos {pos},
            m_space {pbc_space} {
    }

    int Particle::get_index() {
        return m_index;
    }

    int Particle::get_type() {
        return m_type;
    }

    vecT& Particle::get_pos(CoorSet coorset) {
        if (coorset == CoorSet::current) {
            return m_pos;
        }
        else {
            return m_trial_pos;
        }
    }

    Orientation& Particle::get_ore(CoorSet coorset) {
        if (coorset == CoorSet::current) {
            return m_ore;
        }
        else {
            return m_trial_ore;
        }
    }

    void Particle::set_pos(vecT pos) {
        m_pos = pos;
    }

    void Particle::translate(vecT& disv) {
        m_trial_pos = m_pos + disv;
        m_trial_pos = m_space.wrap(m_trial_pos);
    }

    void Particle::rotate(vecT& rot_c, rotMatT& rot_mat) {
        m_trial_pos -= rot_c;
        m_trial_pos = rot_mat * m_trial_pos;
        m_trial_pos += rot_c;
        m_trial_pos = m_space.wrap(m_trial_pos);
    }

    void Particle::trial_to_current() {
        m_pos = m_trial_pos;
        m_ore = m_trial_ore;
    }

    void Particle::current_to_trial() {
        m_trial_pos = m_pos;
        m_trial_ore = m_ore;
    }

    PatchyParticle::PatchyParticle(int index, int type, vecT pos,
            Orientation ore, CuboidPBC& pbc_space):
            Particle {index, type, pos, ore, pbc_space} {
    }

    void PatchyParticle::rotate(vecT& rot_c, rotMatT& rot_mat) {
        Particle::rotate(rot_c, rot_mat);
        m_trial_ore.patch_norm = rot_mat * m_ore.patch_norm;
    }

    OrientedPatchyParticle::OrientedPatchyParticle(int index, int type,
            vecT pos, Orientation ore, CuboidPBC& pbc_space):
            PatchyParticle {index, type, pos, ore, pbc_space} {
    }

    void OrientedPatchyParticle::rotate(vecT& rot_c, rotMatT& rot_mat) {
        PatchyParticle::rotate(rot_c, rot_mat);
        m_trial_ore.patch_orient = rot_mat * m_ore.patch_orient;
    }

    DoubleOrientedPatchyParticle::DoubleOrientedPatchyParticle(int index,
            int type, vecT pos, Orientation ore, CuboidPBC& pbc_space):
            PatchyParticle {index, type, pos, ore, pbc_space} {
    }

    void DoubleOrientedPatchyParticle::rotate(vecT& rot_c, rotMatT& rot_mat) {
        PatchyParticle::rotate(rot_c, rot_mat);
        m_trial_ore.patch_orient = rot_mat * m_ore.patch_orient;
        m_trial_ore.patch_orient2 = rot_mat * m_ore.patch_orient2;
    }
}
