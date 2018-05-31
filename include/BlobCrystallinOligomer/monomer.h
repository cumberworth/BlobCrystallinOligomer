// monomer.h

#ifndef MONOMER_H
#define MONOMER_H

#include <memory>
#include <vector>

#include "BlobCrystallinOligomer/ifile.h"
#include "BlobCrystallinOligomer/particle.h"
#include "BlobCrystallinOligomer/shared_types.h"
#include "BlobCrystallinOligomer/space.h"

namespace monomer {

    using ifile::MonomerData;
    using ifile::ParticleData;
    using particle::Particle;
    using shared_types::CoorSet;
    using shared_types::distT;
    using shared_types::rotMatT;
    using shared_types::vecT;
    using space::CuboidPBC;
    using std::vector;
    using std::reference_wrapper;
    using std::unique_ptr;

    typedef vector<reference_wrapper<Particle>> particleArrayT;

    // Consider making multiple classes for different versions of the alphaB
    // monomer model, or in the distant future, alphaA monomers

    /** alphB cyrstallin coarse grained monomer
      *
      * Contains the particles that make it up and an interface for manipulating
      * the configuration.
      */
    class Monomer {
        public:
            Monomer(MonomerData m_data, CuboidPBC& pbc_space);

            /** Unique index */
            int get_index();

            /** Conformer (two NTD configs) */
            int get_conformer(CoorSet coorset);

            /** Get specified particle */
            Particle& get_particle(int particle_i);

            /** Get all particles */
            particleArrayT get_particles();

            int get_num_particles();

            /** Get geometric center of all particles */
            vecT get_center(CoorSet coorset);

            /** Get maximum length from monomer center to particle center */
            distT get_radius();

            /** Translate monomer by given vector */
            void translate(vecT disv);

            /** Rotate monomer by given ? around given origin */
            void rotate(vecT rot_c, rotMatT rot_mat);

            /** Unwrap monomer relative to reference position */
            void unwrap(vecT ref_pos);

            /** Flip conformation */
            void flip_conformation();

            /** Make trial configuration current configuration */
            void trial_to_current();

            /** Reset trial to current */
            void current_to_trial();

        private:
            int m_index; // Unique monomer index
            int m_trial_conformer;
            int m_conformer;
            CuboidPBC& m_space;

            vector<unique_ptr<Particle>> m_particles;
            particleArrayT m_particle_refs;
            int m_num_particles;
            distT m_r;

            void create_particles(vector<ParticleData> p_datas,
                    CuboidPBC& pbc_space);
            void calc_monomer_radius();
    };
}

#endif // MONOMER_H
