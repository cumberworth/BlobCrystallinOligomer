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
    using shared_types::rotMatT;
    using shared_types::vecT;
    using space::CuboidPBC;
    using std::vector;
    using std::reference_wrapper;
    using std::unique_ptr;

    typedef vector<reference_wrapper<Particle>> particleArrayT;

    // Consider making multiple classes for different versions of the alphaB
    // monomer model, or in the distant future, alphaA monomers
    class Monomer {
        public:
            Monomer(MonomerData m_data, CuboidPBC& pbc_space);

            // Accessors
            int get_index();
            /*  Unique index */

            particleArrayT get_particles();
            /*  Get all particles */

            int get_num_particles();

            vecT get_center();
            /* Get geometric center of all particles */

            // Configuration manipulation
            void translate(vecT disv);
            /*  Translate monomer by given vector */

            void rotate(vecT rot_c, rotMatT rot_mat);
            /*  Rotate monomer by given ? around given origin */

            void trial_to_current();
            /*  Make trial configuration current configuration */

        private:
            int m_index; // Unique monomer index
            vector<unique_ptr<Particle>> m_particles;
            particleArrayT m_particle_refs;
            int m_num_particles;

            void create_particles(vector<ParticleData> p_datas,
                    CuboidPBC& pbc_space);
    };
}

#endif // MONOMER_H
