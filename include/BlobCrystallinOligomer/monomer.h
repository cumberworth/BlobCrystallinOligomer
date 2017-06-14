// monomer.h

#ifndef MONOMER_H
#define MONOMER_H

#include <memory>
#include <vector>

#include "BlobCrystallinOligomer/file.h"
#include "BlobCrystallinOligomer/particle.h"
#include "BlobCrystallinOligomer/shared_types.h"

namespace monomer {

    using file::MonomerData;
    using file::ParticleData;
    using std::vector;
    using std::reference_wrapper;
    using std::unique_ptr;
    using shared_types::vecT;
    using particle::Particle;

    typedef vector<reference_wrapper<Particle>> particleArrayT;

    // Consider making multiple classes for different versions of the alphaB
    // monomer model, or in the distant future, alphaA monomers
    class Monomer {
        public:
            Monomer(MonomerData m_data);


            // Accessors
            particleArrayT get_particles();
            /*  Get all particles */

            // Configuration manipulation
            void translate(vecT disV);
            /*  Translate monomer by given vector */

            void rotate();
            /*  Rotate monomer by given ? around given origin */

        private:
            int m_index; // Unique monomer index
            vector<unique_ptr<Particle>> m_particles;
            particleArrayT m_particle_refs;

            void create_particles(vector<ParticleData> p_datas);
    };
}

#endif // MONOMER_H
