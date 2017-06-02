// monomer.h

#ifndef MONOMER_H
#define MONOMER_H

#include <vector>

#include "BlobCrystallinOligomer/shared_types.h"
#include "BlobCrystallinOligomer/particle.h"

namespace monomer {

    using std::vector;
    using std::reference_wrapper;
    using shared_types::vecT;
    using particle::Particle;

    typedef vector<reference_wrapper<Particle>> particleArrayT;

    // Consider making multiple classes for different versions of the alphaB
    // monomer model, or in the distant future, alphaA monomers

    class Monomer {
        public:
            Monomer();
            void translate(vecT disV);
            // Translate monomer by given vector

            void rotate();
            // Rotate monomer by given ? around given origin

        private:
            int m_index; // Unique monomer index
            particleArrayT m_particles;

    };
}

#endif // MONOMER_H
