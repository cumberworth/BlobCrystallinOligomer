// particle.h

#ifndef PARTICLE_H
#define PARTICLE_H

#include "shared_types.h"

namespace particle {

    using shared_types::vecT;

    class Particle {
        public:
            int get_index();
            int get_type();
            vecT get_pos();

            void translate();
            // Translate particle by given vector

            void rotate();
            // Rotate particle by given ? around given origin

        private:
            int m_index;
            int m_type;
            vecT m_pos;
    };

    class PatchyParticle: public Particle {
        public:
            vecT get_patch_norm();

        private:
            vecT m_patch_norm;

    };

    class OrientedPatchyParticle: public PatchyParticle {
        public:
            vecT get_patch_orient();

        private:
            vecT m_patch_orient;
}

#endif // PARTICLE_H
