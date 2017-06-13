// particle.h

#ifndef PARTICLE_H
#define PARTICLE_H

#include "BlobCrystallinOligomer/shared_types.h"

namespace particle {

    using shared_types::vecT;
    using shared_types::CoorSet;

    class Particle {
        public:
            Particle(int index, int type, vecT pos);
            virtual ~Particle();
            int get_index();
            int get_type();
            vecT get_pos(CoorSet coorset);

            void translate();
            // Translate particle by given vector

            void rotate();
            // Rotate particle by given ? around given origin

        private:
            int m_index; // Unique sphere index
            int m_type; // Particle type
            vecT m_pos;
    };

    class PatchyParticle: public Particle {
        public:
            PatchyParticle(int index, int type, vecT pos, vecT patch_norm);
            vecT get_patch_norm();

        private:
            vecT m_patch_norm;

    };

    class OrientedPatchyParticle: public PatchyParticle {
        public:
            OrientedPatchyParticle(int index, int type, vecT pos,
                    vecT patch_norm, vecT patch_orient);
            vecT get_patch_orient();

        private:
            vecT m_patch_orient;
    };
}

#endif // PARTICLE_H
