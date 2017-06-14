// particle.h

#ifndef PARTICLE_H
#define PARTICLE_H

#include "BlobCrystallinOligomer/shared_types.h"

namespace particle {

    using shared_types::vecT;
    using shared_types::CoorSet;

    struct Orientation {
        vecT patch_norm;
        vecT patch_orient;
    };

    class Particle {
        public:
            Particle(int index, int type, vecT pos, Orientation ore);
            virtual ~Particle();

            int get_index();
            int get_type();
            vecT get_pos(CoorSet coorset);
            Orientation get_ore(CoorSet coorset);

            void translate();
            /*  Translate particle by given vector and store as trial*/

            virtual void rotate();
            /*  Rotate particle by given around given origin and store as trial*/

            void trial_to_current();
            /*  Make trial configuration current configuration */

        protected:
            Orientation m_ore;
            Orientation m_trial_ore;

        private:
            int m_index; // Unique sphere index
            int m_type; // Particle type
            vecT m_pos;
            vecT m_trial_pos;
    };

    class PatchyParticle: public Particle {
        public:
            virtual void rotate();
    };

    class OrientedPatchyParticle: public PatchyParticle {
        public:
            void rotate();
    };
}

#endif // PARTICLE_H
