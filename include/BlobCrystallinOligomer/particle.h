// particle.h

#ifndef PARTICLE_H
#define PARTICLE_H

#include "BlobCrystallinOligomer/shared_types.h"
#include "BlobCrystallinOligomer/space.h"

namespace particle {

    using shared_types::vecT;
    using shared_types::CoorSet;
    using space::CuboidPBC;

    struct Orientation {
        vecT patch_norm;
        vecT patch_orient;
    };

    class Particle {
        public:
            Particle(int index, int type, vecT pos, Orientation ore,
                    CuboidPBC& pbc_space);
            virtual ~Particle();

            int get_index();
            int get_type();
            vecT get_pos(CoorSet coorset);
            Orientation get_ore(CoorSet coorset);

            void translate(vecT disv);
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
            CuboidPBC& m_space;
    };

    class PatchyParticle: public Particle {
        public:
            PatchyParticle(int index, int type, vecT pos, Orientation ore,
                    CuboidPBC& pbc_space);
            virtual void rotate();
    };

    class OrientedPatchyParticle: public PatchyParticle {
        public:
            OrientedPatchyParticle(int index, int type, vecT pos,
                    Orientation ore, CuboidPBC& pbc_space);
            void rotate();
    };
}

#endif // PARTICLE_H
