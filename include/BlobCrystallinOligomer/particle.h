// particle.h

#ifndef PARTICLE_H
#define PARTICLE_H

#include "BlobCrystallinOligomer/shared_types.h"
#include "BlobCrystallinOligomer/space.h"

namespace particle {

    using shared_types::rotMatT;
    using shared_types::vecT;
    using shared_types::CoorSet;
    using space::CuboidPBC;

    /** Particle orientation
      *
      * Explicity stores the current configuration of the patch vectors
      */
    struct Orientation {
        vecT patch_norm {0, 0, 0};
        vecT patch_orient {0, 0, 0};
    };

    /** General class for particles, the basic building blocks of structures
      *
      * Contains position vector and type (for deciding interaction potentials.
      * There is both a current position and a trial position for easily reverting
      * changes to the configurations if moves are rejected. Also provides shared interface
      * and implementation for more complex particle types with patches.
      */
    class Particle {
        public:
            Particle(int index, int type, vecT pos, Orientation ore,
                    CuboidPBC& pbc_space);
            virtual ~Particle() {}

            int get_index();
            int get_type();
            vecT& get_pos(CoorSet coorset);
            Orientation& get_ore(CoorSet coorset);

            /** Set the current position */
            void set_pos(vecT& pos);

            /** Translate particle by given vector and store as trial*/
            void translate(vecT& disv);

            /** Rotate particle around given origin and store as trial */
            virtual void rotate(vecT& rot_c, rotMatT& rot_mat);

            /** Make trial configuration current configuration */
            void trial_to_current();

            /** Reset trial with current */
            void current_to_trial();

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

    /** Particle class with one directional patch */
    class PatchyParticle: public Particle {
        public:
            PatchyParticle(int index, int type, vecT pos, Orientation ore,
                    CuboidPBC& pbc_space);
            virtual void rotate(vecT& rot_c, rotMatT& rot_mat);
    };

    /** Particle class with one directional patch and one orientational patch */
    class OrientedPatchyParticle: public PatchyParticle {
        public:
            OrientedPatchyParticle(int index, int type, vecT pos,
                    Orientation ore, CuboidPBC& pbc_space);
            void rotate(vecT& rot_c, rotMatT& rot_mat);
    };
}

#endif // PARTICLE_H
