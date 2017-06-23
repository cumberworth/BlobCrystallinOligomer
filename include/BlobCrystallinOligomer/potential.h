// potential.h

#ifndef POTENTIAL_H
#define POTENTIAL_H

#include "BlobCrystallinOligomer/particle.h"
#include "BlobCrystallinOligomer/shared_types.h"

namespace potential {

    using particle::Orientation;
    using shared_types::distT;
    using shared_types::eneT;
    using shared_types::vecT;

    // Should do a fast version that takes presquared values
    eneT guassian(distT theta, distT sig);

    class PairPotential {
        public:
            PairPotential(distT rcut);
            virtual eneT calc_energy(distT rdist, vecT& p_diff, Orientation& ore1,
                    Orientation& ore2) = 0;
            bool particles_interacting(distT rdist);

        private:
            distT m_rcut;
    };

    class HardSpherePotential: public PairPotential {
        public:
            HardSpherePotential(distT sigh);
            eneT calc_energy(distT rdist, vecT&, Orientation&, Orientation&);
        private:
            distT m_sigh; // Sphere radius
    };

    class ShiftedLJPotential: public PairPotential {
        public:
            ShiftedLJPotential(eneT eps, distT sigl, distT rcut);
            eneT calc_energy(distT rdist, vecT&, Orientation& ore1, Orientation& ore2);

        private:
            eneT m_eps; // Well depth
            eneT m_four_eps; // Well depth premultiplied by 4
            distT m_sigl; // Zero point
            distT m_rcut; // Cutoff
            eneT m_shift; // Shift
    };

    class PatchyPotential: public PairPotential {
        public:
            PatchyPotential(eneT eps, distT sigl, distT rcut, distT siga1,
                    distT siga2);
            eneT calc_energy(distT rdist, vecT& p_diff, Orientation& ore1,
                    Orientation& ore2);

        private:
            ShiftedLJPotential m_lj; // Radial component
            distT m_siga1; // Patch width of particle 1
            distT m_siga2; // Patch width of particle 2
    };

    class OrientedPatchyPotential: public PairPotential {
        public:
            OrientedPatchyPotential(eneT eps, distT sigl, distT rcut,
                    distT siga1, distT siga2, distT sigt);
            eneT calc_energy(distT rdist, vecT& p_diff, Orientation& ore1,
                    Orientation& ore2);

        private:
            PatchyPotential m_patchy; // Unoriented-patchy potential
            distT m_sigt; // Orientation width
    };
}

#endif // POTENTIAL_H
