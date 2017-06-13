// potential.h

#ifndef POTENTIAL_H
#define POTENTIAL_H

#include "BlobCrystallinOligomer/shared_types.h"

namespace potential {

    using shared_types::distT;
    using shared_types::eneT;
    using shared_types::vecT;

    class PairPotential {
        public:
            virtual bool particles_interacting(distT rdist) = 0;
    };

    class DistPairPotential {
        public:
            virtual eneT calc_energy(distT rdist) = 0; // What about orientation?
    };

    class HardSpherePotential: public DistPairPotential {
        public:
            HardSpherePotential(double sigh);
            bool particles_interacting(distT rdist);
        private:
            double m_sigh; // Sphere radius
    };

    class ShiftedLJPotential: public DistPairPotential {
        public:
            ShiftedLJPotential(double eps, double sigl, double rcut);
            eneT calc_energy(distT rdist);
            bool particles_interacting(distT rdist); 

        private:
            double m_eps; // Well depth
            double m_sigl; // Zero point
            double m_rcut; // Cutoff
    };

    class PatchyPotential: public PairPotential {
        public:
            PatchyPotential(double eps, double sigl, double rcut, double siga1,
                    double siga2);
            eneT calc_energy(distT rdist, vecT patch_norm1, vecT patch_norm2);
            bool particles_interacting(distT rdist);

        private:
            ShiftedLJPotential m_lj; // Radial component
            double m_siga1; // Patch width of particle 1
            double m_siga2; // Patch width of particle 2
    };

    class OrientedPatchyPotential: public PairPotential {
        public:
            OrientedPatchyPotential(double eps, double sigl, double rcut,
                    double siga1, double siga2, double sigt1, double sigt2);
            eneT calc_energy(distT rdist, vecT patch_norm1, vecT patch_norm2,
                    vecT patch_orient1, vecT patch_orient2);
            bool particles_interacting(distT rdist);

        private:
            PatchyPotential m_patchy; // Unoriented-patchy potential
            double m_sigt1; // Orientation width 1
            double m_sigt2; // Orientation width 2
    };
}

#endif // POTENTIAL_H
