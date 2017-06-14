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
    eneT guassian(double theta, double sig);

    class PairPotential {
        public:
            PairPotential(double rcut);
            virtual eneT calc_energy(distT rdist, vecT p_diff, Orientation ore1,
                    Orientation ore2) = 0;
            bool particles_interacting(distT rdist);

        private:
            double m_rcut;
    };

    class HardSpherePotential: public PairPotential {
        public:
            HardSpherePotential(double sigh);
            eneT calc_energy(distT rdist, vecT, Orientation, Orientation);
        private:
            double m_sigh; // Sphere radius
    };

    class ShiftedLJPotential: public PairPotential {
        public:
            ShiftedLJPotential(double eps, double sigl, double rcut);
            eneT calc_energy(distT rdist, vecT, Orientation ore1, Orientation ore2);

        private:
            double m_eps; // Well depth
            double m_sigl; // Zero point
            double m_rcut; // Cutoff
            double m_shift; // Shift
    };

    class PatchyPotential: public PairPotential {
        public:
            PatchyPotential(double eps, double sigl, double rcut, double siga1,
                    double siga2);
            eneT calc_energy(distT rdist, vecT p_diff, Orientation ore1,
                    Orientation ore2);

        private:
            ShiftedLJPotential m_lj; // Radial component
            double m_siga1; // Patch width of particle 1
            double m_siga2; // Patch width of particle 2
    };

    class OrientedPatchyPotential: public PairPotential {
        public:
            OrientedPatchyPotential(double eps, double sigl, double rcut,
                    double siga1, double siga2, double sigt);
            eneT calc_energy(distT rdist, vecT p_diff, Orientation ore1,
                    Orientation ore2);

        private:
            PatchyPotential m_patchy; // Unoriented-patchy potential
            double m_sigt; // Orientation width
    };
}

#endif // POTENTIAL_H
