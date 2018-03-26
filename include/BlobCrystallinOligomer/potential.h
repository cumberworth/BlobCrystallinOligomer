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

    /** Return value of Gaussian function */
    eneT guassian(distT theta, distT sig);

    /** Return dihedral angle */
    distT dihedral(vecT ore1, vecT ore2, vecT p_diff);

    /** Interface and shared implementation to pair potentials */
    class PairPotential {
        public:
            PairPotential(distT rcut);

            /** Calculate pair potential */
            virtual eneT calc_energy(
                    distT rdist,
                    vecT& p_diff,
                    Orientation& ore1,
                    Orientation& ore2) = 0;

            /** Check if pair potential is non-zero */
            bool particles_interacting(distT rdist);

        private:
            distT m_rcut;
    };

    /** No interaction */
    class ZeroPotential:
            public PairPotential {

        public:
            ZeroPotential();
            eneT calc_energy(distT, vecT&, Orientation&, Orientation&);
    };

    /** Hard sphere potential */
    class HardSpherePotential:
            public PairPotential {

        public:
            HardSpherePotential(distT sigh);
            eneT calc_energy(distT rdist, vecT&, Orientation&, Orientation&);

        private:
            distT m_sigh; // Sphere radius
    };

    /** Square well potential */
    class SquareWellPotential:
            public PairPotential {

        public:
            SquareWellPotential(eneT eps, distT rcut);
            eneT calc_energy(distT rdist, vecT&, Orientation&, Orientation&);

        private:
            eneT m_eps;
            distT m_rcut;
    };

    /** Harmonic well potential */
    class HarmonicWellPotential:
            public PairPotential {

        public:
            HarmonicWellPotential(eneT eps, distT rcut);
            eneT calc_energy(
                    distT rdist,
                    vecT& p_diff,
                    Orientation& ore1,
                    Orientation& ore2);

        private:
            eneT m_eps;
            distT m_rcut;
            distT m_a; // steepness parameter
    };

    /** Harmonic well potential */
    class AngularHarmonicWellPotential:
            public PairPotential {

        public:
            AngularHarmonicWellPotential(eneT eps, distT rcut, distT siga);
            eneT calc_energy(distT rdist, vecT&, Orientation&, Orientation&);

        private:
            HarmonicWellPotential m_hwell;
            eneT m_eps;
            distT m_rcut;
            distT m_siga;
    };

    /** Shifted Lennard Jone potential */
    class ShiftedLJPotential:
            public PairPotential {

        public:
            ShiftedLJPotential(eneT eps, distT sigl, distT rcut);
            eneT calc_energy(
                    distT rdist,
                    vecT&,
                    Orientation& ore1,
                    Orientation& ore2);

        private:
            eneT m_eps; // Well depth
            eneT m_four_eps; // Well depth premultiplied by 4
            distT m_sigl; // Zero point
            distT m_rcut; // Cutoff
            eneT m_shift; // Shift
    };

    /** Directional patch potential
      *
      * Taken from [doye paper]
      */
    class PatchyPotential: public PairPotential {
        public:
            PatchyPotential(eneT eps, distT sigl, distT rcut, distT siga1,
                    distT siga2);
            eneT calc_energy(distT rdist, vecT& p_diff, Orientation& ore1,
                    Orientation& ore2);

        private:
            ShiftedLJPotential m_lj; // Radial component
            distT m_sigl;
            distT m_siga1; // Patch width of particle 1
            distT m_siga2; // Patch width of particle 2
    };

    /** Oriented patch potential
      *
      * Taken from [doye paper]
      */
    class OrientedPatchyPotential: public PairPotential {
        public:
            OrientedPatchyPotential(eneT eps, distT sigl, distT rcut,
                    distT siga1, distT siga2, distT sigt);
            eneT calc_energy(distT rdist, vecT& p_diff, Orientation& ore1,
                    Orientation& ore2);

        private:
            PatchyPotential m_patchy; // Unoriented-patchy potential
            distT m_sigl;
            distT m_sigt; // Orientation width
    };

    /** Oriented patch potential
      *
      * One more vector added to oriented patchy
      */
    class DoubleOrientedPatchyPotential: public PairPotential {
        public:
            DoubleOrientedPatchyPotential(eneT eps, distT sigl, distT rcut,
                    distT siga1, distT siga2, distT sigt);
            eneT calc_energy(distT rdist, vecT& p_diff, Orientation& ore1,
                    Orientation& ore2);

        private:
            PatchyPotential m_patchy; // Unoriented-patchy potential
            distT m_sigl;
            distT m_sigt; // Orientation width
    };
}

#endif // POTENTIAL_H
