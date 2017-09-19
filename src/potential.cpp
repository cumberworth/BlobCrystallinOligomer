// potential.cpp

#include <cmath>
#include <iostream>

#include "BlobCrystallinOligomer/potential.h"

namespace potential {

    using shared_types::inf;
    using std::acos;
    using std::cout;
    using std::exp;
    using std::pow;

    eneT gaussian(double theta, double sig) {

        // Should do a fast version that takes presquared values
        return exp(-pow(theta, 2)/(2*pow(sig, 2)));
    }

    PairPotential::PairPotential(double rcut): m_rcut {rcut} {}

    bool PairPotential::particles_interacting(distT rdist) {
        bool interacting {false};
        if (rdist < m_rcut) {
            interacting = true;
        }

        return interacting;
    }

    HardSpherePotential::HardSpherePotential(double sigh):
            PairPotential {sigh}, m_sigh {sigh} {
    }

    eneT HardSpherePotential::calc_energy(distT rdist, vecT&, Orientation&, Orientation&) {
        eneT ene {0};
        if (rdist < m_sigh) {
            ene = inf;
        }

        return ene;
    }

    ShiftedLJPotential::ShiftedLJPotential(double eps, double sigl,
            double rcut):
            PairPotential {rcut}, m_eps {eps}, m_four_eps {4*eps},
            m_sigl {sigl}, m_rcut {rcut} {

        distT sig_r_ratio {m_sigl/rcut};
        m_shift = m_four_eps*(pow(sig_r_ratio, 12) - pow(sig_r_ratio, 6));
    }

    eneT ShiftedLJPotential::calc_energy(distT rdist, vecT&, Orientation&,
            Orientation&) {

        eneT ene;
        if (rdist >= m_rcut) {
            ene = 0;
            return ene;
        }
        distT sig_r_ratio {m_sigl/rdist};
        ene = m_four_eps*(pow(sig_r_ratio, 12) - pow(sig_r_ratio, 6)) - m_shift;

        return ene;
    }

    PatchyPotential::PatchyPotential(double eps, double sigl, double rcut, 
            double siga1, double siga2):
            PairPotential {rcut}, m_lj {eps, sigl, rcut},
            m_sigl {sigl},
            m_siga1 {siga1},
            m_siga2 {siga2} {
    }

    eneT PatchyPotential::calc_energy(distT rdist, vecT& p_diff, Orientation& ore1,
            Orientation& ore2) {

        eneT ene {m_lj.calc_energy(rdist, p_diff, ore1, ore2)};
        if (rdist < m_sigl or ene == 0) {
            return ene;
        }
        vecT p_diff_unit_ij {p_diff / rdist};
        vecT p_diff_unit_ji {-p_diff / rdist};
        distT dot1 {p_diff_unit_ij.dot(ore1.patch_norm)};
        distT dot2 {p_diff_unit_ji.dot(ore2.patch_norm)};
        // TODO fix this shit
        if (dot1 > 1) {
            dot1 = 1;
        }
        if (dot1 < -1) {
            dot1 = -1;
        }
        if (dot2 > 1) {
            dot2 = 1;
        }
        if (dot2 < -1) {
            dot2 = -1;
        }
        distT theta1 {acos(dot1)};
        distT theta2 {acos(dot2)};
        ene *= gaussian(theta1, m_siga1);
        ene *= gaussian(theta2, m_siga2);

        return ene;
    }

    OrientedPatchyPotential::OrientedPatchyPotential(double eps, double sigl,
            double rcut, double siga1, double siga2, double sigt):
            PairPotential {rcut}, m_patchy {eps, sigl, rcut, siga1, siga2},
            m_sigl {sigl},
            m_sigt {sigt} {
    }

    eneT OrientedPatchyPotential::calc_energy(distT rdist, vecT& p_diff,
            Orientation& ore1, Orientation& ore2) {

        eneT ene {m_patchy.calc_energy(rdist, p_diff, ore1, ore2)};
        if (rdist < m_sigl or ene == 0) {
            return ene;
        }
        vecT p_diff_unit_ij {p_diff / p_diff.norm()};
        vecT p_diff_unit_ji {-p_diff_unit_ij};
        vecT proj1 {ore1.patch_orient.dot(p_diff_unit_ij)*p_diff_unit_ij};
        vecT proj2 {ore2.patch_orient.dot(p_diff_unit_ji)*p_diff_unit_ji};
        vecT rej1 {ore1.patch_orient - proj1};
        vecT rej2 {ore2.patch_orient - proj2};
        distT proj_dot {rej1.dot(rej2)};
        // TODO fix this shit
        distT rat {proj_dot / (rej1.norm()*rej2.norm())};
        if (rat > 1) {
            rat = 1;
        }
        if (rat < -1) {
            rat = -1;
        }
        distT theta {acos(rat)};
        ene *= gaussian(theta, m_sigt);

        return ene;
    }
}
