// potential.cpp

#include <cmath>

#include "BlobCrystallinOligomer/potential.h"
#include "BlobCrystallinOligomer/shared_types.h"

namespace potential {

    using shared_types::eneT;
    using shared_types::distT;
    using shared_types::inf;
    using shared_types::vecT;
    using std::acos;
    using std::exp;
    using std::pow;

    eneT gaussian(double theta, double sig) {
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
        distT sig_r_ratio {m_sigl/rdist};
        eneT ene;
        ene = m_four_eps*(pow(sig_r_ratio, 12) - pow(sig_r_ratio, 6)) - m_shift;

        return ene;
    }

    PatchyPotential::PatchyPotential(double eps, double sigl, double rcut, 
            double siga1, double siga2):
            PairPotential {rcut}, m_lj {eps, sigl, rcut}, m_siga1 {siga1},
            m_siga2 {siga2} {
    }

    eneT PatchyPotential::calc_energy(distT rdist, vecT& p_diff, Orientation& ore1,
            Orientation& ore2) {

        eneT ene {0};
        distT dot1 {p_diff.dot(ore1.patch_norm)};
        distT dot2 {p_diff.dot(ore2.patch_norm)};
        distT theta1 {acos(dot1 / (p_diff.norm()*ore1.patch_norm.norm()))};
        distT theta2 {acos(dot2 / (p_diff.norm()*ore2.patch_norm.norm()))};
        ene += m_lj.calc_energy(rdist, p_diff, ore1, ore2);
        ene *= gaussian(theta1, m_siga1);
        ene *= gaussian(theta2, m_siga2);

        return ene;
    }

    OrientedPatchyPotential::OrientedPatchyPotential(double eps, double sigl,
            double rcut, double siga1, double siga2, double sigt):
            PairPotential {rcut}, m_patchy {eps, sigl, rcut, siga1, siga2},
            m_sigt {sigt} {
    }

    eneT OrientedPatchyPotential::calc_energy(distT rdist, vecT& p_diff,
            Orientation& ore1, Orientation& ore2) {

        eneT ene {0};
        vecT unit_p_diff {p_diff / p_diff.norm()};
        vecT proj1 {ore1.patch_orient.dot(unit_p_diff)*unit_p_diff};
        vecT proj2 {ore2.patch_orient.dot(unit_p_diff)*unit_p_diff};
        vecT rej1 {ore1.patch_orient - proj1};
        vecT rej2 {ore2.patch_orient - proj2};
        distT proj_dot {rej1.dot(rej2)};
        distT theta {acos(proj_dot / (rej1.norm()*rej2.norm()))};
        ene += m_patchy.calc_energy(rdist, p_diff, ore1, ore2);
        ene *= gaussian(m_sigt, theta);

        return ene;
    }
}
