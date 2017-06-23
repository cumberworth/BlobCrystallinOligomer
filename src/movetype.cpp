// movetype.cpp

#include <algorithm>
#include <cmath>
#include <vector>

#include "Eigen/Geometry"

#include "BlobCrystallinOligomer/config.h"
#include "BlobCrystallinOligomer/monomer.h"
#include "BlobCrystallinOligomer/movetype.h"
#include "BlobCrystallinOligomer/param.h"
#include "BlobCrystallinOligomer/shared_types.h"

namespace movetype {

    using config::monomerArrayT;
    using Eigen::AngleAxis;
    using monomer::Monomer;
    using param::InputParams;
    using shared_types::distT;
    using shared_types::eneT;
    using shared_types::inf;
    using shared_types::vecT;
    using shared_types::CoorSet;
    using std::exp;
    using std::fmin;
    using std::fmax;
    using std::vector;
    using std::find;

    vecT random_unit_vector(RandomGens& random_num) {
        /*  Taken from Daan's book, which is taken from Allen and Tildesley */
        distT ransq {2};
        distT ran1;
        distT ran2;
        while (ransq >= 1) {
            ran1 = 1 - 2*random_num.uniform_real();
            ran2 = 1 - 2*random_num.uniform_real();
            ransq = ran1*ran1 + ran2*ran2;
        }
        distT ranh {2*sqrt(1 - ransq)};
        distT x {ran1*ranh};
        distT y {ran2*ranh};
        distT z {1 - 2*ransq};

        return {x, y, z};
    }

	distT random_displacement(distT max_disp, RandomGens& random_num) {
        return max_disp * (random_num.uniform_real() - 0.5);
    }

    MCMovetype::MCMovetype(Config& conf, Energy& ene, RandomGens& random_num,
            InputParams params):
            m_config {conf}, m_energy {ene}, m_random_num {random_num},
            m_beta {1/params.m_temp} {
    }

    bool VMMCMovetype::move() {
        Monomer& monomer_seed {m_config.get_random_monomer()};
        m_cluster.emplace_back(monomer_seed);
        generate_movemap(monomer_seed);
        apply_movemap(monomer_seed);
        add_interacting_pairs(monomer_seed);
        while (m_pair_mis.size() != 0) {
            pair<int, int> cur_pair {select_random_pair()};
            Monomer& monomer1 {m_config.get_monomer(cur_pair.first)};
            Monomer& monomer2 {m_config.get_monomer(cur_pair.second)};
            eneT ene_1 {m_energy.calc_monomer_pair_energy(monomer1,
                    CoorSet::current, monomer2, CoorSet::current)};
            eneT ene_2 {m_energy.calc_monomer_pair_energy(monomer1,
                    CoorSet::trial, monomer2, CoorSet::current)};
            double prelink_for_p {calc_prelink_prob(ene_1, ene_2)};
            bool prelink_accepted {accept_prelink(prelink_for_p)};
            if (not prelink_accepted) {
                continue;
            }
            eneT ene_3 {m_energy.calc_monomer_pair_energy(monomer1,
                    CoorSet::trial, monomer2, CoorSet::current)};
            double prelink_rev_p {calc_prelink_prob(ene_1, ene_3)};
            bool link_accepted {accept_link(prelink_for_p, prelink_rev_p)};
            if (not link_accepted) {
                m_frustrated_links++;
                m_frustrated_mis.push_back(monomer2.get_index());
                continue;
            }
            auto frustrated_mi {find(m_frustrated_mis.begin(),
                    m_frustrated_mis.end(), monomer2.get_index())};
            if (frustrated_mi != m_frustrated_mis.end()) {
                m_frustrated_links--;
                m_frustrated_mis.erase(frustrated_mi);
            }
            m_cluster.emplace_back(monomer2);
            add_interacting_pairs(monomer2);
        }
        bool accepted {accept_move()};
        if (accepted) {
            for (Monomer& mono: m_cluster) {
                mono.trial_to_current();
            }
        }
        reset_internal();

        return accepted;
    }

    void VMMCMovetype::add_interacting_pairs(Monomer& monomer1) {
        monomerArrayT monomers {m_energy.get_interacting_monomers(monomer1,
                CoorSet::current)};
        monomerArrayT trial_monomers {m_energy.get_interacting_monomers(monomer1,
                CoorSet::trial)};
        monomers.insert(monomers.end(), trial_monomers.begin(),
                trial_monomers.end());
        set<int> unique_mis {};
        for (Monomer& mono: monomers) {
            int mono_i2 {mono.get_index()};
            auto insert_result {unique_mis.insert(mono_i2)};
            if (insert_result.second) {
                insert_result = m_interacting_mis.insert(mono_i2);
                if (not insert_result.second) {
                    apply_movemap(mono);
                }
            pair<int, int> pair_mis;
            int mono_i1 {monomer1.get_index()};
            if (mono_i2 > mono_i1) {
                pair_mis = {mono_i1, mono_i2};
            }
            else {
                pair_mis = {mono_i2, mono_i1};
            }
            m_pair_mis.insert(pair_mis);
            }
        }
    }

    pair<int, int> VMMCMovetype::select_random_pair() {
        int pair_i {m_random_num.uniform_int(0, m_pair_mis.size())};
        set<pair<int, int>>::const_iterator it {m_pair_mis.begin()};
        std::advance(it, pair_i);

        return *it;
    }

    double VMMCMovetype::calc_prelink_prob(eneT ene1, eneT ene2) {
        // ene1 should never be infinite
        if (ene2 == inf) {
            return 1;
        }
        else {
            return fmax(0, 1 - exp(-m_beta * (ene2 - ene1)));
        }
    }

    bool VMMCMovetype::accept_prelink(double prelink_p) {
        bool accept;
        if (prelink_p == 0) {
            accept = false;
        }
        else {
            if (prelink_p > m_random_num.uniform_real()) {
                accept = true;
            }
            else {
                accept = false;
            }
        }

        return accept;
    }

    bool VMMCMovetype::accept_link(double prelink_for_p, double prelink_rev_p) {
        double p_accept {fmin(1, prelink_rev_p/prelink_for_p)};
        bool accept;
        if (p_accept == 1) {
            accept = true;
        }
        else {
            if (p_accept > m_random_num.uniform_real()) {
                accept = true;
            }
            else {
                accept = false;
            }
        }

        return accept;
    }

    bool VMMCMovetype::accept_move() {
        if (m_frustrated_links != 0) {
            return false;
        }
        else {
            return true;
        }
    }

    void VMMCMovetype::reset_internal() {
        m_cluster.clear();
        m_frustrated_links = 0;
        m_frustrated_mis.clear();
        m_interacting_mis.clear();
        m_pair_mis.clear();
    }

    TranslationVMMCMovetype::TranslationVMMCMovetype(Config& conf, Energy& ene,
            RandomGens& random_num, InputParams params):
            VMMCMovetype {conf, ene, random_num, params},
                     m_max_disp_tc {params.m_max_disp_tc} {
    }

    void TranslationVMMCMovetype::generate_movemap(Monomer&) {
        for (size_t i {0}; i != 3; i++) {
            m_disp_v[i] = random_displacement(m_max_disp_tc, m_random_num);
        }
    }

    void TranslationVMMCMovetype::apply_movemap(Monomer& monomer) {
        monomer.translate(m_disp_v);
    }

    RotationVMMCMovetype::RotationVMMCMovetype(Config& conf, Energy& ene, RandomGens& random_num,
            InputParams params):
            VMMCMovetype {conf, ene, random_num, params},
                     m_max_disp_rc {params.m_max_disp_rc},
                     m_max_disp_a {params.m_max_disp_a} {
    }

    void RotationVMMCMovetype::generate_movemap(Monomer& seed_monomer) {
        vecT rand_v {random_unit_vector(m_random_num)};
        distT scalar {random_displacement(m_max_disp_rc, m_random_num)};
        m_rot_c = seed_monomer.get_center() + scalar*rand_v;
        vecT axis {random_unit_vector(m_random_num)};
        distT theta {random_displacement(m_max_disp_a, m_random_num)};
        AngleAxis<distT> angle_axis {theta, axis};
        m_rot_mat = angle_axis.toRotationMatrix();
    }

    void RotationVMMCMovetype::apply_movemap(Monomer& monomer) {
        monomer.rotate(m_rot_c, m_rot_mat);
    }
}
