// movetype.cpp

#include <algorithm>
#include <cmath>
#include <vector>

#include "BlobCrystallinOligomer/monomer.h"
#include "BlobCrystallinOligomer/movetype.h"
#include "BlobCrystallinOligomer/shared_types.h"

namespace movetype {

    using config::monomerArrayT;
    using monomer::Monomer;
    using shared_types::eneT;
    using shared_types::CoorSet;
    using std::exp;
    using std::fmin;
    using std::fmax;
    using std::vector;
    using std::find;

    MCMovetype::MCMovetype(Config& conf, Energy& ene, RandomGens& random_num,
            eneT beta):
            m_config {conf}, m_energy {ene}, m_random_num {random_num},
            m_beta {beta} {
    }

    bool VMMCMovetype::move() {
        Monomer& monomer_seed {m_config.get_random_monomer()};
        generate_movemap();
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
            double prelink_for_p {calc_prelink_prob(ene_2 - ene_1)};
            bool prelink_accepted {accept_prelink(prelink_for_p)};
            if (not prelink_accepted) {
                continue;
            }
            eneT ene_3 {m_energy.calc_monomer_pair_energy(monomer1,
                    CoorSet::trial, monomer2, CoorSet::current)};
            double prelink_rev_p {calc_prelink_prob(ene_3 - ene_1)};
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

    double VMMCMovetype::calc_prelink_prob(eneT delta_e) {
        return fmax(0, 1 - exp(-m_beta * delta_e));
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

    void TranslationVMMCMovetype::generate_movemap() {
    }

    void TranslationVMMCMovetype::apply_movemap(Monomer& monomer) {
    }

    void RotationVMMCMovetype::generate_movemap() {
    }

    void RotationVMMCMovetype::apply_movemap(Monomer& monomer) {
    }
}
