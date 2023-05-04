// movetype.cpp

#include <algorithm>
#include <cmath>
#include <vector>

#include "Eigen/Geometry"

#include "BlobCrystallinOligomer/movetype.h"

namespace movetype {

using Eigen::AngleAxis;
using shared_types::CoorSet;
using shared_types::inf;
using std::cout;
using std::exp;
using std::find;
using std::fmax;
using std::fmin;

vecT random_unit_vector(RandomGens& random_num) {
    /*  Taken from Daan's book, which is taken from Allen and Tildesley */
    distT ransq {2};
    distT ran1;
    distT ran2;
    while (ransq >= 1) {
        ran1 = 1 - 2 * random_num.uniform_real();
        ran2 = 1 - 2 * random_num.uniform_real();
        ransq = ran1 * ran1 + ran2 * ran2;
    }
    distT ranh {2 * sqrt(1 - ransq)};
    distT x {ran1 * ranh};
    distT y {ran2 * ranh};
    distT z {1 - 2 * ransq};

    return {x, y, z};
}

distT random_displacement(distT max_disp, RandomGens& random_num) {
    // isn't this only giving between -max_disp/2 and max_disp/2?
    return max_disp * (random_num.uniform_real() - 0.5);
}

Movemap::Movemap(RandomGens& random_num): m_random_num {random_num} {}

TranslationMovemap::TranslationMovemap(
        distT max_disp_tc,
        RandomGens& random_num):
        Movemap {random_num}, m_max_disp_tc {max_disp_tc} {}

void TranslationMovemap::generate_movemap(Monomer&) {
    for (size_t i {0}; i != 3; i++) {
        m_disp_v[i] = random_displacement(m_max_disp_tc, m_random_num);
    }
}

void TranslationMovemap::apply_movemap(Monomer& monomer) {
    monomer.translate(m_disp_v);
}

RotationMovemap::RotationMovemap(
        distT max_disp_rc,
        distT max_disp_a,
        RandomGens& random_num):
        Movemap {random_num},
        m_max_disp_rc {max_disp_rc},
        m_max_disp_a {max_disp_a} {}

void RotationMovemap::generate_movemap(Monomer& monomer) {
    vecT rand_v {random_unit_vector(m_random_num)};
    distT scalar {random_displacement(m_max_disp_rc, m_random_num)};
    m_rot_c = monomer.get_center(CoorSet::current) + scalar * rand_v;
    vecT axis {random_unit_vector(m_random_num)};
    distT theta {random_displacement(m_max_disp_a, m_random_num)};
    AngleAxis<distT> angle_axis {theta, axis};
    m_rot_mat = angle_axis.toRotationMatrix();
}

void RotationMovemap::apply_movemap(Monomer& monomer) {
    monomer.rotate(m_rot_c, m_rot_mat);
}

NTDFlipMovemap::NTDFlipMovemap(Config& config, RandomGens& random_num):
        Movemap {random_num}, m_config {config} {

    m_imat << 1, 0, 0, 0, 1, 0, 0, 0, 1;
}

void NTDFlipMovemap::generate_movemap(Monomer& monomer) {
    auto particles = monomer.get_particles();
    vecT plane_normal;
    auto r = m_random_num.uniform_real();
    if (r < 0.25) {
        auto p = particles[0].get();
        plane_normal = p.get_ore(CoorSet::current).patch_orient;
        m_point_in_plane = p.get_pos(CoorSet::current);
    }
    else if (r < 0.5) {
        auto p = particles[2].get();
        plane_normal = p.get_ore(CoorSet::current).patch_norm;
        m_point_in_plane = p.get_pos(CoorSet::current);
    }
    else if (r < 0.75) {
        auto p1 = particles[0].get();
        auto p2 = particles[1].get();
        auto axis = m_config.calc_interparticle_vector(
                p1, CoorSet::current, p2, CoorSet::current);
        axis.normalize();
        plane_normal = p1.get_ore(CoorSet::current).patch_orient;
        AngleAxis<distT> angle_axis {M_PI / 2, axis};
        plane_normal = angle_axis.toRotationMatrix() * plane_normal;
        m_point_in_plane = p1.get_pos(CoorSet::current);
    }
    else {
        auto p1 = particles[2].get();
        auto p2 = particles[3].get();
        auto axis = m_config.calc_interparticle_vector(
                p1, CoorSet::current, p2, CoorSet::current);
        axis.normalize();
        plane_normal = p1.get_ore(CoorSet::current).patch_norm;
        AngleAxis<distT> angle_axis {M_PI / 2, axis};
        plane_normal = angle_axis.toRotationMatrix() * plane_normal;
        m_point_in_plane = p1.get_pos(CoorSet::current);
    }
    plane_normal.normalize();
    m_ref_mat = m_imat - 2 * plane_normal * plane_normal.transpose();
}

void NTDFlipMovemap::apply_movemap(Monomer& monomer) {
    monomer.rotate(m_point_in_plane, m_ref_mat);
    monomer.flip_conformation();
}

MCMovetype::MCMovetype(
        Config& conf,
        Energy& ene,
        RandomGens& random_num,
        InputParams params,
        string label):
        m_config {conf},
        m_energy {ene},
        m_random_num {random_num},
        m_beta {1 / params.m_temp},
        m_label {label} {}

string MCMovetype::get_label() { return m_label; }

MetMCMovetype::MetMCMovetype(
        Config& conf,
        Energy& ene,
        RandomGens& random_num,
        InputParams params,
        string label,
        string movemap_type):
        MCMovetype {conf, ene, random_num, params, label} {

    if (movemap_type == "translation") {
        m_movemap = std::make_unique<TranslationMovemap>(
                params.m_max_disp_tc, random_num);
    }
    else if (movemap_type == "rotation") {
        m_movemap = std::make_unique<RotationMovemap>(
                params.m_max_disp_rc, params.m_max_disp_a, random_num);
    }
    else if (movemap_type == "ntdflip") {
        m_movemap = std::make_unique<NTDFlipMovemap>(conf, random_num);
    }
}

bool MetMCMovetype::move() {
    Monomer& m {m_config.get_random_monomer()};
    m_movemap->generate_movemap(m);
    m_movemap->apply_movemap(m);
    eneT de {m_energy.calc_monomer_diff(m)};
    bool accepted {accept_move(de)};
    if (accepted) {
        m.trial_to_current();
    }
    else {
        m.current_to_trial();
    }

    return accepted;
}

bool MetMCMovetype::accept_move(eneT de) {
    bool accept;
    if (de == inf) {
        accept = false;
    }
    else {
        eneT paccept {fmin(1, exp(-m_beta * de))};
        if (paccept == 1) {
            accept = true;
        }
        else {
            if (paccept > m_random_num.uniform_real()) {
                accept = true;
            }
            else {
                accept = false;
            }
        }
    }

    return accept;
}

VMMCMovetype::VMMCMovetype(
        Config& conf,
        Energy& ene,
        RandomGens& random_num,
        InputParams params,
        string label,
        string movemap_type):
        MCMovetype {conf, ene, random_num, params, label} {

    if (movemap_type == "translation") {
        m_movemap = std::make_unique<TranslationMovemap>(
                params.m_max_disp_tc, random_num);
    }
    else if (movemap_type == "rotation") {
        m_movemap = std::make_unique<RotationMovemap>(
                params.m_max_disp_rc, params.m_max_disp_a, random_num);
    }
}

bool VMMCMovetype::move() {
    Monomer& monomer_seed {m_config.get_random_monomer()};
    m_cluster.emplace_back(monomer_seed);
    m_interacting_mis.insert(monomer_seed.get_index());
    m_movemap->generate_movemap(monomer_seed);
    m_movemap->apply_movemap(monomer_seed);
    add_interacting_pairs(monomer_seed);
    while (m_pair_mis.size() != 0) {
        pair<int, int> cur_pair {pop_random_pair()};
        Monomer& monomer1 {m_config.get_monomer(cur_pair.first)};
        Monomer& monomer2 {m_config.get_monomer(cur_pair.second)};
        bool mono_in_cluster {false};
        for (auto mono_cluster: m_cluster) {
            if (monomer2.get_index() == mono_cluster.get().get_index()) {
                mono_in_cluster = true;
                break;
            }
        }
        if (mono_in_cluster) {
            continue;
        }
        eneT ene_1 {m_energy.calc_monomer_pair_energy(
                monomer1, CoorSet::current, monomer2, CoorSet::current)};
        eneT ene_2 {m_energy.calc_monomer_pair_energy(
                monomer1, CoorSet::trial, monomer2, CoorSet::current)};
        double prelink_for_p {calc_prelink_prob(ene_1, ene_2)};
        bool prelink_accepted {accept_prelink(prelink_for_p)};
        if (not prelink_accepted) {
            continue;
        }
        eneT ene_3 {m_energy.calc_monomer_pair_energy(
                monomer1, CoorSet::current, monomer2, CoorSet::trial)};
        double prelink_rev_p {calc_prelink_prob(ene_1, ene_3)};
        bool link_accepted {accept_link(prelink_for_p, prelink_rev_p)};
        if (not link_accepted) {
            m_frustrated_links++;
            m_frustrated_mis.push_back(monomer2.get_index());
            continue;
        }
        auto frustrated_mi {
                find(m_frustrated_mis.begin(),
                     m_frustrated_mis.end(),
                     monomer2.get_index())};
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
    for (auto i: m_interacting_mis) {
        Monomer& mono {m_config.get_monomer(i)};
        mono.current_to_trial();
    }
    reset_internal();

    return accepted;
}

void VMMCMovetype::add_interacting_pairs(Monomer& monomer1) {

    // Get all monomers that are interacting before and after movemap
    monomerArrayT monomers {
            m_energy.get_interacting_monomers(monomer1, CoorSet::current)};
    monomerArrayT trial_monomers {
            m_energy.get_interacting_monomers(monomer1, CoorSet::trial)};
    monomers.insert(
            monomers.end(), trial_monomers.begin(), trial_monomers.end());

    for (Monomer& mono_any: monomers) {

        // The second monomer should not already be in the cluster
        int mono_i2 {mono_any.get_index()};
        bool mono_in_cluster {false};
        for (auto mono_cluster: m_cluster) {
            if (mono_i2 == mono_cluster.get().get_index()) {
                mono_in_cluster = true;
                break;
            }
        }
        if (mono_in_cluster) {
            continue;
        }

        // The pair should not have been proposed before
        int mono_i1 {monomer1.get_index()};
        pair<int, int> pair_mis {mono_i1, mono_i2};
        auto proposed {not m_proposed_pairs.insert(pair_mis).second};
        if (proposed) {
            continue;
        }

        // Only apply the movemap once TODO only apply after prelink accepted
        bool movemap_applied {not m_interacting_mis.insert(mono_i2).second};
        if (not movemap_applied) {
            m_movemap->apply_movemap(mono_any);
        }
        m_pair_mis.insert(pair_mis);
    }
}

pair<int, int> VMMCMovetype::pop_random_pair() {
    int pair_i {m_random_num.uniform_int(0, m_pair_mis.size() - 1)};
    set<pair<int, int>>::const_iterator it {m_pair_mis.begin()};
    std::advance(it, pair_i);
    pair<int, int> sel_pair {*it};
    m_pair_mis.erase(it);

    return sel_pair;
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
    double p_accept {fmin(1, prelink_rev_p / prelink_for_p)};
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
    m_proposed_pairs.clear();
    m_interacting_mis.clear();
    m_pair_mis.clear();
}
} // namespace movetype
