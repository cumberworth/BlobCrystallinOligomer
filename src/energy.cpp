// energy.cpp

#include <memory>
#include <vector>

#include "BlobCrystallinOligomer/energy.h"

namespace energy {

    using ifile::InputEnergyFile;
    using monomer::particleArrayT;
    using potential::ZeroPotential;
    using potential::HardSpherePotential;
    using potential::SquareWellPotential;
    using potential::HarmonicWellPotential;
    using potential::AngularHarmonicWellPotential;
    using potential::ShiftedLJPotential;
    using potential::PatchyPotential;
    using potential::OrientedPatchyPotential;
    using potential::DoubleOrientedPatchyPotential;
    using shared_types::inf;
    using shared_types::InputError;
    using shared_types::vecT;
    using std::cout;

    Energy::Energy(Config& conf, InputParams& params):
            m_config {conf},
            m_max_cutoff {params.m_max_cutoff} {

        InputEnergyFile energy_file {params.m_energy_filename};
        vector<PotentialData> potentials {energy_file.get_potentials()};
        vector<InteractionData> same_conformers_interactions {
                energy_file.get_same_conformers_interactions()};
        vector<InteractionData> different_conformers_interactions {
                energy_file.get_different_conformers_interactions()};
        create_potentials(potentials, same_conformers_interactions,
                different_conformers_interactions);
        eneT total_ene {calc_total_energy()};
        if (total_ene == inf or total_ene != total_ene) {
            cout << "Bad starting configuration\n";
            throw InputError {};
        }
    }

    eneT Energy::calc_total_energy() {
        monomerArrayT monomers {m_config.get_monomers()};
        eneT total_ene {0};
        for (size_t i {0}; i != monomers.size() - 1; i++) {
            Monomer& monomer1 {monomers[i].get()};
            for (size_t j {i + 1}; j != monomers.size(); j++) {
                Monomer& monomer2 {monomers[j].get()};
                total_ene += calc_monomer_pair_energy(monomer1, CoorSet::current,
                        monomer2, CoorSet::current);
            }
        }

        return total_ene;
    }

    eneT Energy::calc_monomer_pair_energy(Monomer& monomer1, CoorSet coorset1,
            Monomer& monomer2, CoorSet coorset2) {

        eneT pair_ene {0};
        particleArrayT particles1 {monomer1.get_particles()};
        particleArrayT particles2 {monomer2.get_particles()};
        for (Particle& p1: particles1) {
            for (Particle& p2: particles2) {
                eneT part_ene {calc_particle_pair_energy(p1,
                        monomer1.get_conformer(coorset1), coorset1, p2,
                        monomer2.get_conformer(coorset2), coorset2)};
                if (part_ene == inf) {
                    return inf;
                }
                pair_ene += part_ene;
            }
        }

        return pair_ene;
    }

    bool Energy::monomers_interacting(Monomer& monomer1, CoorSet coorset1,
            Monomer& monomer2, CoorSet coorset2) {

        bool m_interacting {false};
        if (not monomers_in_range(monomer1, coorset1, monomer2, coorset2)) {
            return m_interacting;
        }
        particleArrayT particles1 {monomer1.get_particles()};
        particleArrayT particles2 {monomer2.get_particles()};
        for (Particle& p1: particles1) {
            for (Particle& p2: particles2) {
                bool p_interacting {particles_interacting(p1,
                        monomer1.get_conformer(coorset1), coorset1, p2,
                        monomer2.get_conformer(coorset2), coorset2)};
                if (p_interacting) {
                    m_interacting = true;
                    return m_interacting;
                }
            }
        }

        return m_interacting;
    }

    bool Energy::monomers_in_range(Monomer& monomer1, CoorSet coorset1,
            Monomer& monomer2, CoorSet coorset2) {

        bool in_range {false};
        distT monomer1_r {monomer1.get_radius()};
        distT monomer2_r {monomer2.get_radius()};
        distT max_interaction_d {monomer1_r + monomer2_r + m_max_cutoff};
        distT d {m_config.calc_dist(monomer1, coorset1, monomer2, coorset2)};
        if (d <= max_interaction_d) {
            in_range = true;
        }

        return in_range;
    }

    monomerArrayT Energy::get_interacting_monomers(Monomer& monomer1,
            CoorSet coorset1) {

        monomerArrayT monomers {m_config.get_monomers()};
        monomerArrayT interacting_monomers {};
        CoorSet coorset2 {CoorSet::current};
        for (size_t i {0}; i != monomers.size(); i++) {
            Monomer& monomer2 {monomers[i].get()};
            if (monomer1.get_index() == monomer2.get_index()) {
                continue;
            }
            if (monomers_interacting(monomer1, coorset1, monomer2, coorset2)) {
                interacting_monomers.push_back(monomer2);
            }
        }

        return interacting_monomers;
    }

    eneT Energy::calc_monomer_diff(Monomer& mono1) {
        monomerArrayT monos {m_config.get_monomers()};
        eneT de {0};
        for (size_t i {0}; i != monos.size(); i++) {
            Monomer& mono2 {monos[i].get()};
            if (mono1.get_index() == mono2.get_index()) {
                continue;
            }
            eneT ene1 {calc_monomer_pair_energy(mono1, CoorSet::current, mono2,
                    CoorSet::current)};
            eneT ene2 {calc_monomer_pair_energy(mono1, CoorSet::trial, mono2,
                    CoorSet::current)};
            de += ene2 - ene1;
            }

        return de;
    }

    bool Energy::particles_interacting(Particle& particle1, int conformer1,
            CoorSet coorset1, Particle& particle2, int conformer2,
            CoorSet coorset2) {

        distT dist {m_config.calc_dist(particle1, coorset1, particle2,
                coorset2)};
        pair<int, int> key {particle1.get_type(), particle2.get_type()};
        PairPotential* pot;
        if (conformer1 == conformer2) {
            pot = &(m_same_pair_to_pot.at(key).get());
        }
        else {
            pot = &(m_different_pair_to_pot.at(key).get());
        }
        bool interacting {pot->particles_interacting(dist)};

        return interacting;
    }

    eneT Energy::calc_particle_pair_energy(Particle& particle1, int conformer1,
            CoorSet coorset1, Particle& particle2, int conformer2,
            CoorSet coorset2) {

        vecT diff {m_config.calc_interparticle_vector(particle2, coorset2,
                particle1, coorset1)};
        distT dist {diff.norm()};
        auto p1_ore {particle1.get_ore(coorset1)};
        auto p2_ore {particle2.get_ore(coorset2)};
        pair<int, int> key {particle1.get_type(), particle2.get_type()};

        PairPotential* pot;
        if (conformer1 == conformer2) {
            pot = &(m_same_pair_to_pot.at(key).get());
        }
        else {
            pot = &(m_different_pair_to_pot.at(key).get());
        }
        eneT ene {pot->calc_energy(dist, diff, p1_ore, p2_ore)};

        return ene;
    }

    void Energy::create_potentials(vector<PotentialData> potentials,
            vector<InteractionData> same_conformers_interactions,
            vector<InteractionData> different_conformers_interactions) {

        for (auto p_data: potentials) {
            PairPotential* pot;
            if (p_data.form == "Zero") {
                pot = new ZeroPotential {};
            }
            if (p_data.form == "HardSphere") {
                pot = new  HardSpherePotential {p_data.sigh};
            }
            if (p_data.form == "SquareWell") {
                pot = new  SquareWellPotential {p_data.eps, p_data.rcut};
            }
            if (p_data.form == "HarmonicWell") {
                pot = new HarmonicWellPotential {p_data.eps, p_data.rcut};
            }
            if (p_data.form == "AngularHarmonicWell") {
                pot = new AngularHarmonicWellPotential {p_data.eps, p_data.rcut, p_data.siga1};
            }
            else if (p_data.form == "ShiftedLJ") {
                pot = new ShiftedLJPotential {p_data.eps, p_data.sigl,
                        p_data.rcut};
            }
            else if (p_data.form == "Patchy") {
                pot = new PatchyPotential {p_data.eps, p_data.sigl, p_data.rcut,
                        p_data.siga1, p_data.siga2};
            }
            else if (p_data.form == "OrientedPatchy") {
                pot = new OrientedPatchyPotential {p_data.eps, p_data.sigl,
                        p_data.rcut, p_data.siga1, p_data.siga2,
                        p_data.sigt};
            }
            else if (p_data.form == "DoubleOrientedPatchy") {
                pot = new DoubleOrientedPatchyPotential {p_data.eps,
                        p_data.sigl, p_data.rcut, p_data.siga1, p_data.siga2,
                        p_data.sigt};
            }
            m_potentials.emplace_back(pot);
        }

        for (auto i_data: same_conformers_interactions) {
            for (auto p_pair: i_data.particle_pairs) {
                PairPotential& pot {*m_potentials[i_data.potential_index]};
                m_same_pair_to_pot.emplace(p_pair, pot);
                pair<int, int> reversed_p_pair {p_pair.second, p_pair.first};
                m_same_pair_to_pot.emplace(reversed_p_pair, pot);
            }
        }

        for (auto i_data: different_conformers_interactions) {
            for (auto p_pair: i_data.particle_pairs) {
                PairPotential& pot {*m_potentials[i_data.potential_index]};
                m_different_pair_to_pot.emplace(p_pair, pot);
                pair<int, int> reversed_p_pair {p_pair.second, p_pair.first};
                m_different_pair_to_pot.emplace(reversed_p_pair, pot);
            }
        }
    }

}
