// config.cpp

#include <iostream>
#include <memory>
#include <utility>

#include "BlobCrystallinOligomer/config.h"

namespace config {

using std::cout;
using std::make_unique;
using std::pair;

Config::Config(InputParams& params, RandomGens& random_num):
        m_space_store {new CuboidPBC()},
        m_space {*m_space_store},
        m_random_num {random_num} {

    InputConfigFile config_file {params.m_config_filename};
    m_box_len = config_file.get_box_len();
    m_space.set_len(config_file.get_box_len());
    m_radius = config_file.get_radius();
    vector<MonomerData> monomers {config_file.get_monomers()};
    create_monomers(monomers);

    // Create monomer reference array
    for (auto& m: m_monomers) {
        m_monomer_refs.emplace_back(*m);
    }
}

Config::Config(
        vector<MonomerData> monomers,
        RandomGens& random_num,
        distT box_len,
        distT radius):
        m_space_store {new CuboidPBC()},
        m_space {*m_space_store},
        m_random_num {random_num},
        m_box_len {box_len},
        m_radius {radius} {

    m_space.set_len(m_box_len);
    create_monomers(monomers);

    // Create monomer reference array
    for (auto& m: m_monomers) {
        m_monomer_refs.emplace_back(*m);
    }
}

Monomer& Config::get_monomer(int monomer_index) {
    return *m_monomers[monomer_index];
}

Monomer& Config::get_random_monomer() {
    int m_i {m_random_num.uniform_int(0, m_monomers.size() - 1)};
    return *m_monomers[m_i];
}

monomerArrayT Config::get_monomers() { return m_monomer_refs; }

int Config::get_num_particles() {
    int num_parts {0};
    for (Monomer& mono: m_monomer_refs) {
        num_parts += mono.get_num_particles();
    }

    return num_parts;
}

int Config::get_num_monomers() { return m_monomers.size(); }

distT Config::get_box_len() { return m_box_len; }

distT Config::get_radius() { return m_radius; }

vecT Config::calc_interparticle_vector(
        Particle& particle1,
        CoorSet coorset1,
        Particle& particle2,
        CoorSet coorset2) {

    vecT& pos1 {particle1.get_pos(coorset1)};
    vecT& pos2 {particle2.get_pos(coorset2)};

    return m_space.calc_diff(pos1, pos2);
}

distT Config::calc_dist(
        Particle& particle1,
        CoorSet& coorset1,
        Particle& particle2,
        CoorSet& coorset2) {

    vecT& pos1 {particle1.get_pos(coorset1)};
    vecT& pos2 {particle2.get_pos(coorset2)};

    return m_space.calc_dist(pos1, pos2);
}

distT Config::calc_dist(
        Monomer& monomer1,
        CoorSet& coorset1,
        Monomer& monomer2,
        CoorSet& coorset2) {

    vecT pos1 {monomer1.get_center(coorset1)};
    vecT pos2 {monomer2.get_center(coorset2)};

    return m_space.calc_dist(pos1, pos2);
}

void Config::update_config_positions(vector<vector<double>> positions) {
    int pos_i {0};
    for (auto monomer: m_monomer_refs) {
        auto particles {monomer.get().get_particles()};
        for (auto particle: particles) {
            vecT pos;
            for (int i {0}; i != 3; i++) {
                pos[i] = positions[pos_i][i];
            }
            particle.get().set_pos(pos);
            pos_i++;
        }
    }
}

void Config::create_monomers(vector<MonomerData> monomers) {
    for (auto m_data: monomers) {
        m_monomers.emplace_back(make_unique<Monomer>(m_data, m_space));
    }
}
} // namespace config
