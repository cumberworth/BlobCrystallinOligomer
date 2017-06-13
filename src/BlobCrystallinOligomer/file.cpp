// file.cpp

#include <fstream>
#include <string>

#include "Json/json.hpp"

#include "BlobCrystallinOligomer/file.h"

namespace file {

    using nlohmann::json;
    using std::ifstream;
    using std::string;

    InputConfigFile::InputConfigFile(string filename) {
        ifstream file {filename};
        file >> m_config_json;
        parse_json();
    }

    vector<MonomerData> InputConfigFile::get_monomers() {
        return m_monomers;
    }

    void InputConfigFile::parse_json() {
        for (auto json_monomer: m_config_json["cgmonomer"]["config"]) {
            int monomer_index {json_monomer["index"]};
            vector<ParticleData> particles {};
            for (auto json_particle: json_monomer["particles"]) {
                int p_index {json_particle["index"]};
                string p_type {json_particle["index"].get<string>()};
                string domain {json_particle["domain"].get<string>()};
                vecT pos {json_particle["pos"]};
                vecT patch_norm;
                vecT patch_orient;
                if (p_type == "SimpleParticle") {
                }
                if (p_type == "PatchyParticle") {
                    patch_norm = json_particle["patch_norm"];
                }
                if (p_type == "OrientedPatchyParticle") {
                    patch_orient = json_particle["patch_orient"];
                }
                ParticleData p_data {p_index, domain, p_type, pos, patch_norm,
                    patch_orient};
                particles.push_back(p_data);
            }
            MonomerData m_data {monomer_index, particles};
            m_monomers.push_back(m_data);
        }
    }

    InputEnergyFile::InputEnergyFile(string filename) {
        ifstream file {filename};
        file >> m_energy_json;
        parse_json();
    }

    void InputEnergyFile::parse_json() {
        json json_potentials {m_energy_json["cgmonomer"]["energy"]["potentials"]};
        for (auto json_potential: json_potentials) {
            PotentialData pot_data {};
            pot_data.index = json_potential["index"];
            string pot_type {json_potential["type"].get<string>()};
            if (pot_type == "HardSphere") {
                pot_data.sigh = json_potential["parameters"]["sigh"];
            }
            if (pot_type == "ShiftedLJ") {
                pot_data.sigl = json_potential["parameters"]["sigl"];
                pot_data.eps = json_potential["parameters"]["eps"];
            }
            if (pot_type == "Patchy") {
                pot_data.siga1 = json_potential["parameters"]["siga1"];
                pot_data.siga2 = json_potential["parameters"]["siga2"];
            }
            if (pot_type == "OrientedPatchy") {
                pot_data.sigt1 = json_potential["parameters"]["sigt1"];
                pot_data.sigt2 = json_potential["parameters"]["sigt2"];
            }
            m_potentials.push_back(pot_data);
        }
        json json_inters {m_energy_json["cgmonomer"]["energy"]["interactions"]};
        for (auto json_inter: json_inters) {
            vector<vector<int>> raw_pairs {json_inter["pairs"]};
            vector<pair<int, int>> particle_pairs {};
            for (auto p_pair: raw_pairs) {
                particle_pairs.push_back({p_pair[0], p_pair[1]});
            }
            int potential {json_inter["potential"]};
            InteractionData i_data {particle_pairs, potential};
            m_interactions.push_back(i_data);
        }
    }
}
