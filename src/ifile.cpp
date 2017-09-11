// ifile.cpp

#include <fstream>
#include <string>

#include "Json/json.hpp"

#include "BlobCrystallinOligomer/ifile.h"

namespace ifile {

    using shared_types::CoorSet;
    using std::ifstream;

    vecT json2vec(json jvec) {
        vector<distT> vec;
        for (auto comp: jvec) {
            vec.push_back(comp);
        }

        return {vec[0], vec[1], vec[2]};
    }

    InputConfigFile::InputConfigFile(string filename) {
        ifstream file {filename};
        file >> m_config_json;
        parse_json();
    }

    vector<MonomerData> InputConfigFile::get_monomers() {
        return m_monomers;
    }

    distT InputConfigFile::get_box_len() {
        return m_box_len;
    }

    distT InputConfigFile::get_radius() {
        return m_radius;
    }

    void InputConfigFile::parse_json() {
        m_box_len = m_config_json["cgmonomer"]["box_len"];
        m_radius = m_config_json["cgmonomer"]["radius"];
        for (auto json_monomer: m_config_json["cgmonomer"]["config"]) {
            int monomer_index {json_monomer["index"]};
            vector<ParticleData> particles {};
            for (auto json_particle: json_monomer["particles"]) {
                int p_index {json_particle["index"]};
                string p_form {json_particle["form"].get<string>()};
                int p_type {json_particle["type"]};
                string domain {json_particle["domain"].get<string>()};
                vecT pos {json2vec(json_particle["pos"])};
                vecT patch_norm;
                vecT patch_orient;
                if (p_form == "PatchyParticle") {
                    patch_norm = json2vec(json_particle["patch_norm"]);
                }
                if (p_form == "OrientedPatchyParticle") {
                    patch_orient = json2vec(json_particle["patch_orient"]);
                }
                ParticleData p_data {p_index, domain, p_form, p_type, pos, patch_norm,
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

    vector<PotentialData> InputEnergyFile::get_potentials() {
        return m_potentials;
    }

    vector<InteractionData> InputEnergyFile::get_interactions() {
        return m_interactions;
    }

    void InputEnergyFile::parse_json() {
        json json_potentials {m_energy_json["cgmonomer"]["energy"]["potentials"]};
        for (auto json_potential: json_potentials) {
            PotentialData pot_data {};
            pot_data.index = json_potential["index"];
            string pot_form {json_potential["form"].get<string>()};
            pot_data.form = pot_form;
            if (pot_form == "HardSphere") {
                pot_data.sigh = json_potential["parameters"]["sigh"];
            }
            if (pot_form == "ShiftedLJ") {
                pot_data.sigl = json_potential["parameters"]["sigl"];
                pot_data.eps = json_potential["parameters"]["eps"];
            }
            if (pot_form == "Patchy") {
                pot_data.siga1 = json_potential["parameters"]["siga1"];
                pot_data.siga2 = json_potential["parameters"]["siga2"];
            }
            if (pot_form == "OrientedPatchy") {
                pot_data.sigt = json_potential["parameters"]["sigt"];
            }
            m_potentials.push_back(pot_data);
        }
        json json_inters {m_energy_json["cgmonomer"]["energy"]["interactions"]};
        for (auto json_inter: json_inters) {
            auto raw_pairs {json_inter["pairs"]};
            vector<pair<int, int>> particle_pairs {};
            for (auto ipair: raw_pairs) {
                particle_pairs.emplace_back(ipair[0], ipair[1]);
            }
            int potential {json_inter["potential"]};
            InteractionData i_data {particle_pairs, potential};
            m_interactions.push_back(i_data);
        }
    }
}
