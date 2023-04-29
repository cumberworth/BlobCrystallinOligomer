// ifile.cpp

#include <fstream>
#include <string>

#include "Json/json.hpp"

#include "BlobCrystallinOligomer/ifile.h"

namespace ifile {

using std::cout;
using std::ifstream;

using shared_types::CoorSet;
using shared_types::InputError;

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

vector<MonomerData> InputConfigFile::get_monomers() { return m_monomers; }

distT InputConfigFile::get_box_len() { return m_box_len; }

distT InputConfigFile::get_radius() { return m_radius; }

void InputConfigFile::parse_json() {
    m_box_len = m_config_json["cgmonomer"]["box_len"];
    m_radius = m_config_json["cgmonomer"]["radius"];
    for (auto json_monomer: m_config_json["cgmonomer"]["config"]) {
        int monomer_index {json_monomer["index"]};
        int conformer {json_monomer["conformer"]};
        vector<ParticleData> particles {};
        for (auto json_particle: json_monomer["particles"]) {
            int p_index {json_particle["index"]};
            string p_form {json_particle["form"].get<string>()};
            int p_type {json_particle["type"]};
            string domain {json_particle["domain"].get<string>()};
            vecT pos {json2vec(json_particle["pos"])};
            vecT patch_norm;
            vecT patch_orient;
            vecT patch_orient2;
            if (p_form == "PatchyParticle" or
                p_form == "OrientedPatchyParticle" or
                p_form == "DoubleOrientedPatchyParticle" or
                p_form == "AngularHarmonicWellPotential") {
                patch_norm = json2vec(json_particle["patch_norm"]);
            }
            if (p_form == "OrientedPatchyParticle" or
                p_form == "DoubleOrientedPatchyParticle") {
                patch_orient = json2vec(json_particle["patch_orient"]);
            }
            if (p_form == "DoubleOrientedPatchyParticle") {
                patch_orient2 = json2vec(json_particle["patch_orient2"]);
            }
            ParticleData p_data {
                    p_index,
                    domain,
                    p_form,
                    p_type,
                    pos,
                    patch_norm,
                    patch_orient,
                    patch_orient2};
            particles.push_back(p_data);
        }
        MonomerData m_data {monomer_index, conformer, particles};
        m_monomers.push_back(m_data);
    }
}

InputEnergyFile::InputEnergyFile(string filename) {
    ifstream file {filename};
    file >> m_energy_json;
    parse_json();
}

vector<PotentialData> InputEnergyFile::get_potentials() { return m_potentials; }

vector<InteractionData> InputEnergyFile::get_same_conformers_interactions() {
    return m_same_conformers_interactions;
}

vector<InteractionData> InputEnergyFile::
        get_different_conformers_interactions() {
    return m_different_conformers_interactions;
}

void InputEnergyFile::parse_json() {
    json json_potentials {m_energy_json["cgmonomer"]["energy"]["potentials"]};
    //        for (auto json_potential: json_potentials) {
    for (auto json_potential: json_potentials[0]) {
        PotentialData pot_data {};
        pot_data.index = json_potential["index"];
        string pot_form {json_potential["form"].get<string>()};
        pot_data.form = pot_form;
        if (pot_form == "Zero") {}
        else if (pot_form == "HardSphere") {
            pot_data.sigh = json_potential["parameters"]["sigh"];
        }
        else if (pot_form == "SquareWell") {
            pot_data.eps = json_potential["parameters"]["eps"];
            pot_data.rcut = json_potential["parameters"]["rcut"];
        }
        else if (pot_form == "HarmonicWell") {
            pot_data.eps = json_potential["parameters"]["eps"];
            pot_data.rcut = json_potential["parameters"]["rcut"];
        }
        else if (pot_form == "AngularHarmonicWell") {
            pot_data.eps = json_potential["parameters"]["eps"];
            pot_data.rcut = json_potential["parameters"]["rcut"];
            pot_data.siga1 = json_potential["parameters"]["siga1"];
        }
        else if (pot_form == "ShiftedLJ") {
            pot_data.sigl = json_potential["parameters"]["sigl"];
            pot_data.eps = json_potential["parameters"]["eps"];
            pot_data.rcut = json_potential["parameters"]["rcut"];
        }
        else if (pot_form == "Patchy") {
            pot_data.sigl = json_potential["parameters"]["sigl"];
            pot_data.eps = json_potential["parameters"]["eps"];
            pot_data.rcut = json_potential["parameters"]["rcut"];
            pot_data.siga1 = json_potential["parameters"]["siga1"];
            pot_data.siga2 = json_potential["parameters"]["siga2"];
        }
        else if (
                pot_form == "OrientedPatchy" or
                pot_form == "DoubleOrientedPatchy") {
            pot_data.sigl = json_potential["parameters"]["sigl"];
            pot_data.eps = json_potential["parameters"]["eps"];
            pot_data.rcut = json_potential["parameters"]["rcut"];
            pot_data.siga1 = json_potential["parameters"]["siga1"];
            pot_data.siga2 = json_potential["parameters"]["siga2"];
            pot_data.sigt = json_potential["parameters"]["sigt"];
        }
        else {
            cout << "No such potential\n";
            throw InputError {};
        }
        m_potentials.push_back(pot_data);
    }
    json json_inters {m_energy_json["cgmonomer"]["energy"]["interactions"]};
    //        for (auto json_inter: json_inters) {
    for (auto json_inter: json_inters[0]) {
        auto raw_pairs {json_inter["pairs"]};
        vector<pair<int, int>> particle_pairs {};
        for (auto ipair: raw_pairs[0]) {
            particle_pairs.emplace_back(ipair[0], ipair[1]);
        }
        int potential {json_inter["potential"]};
        string conformers {json_inter["conformers"].get<string>()};
        InteractionData i_data {particle_pairs, potential};
        if (conformers == "any") {
            m_same_conformers_interactions.push_back(i_data);
            m_different_conformers_interactions.push_back(i_data);
        }
        else if (conformers == "same") {
            m_same_conformers_interactions.push_back(i_data);
        }
        else if (conformers == "different") {
            m_different_conformers_interactions.push_back(i_data);
        }
    }
}
} // namespace ifile
