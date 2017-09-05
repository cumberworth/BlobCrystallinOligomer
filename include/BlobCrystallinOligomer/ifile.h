// ifile.h

#ifndef IFILE_H
#define IFILE_H

#include <fstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "BlobCrystallinOligomer/shared_types.h"
#include "Json/json.hpp"

namespace ifile {

    using nlohmann::json;
    using shared_types::vecT;
    using shared_types::distT;
    using shared_types::stepT;
    using std::pair;
    using std::unordered_map;
    using std::string;
    using std::vector;

    /** Convert a 3D json vector to an eigen vector.
      *
      * Probably a better way to do this than having a custom function
      */
    vecT json2vec(json jvec);

    /** For passing particle data to the config class */
    struct ParticleData {
        int index;
        string domain;
        string form;
        int type;
        vecT pos;
        vecT patch_norm;
        vecT patch_orient;
    };

    /** For passing monomer data to the config class */
    struct MonomerData {
        int index;
        vector<ParticleData> particles;
    };

    /** JSON file format for input topology and configuration */
    class InputConfigFile {
        public:
            InputConfigFile(string filename);
            vector<MonomerData> get_monomers();
            distT get_box_len();
            distT get_radius();

        private:
            json m_config_json;
            vector<MonomerData> m_monomers;
            distT m_box_len;
            distT m_radius;

            void parse_json();

    };

    /** For passing potential data to the energy class */
    struct PotentialData {
        string form;
        int index;
        double sigh;
        double sigl;
        double siga1;
        double siga2;
        double sigt;
        double eps;
        double rcut;
    };

    /** For passing interaction pair type data to the energy class */
    struct InteractionData {
        vector<pair<int, int>> particle_pairs;
        int potential_index;
    };

    /** JSON file format for input energy specifications */
    class InputEnergyFile {
        public:
            InputEnergyFile(string filename);
            vector<PotentialData> get_potentials();
            vector<InteractionData> get_interactions();

        private:
            json m_energy_json;
            vector<PotentialData> m_potentials {};
            vector<InteractionData> m_interactions {};

            void parse_json();

    };
}

#endif // IFILE_H
