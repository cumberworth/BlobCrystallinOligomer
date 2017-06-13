// file.h

#ifndef FILE_H
#define FILE_H

#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "BlobCrystallinOligomer/shared_types.h"
#include "Json/json.hpp"

namespace file {

    using shared_types::vecT;
    using std::pair;
    using std::unordered_map;
    using std::string;
    using std::vector;

    typedef nlohmann::json json;

    struct ParticleData {
        int index;
        string domain;
        string type;
        vecT pos;
        vecT patch_norm;
        vecT patch_orient;
    };

    struct MonomerData {
        int index;
        vector<ParticleData> particles;
    };

    class InputConfigFile {
        /*  JSON file format for input topology and configuration */
        public:
            InputConfigFile(string filename);
            vector<MonomerData> get_monomers();

        private:
            json m_config_json;
            vector<MonomerData> m_monomers;

            void parse_json();

    };

    struct PotentialData {
        string type;
        int index;
        double sigh;
        double sigl;
        double siga1;
        double siga2;
        double sigt1;
        double sigt2;
        double eps;
    };

    struct InteractionData {
        vector<pair<int, int>> particle_pairs;
        int potential_label;
    };

    class InputEnergyFile {
        /*  JSON file format for input energy specifications */
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

    class OutputConfigFile {
    };
}

#endif // FILE_H
