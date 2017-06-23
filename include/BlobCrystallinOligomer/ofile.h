// ofile.h

#ifndef OFILE_H
#define OFILE_H

#include <fstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "BlobCrystallinOligomer/config.h"
#include "BlobCrystallinOligomer/shared_types.h"
#include "Json/json.hpp"

namespace ofile {

    using config::Config;
    using shared_types::distT;
    using shared_types::stepT;
    using std::unordered_map;
    using std::string;

    typedef nlohmann::json json;

    class OutputConfigFile {
        /*  JSON file format for single configuration output */
    };

    class OutputConfigsFile {
        /*  VTF/Tcl file format for "trajectory" configuration output
       
            The VTF format allows the particle positions to be read by VMD.
            However, it cannot patch vectors, so these are written to a Tcl
            script.
       */
        public:
            OutputConfigsFile(string filename, Config& conf);

            void write_structure(Config& conf);
            /*  Write system properties and topology.

                For now this is just the box size and the particle radii
            */

            void write_timestep(Config& config, stepT step);

        private:
            std::ofstream m_file;
    };
}

#endif // OFILE_H
