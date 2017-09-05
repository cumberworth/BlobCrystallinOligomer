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

    /** JSON file format for single configuration output */
    class OutputConfigFile {
    };

    /** VTF/Tcl file format for "trajectory" configuration output
      *
      * The VTF format allows the particle positions to be read by VMD.
      * However, it cannot read patch vectors, so these are written to a Tcl
      * script.
      */
    class OutputConfigsFile {
        public:
            OutputConfigsFile(string filename, Config& conf);

            /** Write system properties and topology.
              *
              * For now this is just the box size and the particle radii
              */
            void write_structure(Config& conf);

            /** Write configuration at current step */
            void write_step(Config& config, stepT step);

        private:
            std::ofstream m_file;
    };
}

#endif // OFILE_H
