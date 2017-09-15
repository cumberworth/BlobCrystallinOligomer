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

    class OutputFile {
        public:
            OutputFile();
            OutputFile(string filename);
            std::ofstream m_file;
    };

    /** VSF file format for topology output */
    class VSFOutputFile:
            virtual private OutputFile {
        public:
            VSFOutputFile();
            VSFOutputFile(string filename, Config& conf);

            /** Write system properties and topology.
              *
              * For now this is just the box size and the particle radii
              */
            void write_structure(Config& conf);

        private:
            std::ofstream m_file;
    };

    /** VCF file format for position frame output */
    class VCFOutputFile:
            virtual private OutputFile {
        public:
            VCFOutputFile();
            VCFOutputFile(string filename);

            /** Write configuration at current step */
            void write_step(Config& config, stepT step);

        private:
    };

    /** VTF file format for topology and position frame output */
    class VTFOutputFile: 
            public VSFOutputFile,
            public VCFOutputFile,
            virtual private OutputFile {
        public:
            VTFOutputFile(string filename, Config& conf);

        private:
            std::ofstream m_file;
    };

    /** Simple output format for patch vectors
      *
      * All patch vectors for a given frame are written on one line. All vectors
      * for a given particle are written consecutively.
      */
    class PatchOutputFile:
            virtual private OutputFile {
        public:
            PatchOutputFile(string filename);
            void write_step(Config& conf);
    };
}

#endif // OFILE_H
