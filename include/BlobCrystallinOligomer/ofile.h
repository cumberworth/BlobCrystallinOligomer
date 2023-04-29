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
using std::string;
using std::unordered_map;

typedef nlohmann::json json;

/** Base class for output files */
class OutputFile {
  public:
    OutputFile();
    OutputFile(string filename);
    virtual ~OutputFile() {}
    void close();

  protected:
    std::ofstream m_file;
    string m_filename;
};

/** VSF file format for topology output */
class VSFOutputFile: virtual public OutputFile {
  public:
    VSFOutputFile();
    VSFOutputFile(string filename, Config& conf);

    /** Write system properties and topology.
     *
     * For now this is just the box size and the particle radii
     */
    void write_structure(Config& conf);
};

/** VCF file format for position frame output */
class VCFOutputFile: virtual public OutputFile {
  public:
    VCFOutputFile();
    VCFOutputFile(string filename);

    /** Write configuration at current step */
    void write_step(Config& config, stepT step);
    void open_write_step_close(Config& config, stepT step);
};

/** VTF file format for topology and position frame output */
class VTFOutputFile:
        public VSFOutputFile,
        public VCFOutputFile,
        virtual public OutputFile {
  public:
    VTFOutputFile(string filename, Config& conf);
    void open_write_close(Config& config, stepT step);
};

/** Simple output format for patch vectors
 *
 * All patch vectors for a given frame are written on one line. All vectors
 * for a given particle are written consecutively.
 */
class PatchOutputFile: virtual public OutputFile {
  public:
    PatchOutputFile(string filename);
    void write_step(Config& conf);
    void open_write_step_close(Config& conf);
};
} // namespace ofile

#endif // OFILE_H
