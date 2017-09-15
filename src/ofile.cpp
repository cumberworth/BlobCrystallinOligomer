// ofile.cpp

#include "BlobCrystallinOligomer/monomer.h"
#include "BlobCrystallinOligomer/ofile.h"
#include "BlobCrystallinOligomer/particle.h"

namespace ofile {

    using monomer::Monomer;
    using nlohmann::json;
    using particle::Particle;
    using shared_types::distT;
    using shared_types::CoorSet;
    using std::ifstream;
    using std::string;

    OutputFile::OutputFile() {
    }

    OutputFile::OutputFile(string filename):
            m_file {filename} {
    }

    VSFOutputFile::VSFOutputFile() {
    }

    VSFOutputFile::VSFOutputFile(string filename, Config& conf):
            OutputFile {filename} {

         write_structure(conf);
    }

    void VSFOutputFile::write_structure(Config& conf) {
        // This is a bit of a hack for the radius
        m_file << "0" << ":" << conf.get_num_particles() - 1 << " ";
        m_file << "radius " << conf.get_radius() << "\n";
        m_file << "\n";
        distT x {conf.get_box_len()};
        m_file << "pbc " << x << " " << x << " " << x << "\n";
        m_file << "\n";
    }

    VCFOutputFile::VCFOutputFile() {
    }

    VCFOutputFile::VCFOutputFile(string filename):
            OutputFile {filename} {
    }

    void VCFOutputFile::write_step(Config& conf, stepT) {
        m_file << "t" << "\n";
        for (Monomer& mono: conf.get_monomers()) {
            for (Particle& part: mono.get_particles()) {
                auto pos {part.get_pos(CoorSet::current)};
                m_file << pos[0] << " " << pos[1] << " " << pos[2] << "\n";
            }
        }
        m_file << "\n";
    }

    VTFOutputFile::VTFOutputFile(string filename, Config& conf):
            OutputFile {filename} {

         write_structure(conf);
    }

    PatchOutputFile::PatchOutputFile(string filename):
            OutputFile {filename} {
    }

    void PatchOutputFile::write_step(Config& conf) {
        for (Monomer& mono: conf.get_monomers()) {
            for (Particle& part: mono.get_particles()) {
                auto ore = part.get_ore(CoorSet::current);
                for (int i {0}; i != 3; i++) {
                    m_file << ore.patch_norm[i] << " ";
                }
                for (int i {0}; i != 3; i++) {
                    m_file << ore.patch_orient[i] << " ";
                }
            }
        }
        m_file << "\n";
    }
}
