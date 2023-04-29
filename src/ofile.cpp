// ofile.cpp

#include "BlobCrystallinOligomer/ofile.h"
#include "BlobCrystallinOligomer/monomer.h"
#include "BlobCrystallinOligomer/particle.h"

namespace ofile {

using monomer::Monomer;
using nlohmann::json;
using particle::Particle;
using shared_types::CoorSet;
using shared_types::distT;
using std::ifstream;
using std::string;

OutputFile::OutputFile() {}

OutputFile::OutputFile(string filename):
        m_file {filename}, m_filename {filename} {}

void OutputFile::close() { m_file.close(); }

VSFOutputFile::VSFOutputFile() {}

VSFOutputFile::VSFOutputFile(string filename, Config& conf):
        OutputFile {filename} {

    write_structure(conf);
}

void VSFOutputFile::write_structure(Config& conf) {
    // This is a bit of a hack for the radius
    int i {0};
    for (auto m: conf.get_monomers()) {
        for (auto p: m.get().get_particles()) {
            m_file << "atom " << i << " ";
            m_file << "type " << p.get().get_type() << " ";
            m_file << "resid " << m.get().get_index() << " ";
            m_file << "radius " << conf.get_radius() << "\n";
            i++;
        }
    }
    m_file << "\n";
    distT x {conf.get_box_len()};
    m_file << "pbc " << x << " " << x << " " << x << "\n";
    m_file << "\n";
}

VCFOutputFile::VCFOutputFile() {}

VCFOutputFile::VCFOutputFile(string filename): OutputFile {filename} {}

void VCFOutputFile::write_step(Config& conf, stepT) {
    m_file << "t"
           << "\n";
    for (Monomer& mono: conf.get_monomers()) {
        for (Particle& part: mono.get_particles()) {
            auto pos {part.get_pos(CoorSet::current)};
            m_file << pos[0] << " " << pos[1] << " " << pos[2] << "\n";
        }
    }
    m_file << "\n";
}

void VCFOutputFile::open_write_step_close(Config& conf, stepT step) {
    m_file.open(m_filename);
    write_step(conf, step);
    m_file.close();
}

VTFOutputFile::VTFOutputFile(string filename, Config& conf):
        OutputFile {filename} {

    write_structure(conf);
}

void VTFOutputFile::open_write_close(Config& conf, stepT step) {
    m_file.open(m_filename);
    write_structure(conf);
    write_step(conf, step);
    m_file.close();
}

PatchOutputFile::PatchOutputFile(string filename): OutputFile {filename} {}

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
            for (int i {0}; i != 3; i++) {
                m_file << ore.patch_orient2[i] << " ";
            }
        }
    }
    m_file << "\n";
}

void PatchOutputFile::open_write_step_close(Config& conf) {
    m_file.open(m_filename);
    write_step(conf);
    m_file.close();
}
} // namespace ofile
