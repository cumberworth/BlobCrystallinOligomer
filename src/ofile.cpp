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

    OutputConfigsFile::OutputConfigsFile(string filename, Config& conf):
            m_file {filename} {

         write_structure(conf);
    }

    void OutputConfigsFile::write_structure(Config& conf) {
        // This is a bit of a hack for the radius
        m_file << "0" << ":" << conf.get_num_particles() << " ";
        m_file << "radius " << conf.get_radius() << "\n";
        m_file << "\n";
        distT x {conf.get_box_len()};
        m_file << "pbc " << x << " " << x << " " << x << "\n";
        m_file << "\n";
    }

    void OutputConfigsFile::write_step(Config& conf, stepT) {
        m_file << "t" << "\n";
        for (Monomer& mono: conf.get_monomers()) {
            for (Particle& part: mono.get_particles()) {
                auto pos {part.get_pos(CoorSet::current)};
                m_file << pos[0] << " " << pos[1] << " " << pos[2] << "\n";
            }
        }
        m_file << "\n";
    }
}
