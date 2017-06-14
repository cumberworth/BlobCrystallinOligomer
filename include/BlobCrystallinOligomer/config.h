// config.h

#ifndef CONFIG_H
#define CONFIG_H

#include <memory>
#include <vector>

#include "BlobCrystallinOligomer/file.h"
#include "BlobCrystallinOligomer/monomer.h"
#include "BlobCrystallinOligomer/param.h"
#include "BlobCrystallinOligomer/particle.h"
#include "BlobCrystallinOligomer/shared_types.h"
#include "BlobCrystallinOligomer/space.h"

namespace config {

    using file::InputConfigFile;
    using file::MonomerData;
    using file::ParticleData;
    using monomer::Monomer;
    using param::InputParams;
    using particle::Particle;
    using shared_types::distT;
    using shared_types::CoorSet;
    using shared_types::vecT;
    using space::CuboidPBC;
    using std::reference_wrapper;
    using std::unique_ptr;
    using std::vector;

    typedef vector<reference_wrapper<Monomer>> monomerArrayT;

    class Config {
        public:
            Config(InputParams params);

            Monomer& get_monomer(int monomer_index);
            
            monomerArrayT get_monomers();

            vecT calc_interparticle_vector(
                    Particle& particle1,
                    CoorSet coorset1,
                    Particle& particle2,
                    CoorSet coorset2);
            /*  Calculate vector from particle 1 to particle 2 */

            distT calc_dist(
                    Particle& particle1,
                    CoorSet coorset1,
                    Particle& particle2,
                    CoorSet coorset2);
            /*  Calculate distance between two particles */

        private:
            vector<unique_ptr<Monomer>> m_monomers;
            monomerArrayT m_monomer_refs;
            CuboidPBC m_space;

            void create_monomers(vector<MonomerData>);
    };
}

#endif // CONFIG_H
