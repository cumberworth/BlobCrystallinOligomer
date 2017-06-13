// config.h

#ifndef CONFIG_H
#define CONFIG_H

#include <vector>

#include "BlobCrystallinOligomer/shared_types.h"
#include "BlobCrystallinOligomer/param.h"
#include "BlobCrystallinOligomer/file.h"
#include "BlobCrystallinOligomer/space.h"
#include "BlobCrystallinOligomer/particle.h"
#include "BlobCrystallinOligomer/monomer.h"

namespace config {

    using file::InputConfigFile;
    using file::MonomerData;
    using file::ParticleData;
    using monomer::Monomer;
    using param::InputParams;
    using particle::Particle;
    using shared_types::distT;
    using shared_types::CoorSet;
    using space::CuboidPBC;
    using std::vector;
    using std::reference_wrapper;

    typedef vector<reference_wrapper<Monomer>> monomerArrayT;

    class Config {
        public:
            Config(InputParams params);
            void extract_config_from_monomers(vector<MonomerData>);
            /*  From standard config transfer format.
                Must be called to create
                the internal config data structure.
            */

            Monomer& get_monomer(int monomer_index);

            // Consider replacing with calc_dist and calc_trial dist, as will
            // only care about current case or the one particle trial case
            distT calc_dist(
                    Particle& particle1,
                    CoorSet coorset1,
                    Particle& particle2,
                    CoorSet coorset2);
            /*  Calculate distance between two particles */

        private:
            monomerArrayT m_monomers {};
            CuboidPBC m_space;

    };
}

#endif // CONFIG_H
