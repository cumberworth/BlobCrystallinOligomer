// system.h

#ifndef CONFIG_H
#define CONFIG_H

#include <vector>

#include "BlobCrystallinOligomer/shared_types.h"
#include "BlobCrystallinOligomer/params.h"
#include "BlobCrystallinOligomer/file.h"
#include "BlobCrystallinOligomer/space.h"
#include "BlobCrystallinOligomer/particle.h"
#include "BlobCrystallinOligomer/monomer.h"

namespace config {

    using std::vector;
    using std::reference_wrapper;
    using shared_types::distT;
    using shared_types::CoorSet;
    using params::InputParams;
    using file::InputConfigFile;
    using particle::Particle;
    using monomer::Monomer;

    typedef vector<reference_wrapper<Monomer>> monomerArrayT;

    template<typename spaceT>
    class Config {
        public:
            Config(InputConfigFile& input_file, InputParams params);
            Monomer& get_monomer(int monomer_index);

            bool monomers_interacting(
                    Monomer& monomer1,
                    CoorSet coorset1,
                    Monomer& monomer2,
                    CoorSet coorset2);
            // Check if monomers within range to have non-zero pair potential

            // Consider replacing with calc_dist and calc_trial dist, as will
            // only care about current case or the one particle trial case
            distT calc_dist(
                    Particle& particle1,
                    CoorSet coorset1,
                    Particle& particle2,
                    CoorSet coorset2);
            // Calculate distance between two particles

        private:
            monomerArrayT m_monomers;
            spaceT m_space;

    };
}

#endif // CONFIG_H
