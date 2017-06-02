// energy.h

#ifndef ENERGY_H
#define ENERGY_H

#include "BlobCrystallinOligomer/shared_types.h"
#include "BlobCrystallinOligomer/params.h"
#include "BlobCrystallinOligomer/particle.h"
#include "BlobCrystallinOligomer/monomer.h"
#include "BlobCrystallinOligomer/config.h"

namespace energy {

    using shared_types::eneT;
    using shared_types::CoorSet;
    using params::InputParams;
    using particle::Particle;
    using monomer::Monomer;
    using config::Config;

    template<typename configT>
    class Energy {
        public:
            Energy(configT& config, InputParams params);
            Energy(InputParams params);
            eneT calc_monomer_pair_energy(
                    Monomer& monomer1,
                    CoorSet coorset1,
                    Monomer& monomer2,
                    CoorSet coorset2);
            eneT calc_sphere_pair_energy(
                    Particle& particle1,
                    CoorSet coorset1,
                    Particle& particle2,
                    CoorSet coorset2);

        private:
            configT& m_config;
            // map of particle type pairs to potentials
    };
}

#endif // ENERGY_H
