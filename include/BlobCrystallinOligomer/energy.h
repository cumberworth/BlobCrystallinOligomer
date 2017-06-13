// energy.h

#ifndef ENERGY_H
#define ENERGY_H

#include <vector>

#include "BlobCrystallinOligomer/config.h"
#include "BlobCrystallinOligomer/file.h"
#include "BlobCrystallinOligomer/monomer.h"
#include "BlobCrystallinOligomer/param.h"
#include "BlobCrystallinOligomer/particle.h"
#include "BlobCrystallinOligomer/shared_types.h"

namespace energy {

    using config::Config;
    using file::InteractionData;
    using file::PotentialData;
    using monomer::Monomer;
    using param::InputParams;
    using particle::Particle;
    using shared_types::CoorSet;
    using shared_types::eneT;
    using std::vector;

    class Energy {
        public:
            Energy(Config& conf, InputParams params);

            eneT calc_monomer_pair_energy(
                    Monomer& monomer1,
                    CoorSet coorset1,
                    Monomer& monomer2,
                    CoorSet coorset2);
            /*  Calculate pair energy between two monomers.*/

            bool monomers_interacting(
                    Monomer& monomer1,
                    CoorSet coorset1,
                    Monomer& monomer2,
                    CoorSet coorset2);
            /* Check if monomers within range to have non-zero pair potential */

            eneT calc_particle_pair_energy(
                    Particle& particle1,
                    CoorSet coorset1,
                    Particle& particle2,
                    CoorSet coorset2);
            /*  Calculate pair energy between two particles.*/

        private:
            Config& m_config;
            // map of particle type pairs to potentials

            void create_potentials(vector<PotentialData> potentials,
                    vector<InteractionData> interactions);
    };
}

#endif // ENERGY_H
