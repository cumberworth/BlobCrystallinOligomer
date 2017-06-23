// energy.h

#ifndef ENERGY_H
#define ENERGY_H

#include <memory>
#include <unordered_map>
#include <utility>
#include <vector>

#include "BlobCrystallinOligomer/config.h"
#include "BlobCrystallinOligomer/ifile.h"
#include "BlobCrystallinOligomer/hash.h"
#include "BlobCrystallinOligomer/monomer.h"
#include "BlobCrystallinOligomer/param.h"
#include "BlobCrystallinOligomer/particle.h"
#include "BlobCrystallinOligomer/potential.h"
#include "BlobCrystallinOligomer/shared_types.h"

namespace energy {

    using config::Config;
    using config::monomerArrayT;
    using ifile::InteractionData;
    using ifile::PotentialData;
    using monomer::Monomer;
    using param::InputParams;
    using particle::Particle;
    using potential::PairPotential;
    using shared_types::CoorSet;
    using shared_types::eneT;
    using std::pair;
    using std::reference_wrapper;
    using std::unique_ptr;
    using std::unordered_map;
    using std::vector;

    class Energy {
        public:
            Energy(Config& conf, InputParams params);
            Energy(Config& conf, vector<PotentialData>,
                    vector<InteractionData>);

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
            /*  Check if monomers within range to have non-zero pair potential */

            monomerArrayT get_interacting_monomers(Monomer& monomer1,
                    CoorSet coorset1);
            /*  Create list of monomers interacting with given */

            bool particles_interacting(
                    Particle& particle1,
                    CoorSet coorset1,
                    Particle& particle2,
                    CoorSet coorset2);
            /*  Check if particles within range to have non-zero pair potential */

            eneT calc_particle_pair_energy(
                    Particle& particle1,
                    CoorSet coorset1,
                    Particle& particle2,
                    CoorSet coorset2);
            /*  Calculate pair energy between two particles.*/

        private:
            Config& m_config;
            vector<unique_ptr<PairPotential>> m_potentials;
            //unordered_map<pair<int, int>, reference_wrapper<PairPotential>> m_pair_to_pot;
            unordered_map<pair<int, int>, reference_wrapper<PairPotential>> m_pair_to_pot;

            void create_potentials(vector<PotentialData> potentials,
                    vector<InteractionData> interactions);
    };
}

#endif // ENERGY_H
