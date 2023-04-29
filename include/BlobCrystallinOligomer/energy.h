// energy.h

#ifndef ENERGY_H
#define ENERGY_H

#include <memory>
#include <unordered_map>
#include <utility>
#include <vector>

#include "BlobCrystallinOligomer/config.h"
#include "BlobCrystallinOligomer/hash.h"
#include "BlobCrystallinOligomer/ifile.h"
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
using shared_types::distT;
using shared_types::eneT;
using std::pair;
using std::reference_wrapper;
using std::unique_ptr;
using std::unordered_map;
using std::vector;

/** System energy
 *
 * Contains all potentials present in system and maps from pairs of
 * particles to their interaction potential type. Responsible for
 * instantiating the potentials.
 */
class Energy {
  public:
    Energy(Config& conf, InputParams& params);
    Energy(Config& conf, vector<PotentialData>, vector<InteractionData>);

    /** Calculate total system energy */
    eneT calc_total_energy();

    /** Calculate pair energy between two monomers */
    eneT calc_monomer_pair_energy(
            Monomer& monomer1,
            CoorSet coorset1,
            Monomer& monomer2,
            CoorSet coorset2);

    /** Check if monomers within range to have non-zero pair potential */
    bool monomers_interacting(
            Monomer& monomer1,
            CoorSet coorset1,
            Monomer& monomer2,
            CoorSet coorset2);

    /** Create list of monomers interacting with given monomer */
    monomerArrayT get_interacting_monomers(Monomer& monomer1, CoorSet coorset1);

    /** Calculate energy difference between current and trial */
    eneT calc_monomer_diff(Monomer& monomer);

    /** Check if particles within range to have non-zero pair potential */
    bool particles_interacting(
            Particle& particle1,
            int conformer1,
            CoorSet coorset1,
            Particle& particle2,
            int conformer2,
            CoorSet coorset2);

    /** Calculate pair energy between two particles.*/
    eneT calc_particle_pair_energy(
            Particle& particle1,
            int conformer1,
            CoorSet coorset1,
            Particle& particle2,
            int conformer2,
            CoorSet coorset2);

  private:
    Config& m_config;
    vector<unique_ptr<PairPotential>> m_potentials;
    unordered_map<pair<int, int>, reference_wrapper<PairPotential>>
            m_same_pair_to_pot;
    unordered_map<pair<int, int>, reference_wrapper<PairPotential>>
            m_different_pair_to_pot;
    distT m_max_cutoff;

    void create_potentials(
            vector<PotentialData> potentials,
            vector<InteractionData> same_conformers_interactions,
            vector<InteractionData> different_conformers_interactions);
    bool monomers_in_range(
            Monomer& monomer1,
            CoorSet coorset1,
            Monomer& monomer2,
            CoorSet coorset2);
};
} // namespace energy

#endif // ENERGY_H
