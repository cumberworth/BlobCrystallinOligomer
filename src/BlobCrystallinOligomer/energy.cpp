// energy.cpp

#include <vector>

#include "BlobCrystallinOligomer/config.h"
#include "BlobCrystallinOligomer/energy.h"
#include "BlobCrystallinOligomer/param.h"

namespace energy {

    using config::Config;
    using file::InputEnergyFile;
    using file::InteractionData;
    using file::PotentialData;
    using param::InputParams;
    using std::vector;

    Energy::Energy(Config& conf, InputParams params):
            m_config {conf} {

        InputEnergyFile energy_file {params.m_energy_filename};
        vector<PotentialData> potentials {energy_file.get_potentials()};
        vector<InteractionData> interactions {energy_file.get_interactions()};
        create_potentials(potentials, interactions);
    }

    void Energy::create_potentials(vector<PotentialData> potentials,
            vector<InteractionData> interactions) {
    }
}
