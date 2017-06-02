// energy.h

#ifndef ENERGY_H
#define ENERGY_H

#include "shared_types.h"
#include "params.h"

namespace energy {

    using shared_types::eneT;
    using params::InputParams;

    class Energy {
        public:
            Energy(InputParams params);
            eneT calc_monomer_pair();
            eneT calc_sphere_pair();

        private:
    };
}

#endif // ENERGY_H
