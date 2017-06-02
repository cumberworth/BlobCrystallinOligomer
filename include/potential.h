// potential.h

#ifndef POTENTIAL_H
#define POTENTIAL_H

#include "shared_types.h"

namespace potential {

    using shared_types::distT;
    using shared_types::eneT;

    class PairPotential {
        public:
            virtual eneT calc_energy(?, distT rdist) = 0;

    };

    // Do these contain the parameters, or are they passed them?
    // Am I making one for each monomer, one for each monomer type (so 1 total)?
    class LJPotential: public PairPotential {
        public:
            eneT calc_energy(?, distT rdist);

    };

    class AngularPotential: public PairPotential {
        public:
            eneT calc_energy(?, distT theta);
    };
}

#endif // POTENTIAL_H
