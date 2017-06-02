// system.h

#ifndef CONFIG_H
#define CONFIG_H

#include <vector>

#include "monomer.h"

namespace config {

    using std::vector;
    using std::reference_wrapper;

    using monomer::Monomer;

    typedef vector<reference_wrapper<Monomer>> monomerArrayT;

    class Config {
        public:
             monomerArrayT calc_interacting_monomers(Monomer m1);

        private:
             monomerArrayT m_monomers;

    };
}

#endif // CONFIG_H
