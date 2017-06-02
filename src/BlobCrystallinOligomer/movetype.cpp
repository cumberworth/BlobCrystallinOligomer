// movetype.cpp

#include "movetype.h"

namespace movetype {

    bool VMMCMovetype::move() {
        // generate movemap
        // select monomer1
        // while list not empty
            // apply movemap
            // generate list of interacting pairs
            // randomly select pair (and remove)
            // apply movemap to monomers if applicable
            // calculate energy after applying movemap to first monomer
            // test prelink
            // if fail, continue
            // calculate energy after applying movemap to second monomer
            // test link
            // if fail, record frustrated link, continue
            // check if monomer involved in frustrated link, record if so
            // set new monomer1
        // If frustrated links remain, reject
        // else apply movemap to all linked monomers
    }
}
