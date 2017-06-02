// monomer.h

#ifndef MONOMER_H
#define MONOMER_H

#include <vector>

#include <Eigen/Dense>

#include "shared_types.h"

namespace monomer {

    using std::vector;
    using shared_types::vecT;

    // Consider making multiple classes for different versions of the alphaB
    // monomer model, or in the distant future, alphaA monomers

    class Monomer {
        // Implementation of alphaB monomer structure
        public:
            void translate(vecT disV);
            // Translate monomer by given vector

            void rotate();
            // Rotate monomer by given ? around given origin

        private:
            int m_index; // Unique monomer index
            vector<int> m_sphere_indices; // Unique sphere indices

            // Patch vectors
            vecT acd1_patch_norm; // Vector normal to ACD1 patch
            vecT acd1_trial_patch_norm; // Trial vector normal to ACD1 patch
            // This is defined as being the normalized vector of the position
            // of the ACD2 - ACD1

            vecT acd1_patch_ref; // Reference vector for torsional alignment
            vecT acd1_trial_patch_ref; // Trieal reference vector for torsional alignment
            // This is defined as being in the plane and direction of the NTD
            // spheres and perpindicular to the patch norm

            vecT ntd1_trial_patch_norm; // Trial vector normal to ACD patch
            // This is defined as being phi degrees plus pi to the ACD patch
            // norm

    };
}

#endif // MONOMER_H
