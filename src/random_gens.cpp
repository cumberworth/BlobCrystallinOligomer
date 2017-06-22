// random_gens.cpp

#include <memory>
#include <random>

#include "BlobCrystallinOligomer/random_gens.h"

namespace random_gens {

    using std::uniform_int_distribution;
    using std::unique_ptr;

    RandomGens::RandomGens() {
        std::random_device true_random_engine {}; 
        auto seed {true_random_engine()}; 
        m_random_engine.seed(seed); 
    }

    double RandomGens::uniform_real() {
        return m_uniform_real_dist(m_random_engine);
    }

    int RandomGens::uniform_int(int lower, int upper) {

        // Check if distribution used previously
        pair<int, int> key {lower, upper};
        if (m_uniform_int_dists.find(key) != m_uniform_int_dists.end()) {
            auto dist {m_uniform_int_dists.at(key)};
            return dist(m_random_engine);
        }

        // If not, make it and store it
        else {
            uniform_int_distribution<int>* dist {
                    new uniform_int_distribution<int> {lower, upper}};
            int random_int {(*dist)(m_random_engine)};
            m_int_dists.emplace_back(dist);
            m_uniform_int_dists.insert({key, *m_int_dists.back()});
            return random_int;
        }
    }
}
