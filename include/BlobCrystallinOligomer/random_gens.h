// random_gens.h

#ifndef RANDOM_GENS_H
#define RANDOM_GENS_H

#include <memory>
#include <random>
#include <unordered_map>
#include <utility>

#include "hash.h"

namespace random_gens {

    using std::mt19937_64;
    using std::pair;
    using std::unordered_map;
    using std::uniform_int_distribution;
    using std::uniform_real_distribution;
    using std::unique_ptr;
    using std::vector;

    class RandomGens {
        public:
            /*  Class for storing and using instances of random number generators */

            RandomGens();

            double uniform_real() ;
            /*  Draw a real number uniformly from 0 < x < 1 */

            int uniform_int(int lower, int upper);
            /*  Draw an integer uniformly from lower <= x <= upper. */

        private:
            mt19937_64 m_random_engine {};
            uniform_real_distribution<double> m_uniform_real_dist;
            vector<unique_ptr<uniform_int_distribution<int>>> m_int_dists;
            unordered_map<pair<int, int>, uniform_int_distribution<int>&>
                    m_uniform_int_dists;

    };
}

#endif // RANDOM_GENS_H
