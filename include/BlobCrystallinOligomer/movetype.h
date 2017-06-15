// movetype.h

#ifndef MOVETYPE_H
#define MOVETYPE_H

#include <set>
#include <utility>
#include <vector>

#include "BlobCrystallinOligomer/config.h"
#include "BlobCrystallinOligomer/energy.h"
#include "BlobCrystallinOligomer/hash.h"
#include "BlobCrystallinOligomer/monomer.h"
#include "BlobCrystallinOligomer/random_gens.h"

namespace movetype {

    using config::Config;
    using config::Monomer;
    using config::monomerArrayT;
    using energy::Energy;
    using random_gens::RandomGens;
    using shared_types::eneT;
    using std::pair;
    using std::set;
    using std::vector;

    class MCMovetype {
        // MC movetype interface
        public:
            MCMovetype(Config& conf, Energy& ene, RandomGens& random_num,
                    eneT m_beta);
            virtual bool move() = 0;

        protected:
            Config& m_config;
            Energy& m_energy;
            RandomGens& m_random_num;
            eneT m_beta;
    };

    class VMMCMovetype: public MCMovetype {
        public:
            using MCMovetype::MCMovetype;
            bool move();

        private:
            monomerArrayT m_cluster;
            int m_frustrated_links {0};
            vector<int> m_frustrated_mis {};
            set<int> m_interacting_mis;
            set<pair<int, int>> m_pair_mis;

            void virtual generate_movemap() = 0;
            void virtual apply_movemap(Monomer& monomer) = 0;

            void add_interacting_pairs(Monomer& monomer1);
            pair<int, int> select_random_pair();
            double calc_prelink_prob(eneT delta_e);
            bool accept_prelink(double prelink_p);
            bool accept_link(double prelink_for_p, double prelink_rev_p);
            bool accept_move();
            void reset_internal();
    };

    class TranslationVMMCMovetype: public VMMCMovetype {
        public:
            using VMMCMovetype::VMMCMovetype;
        private:
            void generate_movemap();
            void apply_movemap(Monomer& monomer);
    };

    class RotationVMMCMovetype: public VMMCMovetype {
        public:
            using VMMCMovetype::VMMCMovetype;
        private:
            void generate_movemap();
            void apply_movemap(Monomer& monomer);
    };
}

#endif // MOVETYPE_H
