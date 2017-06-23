// movetype.h

#ifndef MOVETYPE_H
#define MOVETYPE_H

#include <cmath>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "BlobCrystallinOligomer/config.h"
#include "BlobCrystallinOligomer/energy.h"
#include "BlobCrystallinOligomer/hash.h"
#include "BlobCrystallinOligomer/monomer.h"
#include "BlobCrystallinOligomer/param.h"
#include "BlobCrystallinOligomer/random_gens.h"

namespace movetype {

    using config::Config;
    using config::Monomer;
    using config::monomerArrayT;
    using energy::Energy;
    using param::InputParams;
    using random_gens::RandomGens;
    using shared_types::eneT;
    using shared_types::distT;
    using shared_types::rotMatT;
    using shared_types::vecT;
    using std::pair;
    using std::set;
    using std::sqrt;
    using std::string;
    using std::vector;

    vecT random_unit_vector(RandomGens& random_num);
    distT random_displacement(distT max_disp, RandomGens& random_num);

    class MCMovetype {
        // MC movetype interface
        public:
            MCMovetype(Config& conf, Energy& ene, RandomGens& random_num,
                    InputParams params);
            virtual bool move() = 0;
            virtual string label() {return "MCMovetype";}

        protected:
            Config& m_config;
            Energy& m_energy;
            RandomGens& m_random_num;
            eneT m_beta;
    };

    class NTDFlipMCMovetype: public MCMovetype {
        public:
            using MCMovetype::MCMovetype;
            bool move() {}
            string label() {return "NTDFlipMCMovetype";}
    };

    class VMMCMovetype: public MCMovetype {
        public:
            using MCMovetype::MCMovetype;
            bool move();
            string label() {return "VMMCMovetype";}

        private:
            monomerArrayT m_cluster;
            int m_frustrated_links {0};
            vector<int> m_frustrated_mis {};
            set<int> m_interacting_mis;
            set<pair<int, int>> m_pair_mis;

            void virtual generate_movemap(Monomer& seed_monomer) = 0;
            void virtual apply_movemap(Monomer& monomer) = 0;

            void add_interacting_pairs(Monomer& monomer1);
            pair<int, int> select_random_pair();
            double calc_prelink_prob(eneT ene1, eneT ene2);
            bool accept_prelink(double prelink_p);
            bool accept_link(double prelink_for_p, double prelink_rev_p);
            bool accept_move();
            void reset_internal();
    };

    class TranslationVMMCMovetype: public VMMCMovetype {
        public:
            TranslationVMMCMovetype(Config& conf, Energy& ene, RandomGens& random_num,
                    InputParams params);
            string label() {return "TranslationVMMCMovetype";}
        private:
            void generate_movemap(Monomer&);
            void apply_movemap(Monomer& monomer);

            distT m_max_disp_tc;
            vecT m_disp_v;
    };

    class RotationVMMCMovetype: public VMMCMovetype {
        public:
            RotationVMMCMovetype(Config& conf, Energy& ene, RandomGens& random_num,
                    InputParams params);
            string label() {return "RotationVMMCMovetype";}
        private:
            void generate_movemap(Monomer& seed_monomer);
            void apply_movemap(Monomer& monomer);

            distT m_max_disp_rc;
            distT m_max_disp_a;
            vecT m_rot_c; // Centre of rotation
            rotMatT m_rot_mat; // Rotation matrix
            
    };
}

#endif // MOVETYPE_H
