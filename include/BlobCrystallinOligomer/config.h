// config.h

#ifndef CONFIG_H
#define CONFIG_H

#include <memory>
#include <vector>

#include "BlobCrystallinOligomer/ifile.h"
#include "BlobCrystallinOligomer/monomer.h"
#include "BlobCrystallinOligomer/param.h"
#include "BlobCrystallinOligomer/particle.h"
#include "BlobCrystallinOligomer/random_gens.h"
#include "BlobCrystallinOligomer/shared_types.h"
#include "BlobCrystallinOligomer/space.h"

namespace config {

    using ifile::InputConfigFile;
    using ifile::MonomerData;
    using ifile::ParticleData;
    using monomer::Monomer;
    using param::InputParams;
    using particle::Particle;
    using random_gens::RandomGens;
    using shared_types::distT;
    using shared_types::CoorSet;
    using shared_types::vecT;
    using space::CuboidPBC;
    using std::reference_wrapper;
    using std::unique_ptr;
    using std::vector;

    typedef vector<reference_wrapper<Monomer>> monomerArrayT;

    class Config {
        public:
            Config(InputParams params, RandomGens& random_num);

            Monomer& get_monomer(int monomer_index);

            Monomer& get_random_monomer();
            /*  Draw monomer with uniform probability */
            
            monomerArrayT get_monomers();

            int get_num_particles();

            distT get_box_len();

            distT get_radius();

            vecT calc_interparticle_vector(
                    Particle& particle1,
                    CoorSet coorset1,
                    Particle& particle2,
                    CoorSet coorset2);
            /*  Calculate vector from particle 1 to particle 2 */

            distT calc_dist(
                    Particle& particle1,
                    CoorSet coorset1,
                    Particle& particle2,
                    CoorSet coorset2);
            /*  Calculate distance between two particles */

        private:
            vector<unique_ptr<Monomer>> m_monomers;
            monomerArrayT m_monomer_refs;
            unique_ptr<CuboidPBC> m_space_store;
            CuboidPBC& m_space;
            RandomGens& m_random_num;
            distT m_box_len;
            distT m_radius;

            void create_monomers(vector<MonomerData>);
    };
}

#endif // CONFIG_H
