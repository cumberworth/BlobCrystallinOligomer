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

    /** System configuration container
      *
      * Holds all monomer objects and provides an interface for configuration
      * properties. Responsible for constructing monomers given monomer data.
      */
    class Config {
        public:
            Config(InputParams params, RandomGens& random_num);
            Config(vector<MonomerData> monomers, RandomGens& random_num,
                    distT box_len, distT radius);

            Monomer& get_monomer(int monomer_index);

            /**  Draw monomer with uniform probability */
            Monomer& get_random_monomer();
            
            /** Get all monomers in system */
            monomerArrayT get_monomers();

            int get_num_particles();

            /** This assumes cuboid geometry */
            distT get_box_len();

            /** The radius of a single bead of a monomer
              *
              * Assumes only only bead size.
              */
            distT get_radius();

            /** Calculate vector from particle 1 to particle 2 */
            vecT calc_interparticle_vector(
                    Particle& particle1,
                    CoorSet coorset1,
                    Particle& particle2,
                    CoorSet coorset2);

            /** Calculate distance between two particles */
            distT calc_dist(
                    Particle& particle1,
                    CoorSet coorset1,
                    Particle& particle2,
                    CoorSet coorset2);

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
