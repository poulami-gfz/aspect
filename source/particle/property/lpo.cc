/*
  Copyright (C) 2015 - 2017 by the authors of the ASPECT code.

 This file is part of ASPECT.

 ASPECT is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2, or (at your option)
 any later version.

 ASPECT is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with ASPECT; see the file LICENSE.  If not see
 <http://www.gnu.org/licenses/>.
 */

#include <aspect/particle/property/lpo.h>
#include <world_builder/grains.h>
#include <aspect/geometry_model/interface.h>
#include <world_builder/world.h>

#include <aspect/utilities.h>

namespace aspect
{
  namespace Particle
  {
    namespace Property
    {

      template <int dim>
      unsigned int LPO<dim>::n_grains = 0;

      template <int dim>
      unsigned int LPO<dim>::n_minerals = 0;

      template <int dim>
      LPO<dim>::LPO ()
      {
        permutation_operator_3d[0][1][2]  = 1;
        permutation_operator_3d[1][2][0]  = 1;
        permutation_operator_3d[2][0][1]  = 1;
        permutation_operator_3d[0][2][1]  = -1;
        permutation_operator_3d[1][0][2]  = -1;
        permutation_operator_3d[2][1][0]  = -1;
      }

      template <int dim>
      void
      LPO<dim>::initialize ()
      {
        const unsigned int my_rank = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
        this->random_number_generator.seed(random_number_seed+my_rank);
      }

      template <int dim>
      void
      LPO<dim>::load_particle_data(unsigned int lpo_data_position,
                                   const ArrayView<double> &data,
                                   std::vector<unsigned int> &deformation_type,
                                   std::vector<double> &volume_fraction_mineral,
                                   std::vector<std::vector<double>> &volume_fractions_grains,
                                   std::vector<std::vector<Tensor<2,3> > > &a_cosine_matrices_grains)
      {
        // the layout of the data vector per perticle is the following (note that for this plugin the following dim's are always 3):
        // 1. M mineral times
        //    1.1  olivine deformation type   -> 1 double, at location
        //                                      => data_position + 0 + mineral_i * (n_grains * 10 + 2)
        //    2.1. Mineral volume fraction    -> 1 double, at location
        //                                      => data_position + 1 + mineral_i *(n_grains * 10 + 2)
        //    2.2. N grains times:
        //         2.1. volume fraction grain -> 1 double, at location:
        //                                      => data_position + 2 + i_grain * 10 + mineral_i *(n_grains * 10 + 2), or
        //                                      => data_position + 2 + i_grain * (2 * Tensor<2,3>::n_independent_components+ 2) + mineral_i * (n_grains * 10 + 2)
        //         2.2. a_cosine_matrix grain -> 9 (Tensor<2,dim>::n_independent_components) doubles, starts at:
        //                                      => data_position + 3 + i_grain * 10 + mineral_i * (n_grains * 10 + 2), or
        //                                      => data_position + 3 + i_grain * (2 * Tensor<2,3>::n_independent_components+ 2) + mineral_i * (n_grains * 10 + 2)
        //
        // Last used data entry is data_position + 0 + n_minerals * (n_grains * 10 + 2);
        // An other note is that we store exactly the same amount of all minerals (e.g. olivine and enstatite
        // grains), although the volume may not be the same. This has to do with that we need a minimum amount
        // of grains per tracer to perform reliable statistics on it. This miminum should be the same for both
        // olivine and enstatite.
        const unsigned int n_grains = LPO<dim>::n_grains;
        const unsigned int n_minerals = LPO<dim>::n_minerals;
        deformation_type.resize(n_minerals);
        volume_fraction_mineral.resize(n_minerals);
        volume_fractions_grains.resize(n_minerals);
        a_cosine_matrices_grains.resize(n_minerals);

        for (size_t mineral_i = 0; mineral_i < n_minerals; mineral_i++)
          {
            deformation_type[mineral_i] = data[lpo_data_position + 0 + mineral_i * (n_grains * 10 + 2)];
            volume_fraction_mineral[mineral_i] = data[lpo_data_position + 1 + mineral_i * (n_grains * 10 + 2)];
            volume_fractions_grains[mineral_i].resize(n_grains);
            a_cosine_matrices_grains[mineral_i].resize(n_grains);
            // loop over grains to store the data of each grain
            for (unsigned int grain_i = 0; grain_i < n_grains; ++grain_i)
              {
                // store volume fraction for olivine grains
                volume_fractions_grains[mineral_i][grain_i] = data[lpo_data_position + 2 + grain_i * 10 + mineral_i *(n_grains * 10 + 2)];

                // store a_{ij} for olivine grains
                for (unsigned int i = 0; i < Tensor<2,3>::n_independent_components ; ++i)
                  {
                    const dealii::TableIndices<2> index = Tensor<2,3>::unrolled_to_component_indices(i);
                    a_cosine_matrices_grains[mineral_i][grain_i][index] = data[lpo_data_position + 3 + grain_i * 10 + mineral_i * (n_grains * 10 + 2) + i];
                  }
              }
          }
      }


      template <int dim>
      void
      LPO<dim>::load_particle_data_extended(unsigned int lpo_data_position,
                                            const ArrayView<double> &data,
                                            std::vector<unsigned int> &deformation_type,
                                            std::vector<double> &volume_fraction_mineral,
                                            std::vector<std::vector<double>> &volume_fractions_grains,
                                            std::vector<std::vector<Tensor<2,3> > > &a_cosine_matrices_grains,
                                            std::vector<std::vector<double> > &volume_fractions_grains_derivatives,
                                            std::vector<std::vector<Tensor<2,3> > > &a_cosine_matrices_grains_derivatives) const
      {
        load_particle_data(lpo_data_position,
                           data,
                           deformation_type,
                           volume_fraction_mineral,
                           volume_fractions_grains,
                           a_cosine_matrices_grains);

        // now store the derivatives if needed
        if (this->advection_method == AdvectionMethod::CrankNicolson)
          {
            volume_fractions_grains_derivatives.resize(n_minerals);
            a_cosine_matrices_grains_derivatives.resize(n_minerals);

            for (size_t mineral_i = 0; mineral_i < n_minerals; mineral_i++)
              {
                volume_fractions_grains_derivatives[mineral_i].resize(n_grains);
                a_cosine_matrices_grains_derivatives[mineral_i].resize(n_grains);

                for (unsigned int grain_i = 0; grain_i < n_grains; ++grain_i)
                  {
                    // store volume fraction for olivine grains
                    volume_fractions_grains_derivatives[mineral_i][grain_i] = data[lpo_data_position  + n_minerals * (n_grains * 10 + 2) + mineral_i * (n_grains * 10)  + grain_i * 10 ];

                    // store a_{ij} for olivine grains
                    for (unsigned int iii = 0; iii < Tensor<2,3>::n_independent_components ; ++iii)
                      {
                        const dealii::TableIndices<2> index = Tensor<2,3>::unrolled_to_component_indices(iii);
                        a_cosine_matrices_grains_derivatives[mineral_i][grain_i][index] = data[lpo_data_position + n_minerals * (n_grains * 10 + 2) + mineral_i * (n_grains * 10)  + grain_i * 10  + 1 + iii];
                      }
                  }
              }
          }
      }


      template <int dim>
      void
      LPO<dim>::store_particle_data(unsigned int lpo_data_position,
                                    const ArrayView<double> &data,
                                    std::vector<unsigned int> &deformation_type,
                                    std::vector<double> &volume_fraction_mineral,
                                    std::vector<std::vector<double>> &volume_fractions_grains,
                                    std::vector<std::vector<Tensor<2,3> > > &a_cosine_matrices_grains)
      {
        // the layout of the data vector per perticle is the following (note that for this plugin the following dim's are always 3):
        // 1. M mineral times
        //    1.1  olivine deformation type   -> 1 double, at location
        //                                      => data_position + 0 + mineral_i * (n_grains * 10 + 2)
        //    2.1. Mineral volume fraction    -> 1 double, at location
        //                                      => data_position + 1 + mineral_i *(n_grains * 10 + 2)
        //    2.2. N grains times:
        //         2.1. volume fraction grain -> 1 double, at location:
        //                                      => data_position + 2 + i_grain * 10 + mineral_i *(n_grains * 10 + 2), or
        //                                      => data_position + 2 + i_grain * (2 * Tensor<2,3>::n_independent_components+ 2) + mineral_i * (n_grains * 10 + 2)
        //         2.2. a_cosine_matrix grain -> 9 (Tensor<2,dim>::n_independent_components) doubles, starts at:
        //                                      => data_position + 3 + i_grain * 10 + mineral_i * (n_grains * 10 + 2), or
        //                                      => data_position + 3 + i_grain * (2 * Tensor<2,3>::n_independent_components+ 2) + mineral_i * (n_grains * 10 + 2)
        //
        // Last used data entry is data_position + 0 + n_minerals * (n_grains * 10 + 2);
        // An other note is that we store exactly the same amount of all minerals (e.g. olivine and enstatite
        // grains), although the volume may not be the same. This has to do with that we need a minimum amount
        // of grains per tracer to perform reliable statistics on it. This miminum should be the same for both
        // olivine and enstatite.
        const unsigned int n_grains = LPO<dim>::n_grains;
        const unsigned int n_minerals = LPO<dim>::n_minerals;
        for (size_t mineral_i = 0; mineral_i < n_minerals; mineral_i++)
          {
            Assert(volume_fractions_grains[mineral_i].size() == n_grains, ExcMessage("Internal error: volume_fractions_mineral[mineral_i] is not the same as n_grains."));
            Assert(a_cosine_matrices_grains[mineral_i].size() == n_grains, ExcMessage("Internal error: a_cosine_matrices_mineral[mineral_i] is not the same as n_grains."));
            data[lpo_data_position + 0 + mineral_i * (n_grains * 10 + 2)] = deformation_type[mineral_i];
            data[lpo_data_position + 1 + mineral_i *(n_grains * 10 + 2)] = volume_fraction_mineral[mineral_i];

            // loop over grains to store the data of each grain
            for (unsigned int grain_i = 0; grain_i < n_grains; ++grain_i)
              {
                // store volume fraction for olivine grains
                data[lpo_data_position + 2 + grain_i * 10 + mineral_i *(n_grains * 10 + 2)] = volume_fractions_grains[mineral_i][grain_i];
                // store a_{ij} for olivine grains
                for (unsigned int i = 0; i < Tensor<2,3>::n_independent_components ; ++i)
                  {
                    const dealii::TableIndices<2> index = Tensor<2,3>::unrolled_to_component_indices(i);
                    data[lpo_data_position + 3 + grain_i * 10 + mineral_i * (n_grains * 10 + 2) + i] = a_cosine_matrices_grains[mineral_i][grain_i][index];
                  }
              }
          }
      }



      template <int dim>
      void
      LPO<dim>::store_particle_data_extended(unsigned int lpo_data_position,
                                             const ArrayView<double> &data,
                                             std::vector<unsigned int> &deformation_type,
                                             std::vector<double> &volume_fraction_mineral,
                                             std::vector<std::vector<double>> &volume_fractions_grains,
                                             std::vector<std::vector<Tensor<2,3> > > &a_cosine_matrices_grains,
                                             std::vector<std::vector<double> > &volume_fractions_grains_derivatives,
                                             std::vector<std::vector<Tensor<2,3> > > &a_cosine_matrices_grains_derivatives) const
      {
        store_particle_data(lpo_data_position,
                            data,
                            deformation_type,
                            volume_fraction_mineral,
                            volume_fractions_grains,
                            a_cosine_matrices_grains);
        // now store the derivatives if needed. They are added after all the other data.
        // lpo_data_position + 3 + n_grains * 10 + mineral_i * (n_grains * 10 + 2)
        // data_position + 0 + n_minerals * (n_grains * 10 + 2)
        if (this->advection_method == AdvectionMethod::CrankNicolson)
          {
            for (size_t mineral_i = 0; mineral_i < n_minerals; mineral_i++)
              {
                Assert(volume_fractions_grains_derivatives.size() == n_minerals, ExcMessage("Internal error: volume_fractions_olivine_derivatives is not the same as n_minerals."));
                Assert(a_cosine_matrices_grains_derivatives.size() == n_minerals, ExcMessage("Internal error: a_cosine_matrices_olivine_derivatives is not the same as n_minerals."));

                for (size_t mineral_i = 0; mineral_i < n_minerals; mineral_i++)
                  {
                    Assert(volume_fractions_grains_derivatives[mineral_i].size() == n_grains, ExcMessage("Internal error: volume_fractions_olivine_derivatives is not the same as n_grains."));
                    Assert(a_cosine_matrices_grains_derivatives[mineral_i].size() == n_grains, ExcMessage("Internal error: a_cosine_matrices_olivine_derivatives is not the same as n_grains."));

                    for (unsigned int grain_i = 0; grain_i < n_grains; ++grain_i)
                      {
                        // store volume fraction for olivine grains
                        data[lpo_data_position + n_minerals * (n_grains * 10 + 2) + mineral_i * (n_grains * 10)  + grain_i * 10] = volume_fractions_grains_derivatives[mineral_i][grain_i];

                        // store a_{ij} for olivine grains
                        for (unsigned int iii = 0; iii < Tensor<2,3>::n_independent_components ; ++iii)
                          {
                            const dealii::TableIndices<2> index = Tensor<2,3>::unrolled_to_component_indices(iii);
                            data[lpo_data_position + n_minerals * (n_grains * 10 + 2) + mineral_i * (n_grains * 10)  + grain_i * 10  + 1 + iii] = a_cosine_matrices_grains_derivatives[mineral_i][grain_i][index];
                          }
                      }
                  }
              }
          }
      }

      template <int dim>
      void
      LPO<dim>::initialize_one_particle_property(const Point<dim> &position,
                                                 std::vector<double> &data) const
      {
        // the layout of the data vector per perticle is the following (note that for this plugin the following dim's are always 3):
        // 1. M mineral times
        //    1.1  olivine deformation type   -> 1 double, at location
        //                                      => data_position + 0 + mineral_i * (n_grains * 10 + 2)
        //    2.1. Mineral volume fraction    -> 1 double, at location
        //                                      => data_position + 1 + mineral_i *(n_grains * 10 + 2)
        //    2.2. N grains times:
        //         2.1. volume fraction grain -> 1 double, at location:
        //                                      => data_position + 2 + i_grain * 10 + mineral_i *(n_grains * 10 + 2), or
        //                                      => data_position + 2 + i_grain * (2 * Tensor<2,3>::n_independent_components+ 2) + mineral_i * (n_grains * 10 + 2)
        //         2.2. a_cosine_matrix grain -> 9 (Tensor<2,dim>::n_independent_components) doubles, starts at:
        //                                      => data_position + 3 + i_grain * 10 + mineral_i * (n_grains * 10 + 2), or
        //                                      => data_position + 3 + i_grain * (2 * Tensor<2,3>::n_independent_components+ 2) + mineral_i * (n_grains * 10 + 2)
        //
        // Note is that we store exactly the same amount of all minerals (e.g. olivine and enstatite
        // grains), although the volume may not be the same. This has to do with that we need a minimum amount
        // of grains per tracer to perform reliable statistics on it. This miminum is the same for both olivine
        // and enstatite.

        // fabric. This is determined in the computations, so set it to -1 for now.
        std::vector<double> deformation_type(n_minerals, -1.0);
        std::vector<std::vector<double > >volume_fractions_grains(n_minerals);
        std::vector<std::vector<Tensor<2,3> > > a_cosine_matrices_grains(n_minerals);

        for (size_t mineral_i = 0; mineral_i < n_minerals; mineral_i++)
          {

            volume_fractions_grains[mineral_i].resize(n_grains);
            a_cosine_matrices_grains[mineral_i].resize(n_grains);

            if (use_world_builder == true)
              {
#ifdef ASPECT_WITH_WORLD_BUILDER
                WorldBuilder::grains wb_grains = this->get_world_builder().grains(Utilities::convert_point_to_array(position),
                                                                                  -this->get_geometry_model().height_above_reference_surface(position),
                                                                                  mineral_i,
                                                                                  n_grains);
                double sum_volume_fractions = 0;
                for (unsigned int grain_i = 0; grain_i < n_grains ; ++grain_i)
                  {
                    sum_volume_fractions += wb_grains.sizes[grain_i];
                    volume_fractions_grains[mineral_i][grain_i] = wb_grains.sizes[grain_i];
                    // we are receiving a array<array<double,3>,3> which nees to be unrolled in the correct way
                    // for a tensor<2,3> it just loops first over the second index and than the first index
                    for (unsigned int component_i = 0; component_i < 3 ; ++component_i)
                      {
                        for (unsigned int component_j = 0; component_j < 3 ; ++component_j)
                          {
                            Assert(!std::isnan(wb_grains.rotation_matrices[grain_i][component_i][component_j]), ExcMessage("Error: not a number."));
                            a_cosine_matrices_grains[mineral_i][grain_i][component_i][component_j] = wb_grains.rotation_matrices[grain_i][component_i][component_j];
                          }
                      }

                  }

                AssertThrow(sum_volume_fractions != 0, ExcMessage("Sum of volumes is equal to zero, which is not supporsed to happen. "
                                                                  "Make sure that all parts of the domain which contain particles are covered by the world builder."));
#else
                AssertThrow(false,
                            "The world builder was requested but not provided. Make sure that aspect is "
                            "compiled with the World Builder and that you provide a world builder file in the input.")
#endif
              }
            else
              {
                const double initial_volume_fraction = 1.0/n_grains;
                boost::random::uniform_real_distribution<double> uniform_distribution(0,1);

                // it is 2 times the amount of grains because we have to compute for olivine and enstatite.
                for (unsigned int grain_i = 0; grain_i < n_grains ; ++grain_i)
                  {
                    // set volume fraction
                    volume_fractions_grains[mineral_i][grain_i] = initial_volume_fraction;

                    // set a uniform random a_cosine_matrix per grain
                    // This function is based on an article in Graphic Gems III, written by James Arvo, Cornell University (p 116-120).
                    // The original code can be found on  http://www.realtimerendering.com/resources/GraphicsGems/gemsiii/rand_rotation.c
                    // and is licenend accourding to this website with the following licence:
                    //
                    // "The Graphics Gems code is copyright-protected. In other words, you cannot claim the text of the code as your own and
                    // resell it. Using the code is permitted in any program, product, or library, non-commercial or commercial. Giving credit
                    // is not required, though is a nice gesture. The code comes as-is, and if there are any flaws or problems with any Gems
                    // code, nobody involved with Gems - authors, editors, publishers, or webmasters - are to be held responsible. Basically,
                    // don't be a jerk, and remember that anything free comes with no guarantee.""
                    //
                    // The book saids in the preface the following: "As in the first two volumes, all of the C and C++ code in this book is in
                    // the public domain, and is yours to study, modify, and use."

                    // first generate three random numbers between 0 and 1 and multiply them with 2 PI or 2 for z. Note that these are not the same as phi_1, theta and phi_2.
                    double one = uniform_distribution(this->random_number_generator);
                    double two = uniform_distribution(this->random_number_generator);
                    double three = uniform_distribution(this->random_number_generator);

                    double theta = 2.0 * M_PI * one; // Rotation about the pole (Z)
                    double phi = 2.0 * M_PI * two; // For direction of pole deflection.
                    double z = 2.0 * three; //For magnitude of pole deflection.

                    // Compute a vector V used for distributing points over the sphere
                    // via the reflection I - V Transpose(V).  This formulation of V
                    // will guarantee that if x[1] and x[2] are uniformly distributed,
                    // the reflected points will be uniform on the sphere.  Note that V
                    // has length sqrt(2) to eliminate the 2 in the Householder matrix.

                    double r  = std::sqrt( z );
                    double Vx = std::sin( phi ) * r;
                    double Vy = std::cos( phi ) * r;
                    double Vz = std::sqrt( 2.0 - z );

                    // Compute the row vector S = Transpose(V) * R, where R is a simple
                    // rotation by theta about the z-axis.  No need to compute Sz since
                    // it's just Vz.

                    double st = std::sin( theta );
                    double ct = std::cos( theta );
                    double Sx = Vx * ct - Vy * st;
                    double Sy = Vx * st + Vy * ct;

                    // Construct the rotation matrix  ( V Transpose(V) - I ) R, which
                    // is equivalent to V S - R.

                    a_cosine_matrices_grains[mineral_i][grain_i][0][0] = Vx * Sx - ct;
                    a_cosine_matrices_grains[mineral_i][grain_i][0][1] = Vx * Sy - st;
                    a_cosine_matrices_grains[mineral_i][grain_i][0][2] = Vx * Vz;

                    a_cosine_matrices_grains[mineral_i][grain_i][1][0] = Vy * Sx + st;
                    a_cosine_matrices_grains[mineral_i][grain_i][1][1] = Vy * Sy - ct;
                    a_cosine_matrices_grains[mineral_i][grain_i][1][2] = Vy * Vz;

                    a_cosine_matrices_grains[mineral_i][grain_i][2][0] = Vz * Sx;
                    a_cosine_matrices_grains[mineral_i][grain_i][2][1] = Vz * Sy;
                    a_cosine_matrices_grains[mineral_i][grain_i][2][2] = 1.0 - z;   // This equals Vz * Vz - 1.0

                  }
              }
          }

        for (size_t mineral_i = 0; mineral_i < n_minerals; mineral_i++)
          {
            data.emplace_back(deformation_type[mineral_i]);
            data.emplace_back(volume_fractions_minerals[mineral_i]);
            for (unsigned int grain_i = 0; grain_i < n_grains ; ++grain_i)
              {
                data.emplace_back(volume_fractions_grains[mineral_i][grain_i]);
                for (unsigned int i = 0; i < Tensor<2,3>::n_independent_components ; ++i)
                  {
                    const dealii::TableIndices<2> index = Tensor<2,3>::unrolled_to_component_indices(i);
                    data.emplace_back(a_cosine_matrices_grains[mineral_i][grain_i][index]);
                  }

              }
          }

        if (this->advection_method == AdvectionMethod::CrankNicolson)
          {
            // start with derivatives set to zero
            for (size_t mineral_i = 0; mineral_i < n_minerals; mineral_i++)
              {
                for (unsigned int i_grain = 0; i_grain < n_grains ; ++i_grain)
                  {
                    data.push_back(0.0);
                    for (unsigned int index = 0; index < Tensor<2,3>::n_independent_components; index++)
                      {
                        data.push_back(0.0);
                      }
                  }
              }
          }
      }

      template<int dim>
      std::vector<Tensor<2,3> >
      LPO<dim>::random_draw_volume_weighting(std::vector<double> fv,
                                             std::vector<Tensor<2,3>> matrices,
                                             unsigned int n_output_grains) const
      {
        // Get volume weighted euler angles, using random draws to convert odf
        // to a discrete number of orientations, weighted by volume
        // 1a. Get index that would sort volume fractions AND
        // ix = np.argsort(fv[q,:]);
        // 1b. Get the sorted volume and angle arrays
        std::vector<double> fv_to_sort = fv;
        std::vector<double> fv_sorted = fv;
        std::vector<Tensor<2,3>> matrices_sorted = matrices;

        unsigned int n_grain = fv_to_sort.size();


        /**
         * ...
         */
        for (int i = n_grain-1; i >= 0; --i)
          {
            unsigned int index_max_fv = std::distance(fv_to_sort.begin(),max_element(fv_to_sort.begin(), fv_to_sort.end()));

            fv_sorted[i] = fv_to_sort[index_max_fv];
            matrices_sorted[i] = matrices[index_max_fv];
            /*Assert(matrices[index_max_fv].size() == 3, ExcMessage("matrices vector (size = " + std::to_string(matrices[index_max_fv].size()) +
                                                                ") should have size 3."));
            Assert(matrices_sorted[i].size() == 3, ExcMessage("matrices_sorted vector (size = " + std::to_string(matrices_sorted[i].size()) +
                                                            ") should have size 3."));*/
            fv_to_sort[index_max_fv] = -1;
          }

        // 2. Get cumulative weight for volume fraction
        std::vector<double> cum_weight(n_grains);
        std::partial_sum(fv_sorted.begin(),fv_sorted.end(),cum_weight.begin());
        // 3. Generate random indices
        boost::random::uniform_real_distribution<> dist(0, 1);
        std::vector<double> idxgrain(n_output_grains);
        for (unsigned int grain_i = 0; grain_i < n_output_grains; ++grain_i)
          {
            idxgrain[grain_i] = dist(this->random_number_generator);
          }

        // 4. Find the maximum cum_weight that is less than the random value.
        // the euler angle index is +1. For example, if the idxGrain(g) < cumWeight(1),
        // the index should be 1 not zero)
        std::vector<Tensor<2,3>> matrices_out(n_output_grains);
        for (unsigned int grain_i = 0; grain_i < n_output_grains; ++grain_i)
          {
            unsigned int counter = 0;
            for (unsigned int grain_j = 0; grain_j < n_grains; ++grain_j)
              {
                // find the first cummulative weight which is larger than the random number
                // todo: there are algorithms to do this faster
                if (cum_weight[grain_j] < idxgrain[grain_i])
                  {
                    counter++;
                  }
                else
                  {
                    break;
                  }
              }
            matrices_out[grain_i] = matrices_sorted[counter];
          }
        return matrices_out;
      }

      template <int dim>
      void
      LPO<dim>::update_one_particle_property(const unsigned int data_position,
                                             const Point<dim> &position,
                                             const Vector<double> &solution,
                                             const std::vector<Tensor<1,dim> > &gradients,
                                             const ArrayView<double> &data) const
      {
        // STEP 1: Load data and preprocess it.

        // need access to the pressure, viscosity,
        // get velocity
        Tensor<1,dim> velocity;
        for (unsigned int i = 0; i < dim; ++i)
          velocity[i] = solution[this->introspection().component_indices.velocities[i]];


        // get velocity gradient tensor.
        Tensor<2,dim> velocity_gradient;
        for (unsigned int d=0; d<dim; ++d)
          velocity_gradient[d] = gradients[d];

        // Calculate strain rate from velocity gradients
        const SymmetricTensor<2,dim> strain_rate = symmetrize (velocity_gradient);

        // compute local stress tensor
        const SymmetricTensor<2,dim> compressible_strain_rate
          = (this->get_material_model().is_compressible()
             ?
             strain_rate - 1./3 * trace(strain_rate) * unit_symmetric_tensor<dim>()
             :
             strain_rate);


        // To determine the deformation type of grains based on figure 4 of Karato et al.,
        // 2008 (Geodynamic Significance of seismic anisotropy o fthe Upper Mantle: New
        // insights from laboratory studies), we need to know the stress and water content.
        // The water content is stored on every particle and the stress is computed here.
        // To compute the stress we need the pressure, the compressible_strain_rate and the
        // viscosity at the location of the particle.

        double pressure = solution[this->introspection().component_indices.pressure];
        double temperature = solution[this->introspection().component_indices.temperature];
        // Only assert in debug mode, because it should already be checked during initialization.
        AssertThrow(this->introspection().compositional_name_exists("water"),
                    ExcMessage("Particle property LPO only works if"
                               "there is a compositional field called water."));
        const unsigned int water_idx = this->introspection().compositional_index_for_name("water");
        double water_content = solution[this->introspection().component_indices.compositional_fields[water_idx]];

        const double dt = this->get_timestep();


        // get the composition of the particle
        std::vector<double> compositions;
        for (unsigned int i = 0; i < this->n_compositional_fields(); i++)
          {
            const unsigned int solution_component = this->introspection().component_indices.compositional_fields[i];
            compositions.push_back(solution[solution_component]);
          }

        // even in 2d we need 3d strain-rates and velocity gradient tensors. So we make them 3d by
        // adding an extra dimension which is zero.
        // Todo: for now we just add a thrid row and column
        // and make them zero. We have to check whether that is correct.
        SymmetricTensor<2,3> strain_rate_3d;
        strain_rate_3d[0][0] = strain_rate[0][0];
        strain_rate_3d[0][1] = strain_rate[0][1];
        //sym: strain_rate_nondim_3d[0][0] = strain_rate[1][0];
        strain_rate_3d[1][1] = strain_rate[1][1];

        if (dim == 3)
          {
            strain_rate_3d[0][2] = strain_rate[0][2];
            strain_rate_3d[1][2] = strain_rate[1][2];
            //sym: strain_rate_nondim_3d[0][0] = strain_rate[2][0];
            //sym: strain_rate_nondim_3d[0][1] = strain_rate[2][1];
            strain_rate_3d[2][2] = strain_rate[2][2];
          }
        Tensor<2,3> velocity_gradient_3d;
        velocity_gradient_3d[0][0] = velocity_gradient[0][0];
        velocity_gradient_3d[0][1] = velocity_gradient[0][1];
        velocity_gradient_3d[1][0] = velocity_gradient[1][0];
        velocity_gradient_3d[1][1] = velocity_gradient[1][1];
        if (dim == 3)
          {
            velocity_gradient_3d[0][2] = velocity_gradient[0][2];
            velocity_gradient_3d[1][2] = velocity_gradient[1][2];
            velocity_gradient_3d[2][0] = velocity_gradient[2][0];
            velocity_gradient_3d[2][1] = velocity_gradient[2][1];
            velocity_gradient_3d[2][2] = velocity_gradient[2][2];
          }

        std::vector<unsigned int> deformation_types;
        std::vector<double> volume_fraction_mineral;
        std::vector<std::vector<double>> volume_fractions_grains;
        std::vector<std::vector<Tensor<2,3> > > a_cosine_matrices_grains;
        std::vector<std::vector<double> > volume_fractions_grains_derivatives;
        std::vector<std::vector<Tensor<2,3> > > a_cosine_matrices_grains_derivatives;

        load_particle_data_extended(data_position,
                                    data,
                                    deformation_types,
                                    volume_fraction_mineral,
                                    volume_fractions_grains,
                                    a_cosine_matrices_grains,
                                    volume_fractions_grains_derivatives,
                                    a_cosine_matrices_grains_derivatives);


        for (size_t mineral_i = 0; mineral_i < n_minerals; mineral_i++)
          {

            // Now compute what type of deformation takes place.
            DeformationType deformation_type = DeformationType::OlivineAFabric;

            switch (deformation_type_selector[mineral_i])
              {
                case DeformationTypeSelector::Passive:
                  deformation_type =  DeformationType::Passive;
                  break;
                case DeformationTypeSelector::OlivineAFabric:
                  deformation_type =  DeformationType::OlivineAFabric;
                  break;
                case DeformationTypeSelector::OlivineBFabric:
                  deformation_type =  DeformationType::OlivineBFabric;
                  break;
                case DeformationTypeSelector::OlivineCFabric:
                  deformation_type =  DeformationType::OlivineCFabric;
                  break;
                case DeformationTypeSelector::OlivineDFabric:
                  deformation_type =  DeformationType::OlivineDFabric;
                  break;
               // case DeformationTypeSelector::OlivineAFabric:
               //   deformation_type =  DeformationType::OlivineAFabric;
               //   break;
                case DeformationTypeSelector::Enstatite:
                  deformation_type =  DeformationType::Enstatite;
                  break;
              //  case DeformationTypeSelector::Enstatite:
               //   deformation_type =  DeformationType::EnstatiteBFabric;
                //  break;
               // case DeformationTypeSelector::Enstatite:
               //   deformation_type =  DeformationType::EnstatiteCFabric;
               //   break;
              //  case DeformationTypeSelector::Enstatite:
              //    deformation_type =  DeformationType::EnstatiteDFabric;
               //   break;
                case DeformationTypeSelector::OlivineKarato2008:
                  // construct the material model inputs and outputs
                  // Since this function is only evaluating one particle,
                  // we use 1 for the amount of quadrature points.
                  MaterialModel::MaterialModelInputs<dim> material_model_inputs(1,this->n_compositional_fields());
                  material_model_inputs.position[0] = position;
                  material_model_inputs.temperature[0] = temperature;
                  material_model_inputs.pressure[0] = pressure;
                  material_model_inputs.velocity[0] = velocity;
                  material_model_inputs.composition[0] = compositions;
                  material_model_inputs.strain_rate[0] = strain_rate;

                  MaterialModel::MaterialModelOutputs<dim> material_model_outputs(1,this->n_compositional_fields());
                  this->get_material_model().evaluate(material_model_inputs, material_model_outputs);
                  double eta = material_model_outputs.viscosities[0];

                  const SymmetricTensor<2,dim> stress = 2*eta*compressible_strain_rate +
                                                        pressure * unit_symmetric_tensor<dim>();
                  const std::array< double, dim > eigenvalues = dealii::eigenvalues(stress);
                  double differential_stress = eigenvalues[0]-eigenvalues[dim-1];
                  deformation_type = determine_deformation_type(differential_stress, water_content);

                  break;
              }

            deformation_types[mineral_i] = (unsigned int)deformation_type;

            const std::array<double,4> ref_resolved_shear_stress = reference_resolved_shear_stress_from_deformation_type(deformation_type);

            for (unsigned int grain_i = 0; grain_i < n_grains; ++grain_i)
              {
                Assert(isfinite(volume_fractions_grains[mineral_i][grain_i]),
                       ExcMessage("volume_fractions_grains[" + std::to_string(grain_i) + "] is not finite directly after loading: "
                                  + std::to_string(volume_fractions_grains[mineral_i][grain_i]) + "."));
              }

            for (unsigned int i = 0; i < n_grains; ++i)
              {
                for (size_t j = 0; j < 3; j++)
                  {
                    for (size_t k = 0; k < 3; k++)
                      {
                        Assert(!std::isnan(a_cosine_matrices_grains[mineral_i][i][j][k]), ExcMessage(" a_cosine_matrices_grains[mineral_i] is nan directly after loading."));
                      }
                  }
              }


            for (unsigned int grain_i = 0; grain_i < n_grains; ++grain_i)
              {
                for (size_t i = 0; i < 3; i++)
                  for (size_t j = 0; j < 3; j++)
                    Assert(fabs(a_cosine_matrices_grains[mineral_i][grain_i][i][j]) -1.0 <= 2.0 * std::numeric_limits<double>::epsilon(),
                           ExcMessage("1. a_cosine_matrices_grains[[" + std::to_string(i) + "][" + std::to_string(j) +
                                      "] is larger than one: " + std::to_string(a_cosine_matrices_grains[mineral_i][grain_i][i][j]) + ". rotation_matrix = \n"
                                      + std::to_string(a_cosine_matrices_grains[mineral_i][grain_i][0][0]) + " " + std::to_string(a_cosine_matrices_grains[mineral_i][grain_i][0][1]) + " " + std::to_string(a_cosine_matrices_grains[mineral_i][grain_i][0][2]) + "\n"
                                      + std::to_string(a_cosine_matrices_grains[mineral_i][grain_i][1][0]) + " " + std::to_string(a_cosine_matrices_grains[mineral_i][grain_i][1][1]) + " " + std::to_string(a_cosine_matrices_grains[mineral_i][grain_i][1][2]) + "\n"
                                      + std::to_string(a_cosine_matrices_grains[mineral_i][grain_i][2][0]) + " " + std::to_string(a_cosine_matrices_grains[mineral_i][grain_i][2][1]) + " " + std::to_string(a_cosine_matrices_grains[mineral_i][grain_i][2][2])));
              }

            /**
            * Now we have loaded all the data and can do the actual computation.
            * The computation consitsts of two parts. The first part is computing
            * the derivatives for the directions and grain sizes. Then those
            * derivatives are used to advect the particle properties.
            */
            double sum_volume_mineral = 0;

            std::pair<std::vector<double>, std::vector<Tensor<2,3> > > derivatives_grains = this->compute_derivatives(volume_fractions_grains[mineral_i],
                                                                                            a_cosine_matrices_grains[mineral_i],
                                                                                            strain_rate_3d,
                                                                                            velocity_gradient_3d,
                                                                                            volume_fraction_mineral[mineral_i],
                                                                                            ref_resolved_shear_stress);

            for (unsigned int grain_i = 0; grain_i < n_grains; ++grain_i)
              {
                Assert(isfinite(volume_fractions_grains[mineral_i][grain_i]),
                       ExcMessage("volume_fractions_grains[mineral_i][" + std::to_string(grain_i) + "] is not finite: "
                                  + std::to_string(volume_fractions_grains[mineral_i][grain_i]) + "."));
              }

            switch (advection_method)
              {
                case AdvectionMethod::ForwardEuler:

                  sum_volume_mineral = this->advect_forward_euler(volume_fractions_grains[mineral_i],
                                                                  a_cosine_matrices_grains[mineral_i],
                                                                  derivatives_grains,
                                                                  dt);

                  break;

                case AdvectionMethod::BackwardEuler:
                  sum_volume_mineral = this->advect_backward_euler(volume_fractions_grains[mineral_i],
                                                                   a_cosine_matrices_grains[mineral_i],
                                                                   derivatives_grains,
                                                                   dt);

                  break;

                case AdvectionMethod::CrankNicolson:
                  sum_volume_mineral = this->advect_Crank_Nicolson(volume_fractions_grains[mineral_i],
                                                                   a_cosine_matrices_grains[mineral_i],
                                                                   derivatives_grains,
                                                                   volume_fractions_grains_derivatives[mineral_i],
                                                                   a_cosine_matrices_grains_derivatives[mineral_i],
                                                                   dt);

                  break;
              }

            // normalize both the olivine and enstatite volume fractions
            const double inv_sum_volume_mineral = 1.0/sum_volume_mineral;

            Assert(std::isfinite(inv_sum_volume_mineral),
                   ExcMessage("inv_sum_volume_mineral is not finite. sum_volume_enstatite = "
                              + std::to_string(sum_volume_mineral)));

            for (unsigned int grain_i = 0; grain_i < n_grains; ++grain_i)
              {
                volume_fractions_grains[mineral_i][grain_i] *= inv_sum_volume_mineral;
                Assert(isfinite(volume_fractions_grains[mineral_i][grain_i]),
                       ExcMessage("volume_fractions_grains[mineral_i]" + std::to_string(grain_i) + "] is not finite: "
                                  + std::to_string(volume_fractions_grains[mineral_i][grain_i]) + ", inv_sum_volume_mineral = "
                                  + std::to_string(inv_sum_volume_mineral) + "."));
              }

            for (unsigned int i = 0; i < n_grains; ++i)
              {
                for (size_t j = 0; j < 3; j++)
                  {
                    for (size_t k = 0; k < 3; k++)
                      {
                        Assert(!std::isnan(a_cosine_matrices_grains[mineral_i][i][j][k]), ExcMessage(" a_cosine_matrices_grains is nan before orthoganalization."));
                      }

                  }

              }



            /**
             * Correct direction cosine matrices numerical error (orthnormality) after integration
             * Follows same method as in matlab version from Thissen of finding the nearest orthonormal
             * matrix using the SVD
             */
            for (unsigned int grain_i = 0; grain_i < n_grains; ++grain_i)
              {
                a_cosine_matrices_grains[mineral_i][grain_i] = dealii::project_onto_orthogonal_tensors(a_cosine_matrices_grains[mineral_i][grain_i]);
                for (size_t i = 0; i < 3; i++)
                  for (size_t j = 0; j < 3; j++)
                    {
                      // I don't think this should happen with the projection, but D-Rex
                      // does not do the orthogonal projection, but just clamps the values
                      // to 1 and -1.
                      Assert(std::fabs(a_cosine_matrices_grains[mineral_i][grain_i][i][j])-1.0 <= 2.0 * std::numeric_limits<double>::epsilon(),
                             ExcMessage("The a_cosine_matrices_grains[mineral_i] has a entry asolute larger than 1:" + std::to_string(std::fabs(a_cosine_matrices_grains[mineral_i][grain_i][i][j])) +"."));
                    }
              }

            for (unsigned int i = 0; i < n_grains; ++i)
              {
                for (size_t j = 0; j < 3; j++)
                  {
                    for (size_t k = 0; k < 3; k++)
                      {
                        Assert(!std::isnan(a_cosine_matrices_grains[mineral_i][i][j][k]),
                               ExcMessage(" a_cosine_matrices_grains[mineral_i] is nan after orthoganalization: "
                                          + std::to_string(a_cosine_matrices_grains[mineral_i][i][j][k])));
                      }

                  }

              }


            for (unsigned int grain_i = 0; grain_i < n_grains; ++grain_i)
              {
                for (size_t i = 0; i < 3; i++)
                  for (size_t j = 0; j < 3; j++)
                    Assert(fabs(a_cosine_matrices_grains[mineral_i][grain_i][i][j])-1.0 <= 2.0 * std::numeric_limits<double>::epsilon(),
                           ExcMessage("3. a_cosine_matrices_grains[mineral_i][" + std::to_string(i) + "][" + std::to_string(j) +
                                      "] is larger than one: " + std::to_string(a_cosine_matrices_grains[mineral_i][grain_i][i][j]) + " (" + std::to_string(a_cosine_matrices_grains[mineral_i][grain_i][i][j]-1.0) + "). rotation_matrix = \n"
                                      + std::to_string(a_cosine_matrices_grains[mineral_i][grain_i][0][0]) + " " + std::to_string(a_cosine_matrices_grains[mineral_i][grain_i][0][1]) + " " + std::to_string(a_cosine_matrices_grains[mineral_i][grain_i][0][2]) + "\n"
                                      + std::to_string(a_cosine_matrices_grains[mineral_i][grain_i][1][0]) + " " + std::to_string(a_cosine_matrices_grains[mineral_i][grain_i][1][1]) + " " + std::to_string(a_cosine_matrices_grains[mineral_i][grain_i][1][2]) + "\n"
                                      + std::to_string(a_cosine_matrices_grains[mineral_i][grain_i][2][0]) + " " + std::to_string(a_cosine_matrices_grains[mineral_i][grain_i][2][1]) + " " + std::to_string(a_cosine_matrices_grains[mineral_i][grain_i][2][2])));
              }

            for (unsigned int grain_i = 0; grain_i < n_grains; ++grain_i)
              {
                Assert(isfinite(volume_fractions_grains[mineral_i][grain_i]),
                       ExcMessage("volume_fractions_grains[mineral_i][" + std::to_string(grain_i) + "] is not finite: "
                                  + std::to_string(volume_fractions_grains[mineral_i][grain_i]) + ", inv_sum_volume_grains[mineral_i] = "
                                  + std::to_string(inv_sum_volume_mineral) + "."));
              }

          }

        store_particle_data_extended(data_position,
                                     data,
                                     deformation_types,
                                     volume_fraction_mineral,
                                     volume_fractions_grains,
                                     a_cosine_matrices_grains,
                                     volume_fractions_grains_derivatives,
                                     a_cosine_matrices_grains_derivatives);



      }

      template<int dim>
      void
      LPO<dim>::orthogonalize_matrix(dealii::Tensor<2, 3> &tensor,
                                     double tolerance) const
      {
        if (std::abs(determinant(tensor) - 1.0) > tolerance)
          {
            LAPACKFullMatrix<double> identity_matrix(3);
            for (size_t i = 0; i < 3; i++)
              {
                identity_matrix.set(i,i,1);
              }

            FullMatrix<double> matrix_olivine(3);
            LAPACKFullMatrix<double> lapack_matrix_olivine(3);
            LAPACKFullMatrix<double> result(3);
            LAPACKFullMatrix<double> result2(3);

            // todo: find or add dealii functionallity to copy in one step.
            matrix_olivine.copy_from(tensor);
            lapack_matrix_olivine.copy_from(matrix_olivine);

            // now compute the svd of the matrices
            lapack_matrix_olivine.compute_svd();

            // Use the SVD results to orthogonalize: ((U*I)*V^T)^T
            lapack_matrix_olivine.get_svd_u().mmult(result,identity_matrix);
            result.mmult(result2,(lapack_matrix_olivine.get_svd_vt()));

            // todo: find or add dealii functionallity to copy in one step.
            matrix_olivine = result2;
            matrix_olivine.copy_to(tensor);
          }
      }

      template<int dim>
      DeformationType
      LPO<dim>::determine_deformation_type(const double stress, const double water_content) const
      {
        constexpr double MPa = 1e6;
        constexpr double ec_line_slope = -500./1050.;
        if (stress > (53)*MPa)
          {
            if (stress < (62)*MPa)
              {
                return DeformationType::OlivineDFabric;
              }
            else
              {
                return DeformationType::Passive;
              }
          }
        else
          {
            if (stress > (52)*MPa)
              {
                return DeformationType::OlivineCFabric;
              }
            else
              {
                if (stress > (33)*MPa)
                  {
                    return DeformationType::OlivineBFabric;
                  }
                else
                  {
                    return DeformationType::OlivineAFabric;
                  }
              }
          }
      }

      template<int dim>
      std::array<double,5>
      LPO<dim>::reference_resolved_shear_stress_from_deformation_type(DeformationType deformation_type, double max_value) const
      {
        std::array<double,5> ref_resolved_shear_stress;
        switch (deformation_type)
          {
            // from Kaminski and Ribe, GJI 2004.
            case DeformationType::OlivineAFabric :
              ref_resolved_shear_stress[0] = 1;
              ref_resolved_shear_stress[1] = 4;
              ref_resolved_shear_stress[2] = 6;
              ref_resolved_shear_stress[3] = 3;
              ref_resolved_shear_stress[4] = 6;
              ref_resolved_shear_stress[5] = 4;
              ref_resolved_shear_stress[6] = 6;
              ref_resolved_shear_stress[7] = 6;
              ref_resolved_shear_stress[8] = 30;
              ref_resolved_shear_stress[9] = 25;
              ref_resolved_shear_stress[10] = 25;
              break;

            // from Kaminski and Ribe, GJI 2004.
            case DeformationType::OlivineBFabric :
            ref_resolved_shear_stress[0] = 1;
              ref_resolved_shear_stress[1] = 4;
              ref_resolved_shear_stress[2] = 6;
              ref_resolved_shear_stress[3] = 3;
              ref_resolved_shear_stress[4] = 6;
              ref_resolved_shear_stress[5] = 4;
              ref_resolved_shear_stress[6] = 6;
              ref_resolved_shear_stress[7] = 6;
              ref_resolved_shear_stress[8] = 30;
              ref_resolved_shear_stress[9] = 25;
              ref_resolved_shear_stress[10] = 25;
              break;

            // from Kaminski and Ribe, GJI 2004.
            case DeformationType::OlivineCFabric :
              ref_resolved_shear_stress[0] = 1;
              ref_resolved_shear_stress[1] = 4;
              ref_resolved_shear_stress[2] = 6;
              ref_resolved_shear_stress[3] = 3;
              ref_resolved_shear_stress[4] = 6;
              ref_resolved_shear_stress[5] = 4;
              ref_resolved_shear_stress[6] = 6;
              ref_resolved_shear_stress[7] = 6;
              ref_resolved_shear_stress[8] = 30;
              ref_resolved_shear_stress[9] = 25;
              ref_resolved_shear_stress[10] = 25;
              break;

            // from Kaminski and Ribe, GRL 2002.
            case DeformationType::OlivineDFabric :
              ref_resolved_shear_stress[0] = 1;
              ref_resolved_shear_stress[1] = 4;
              ref_resolved_shear_stress[2] = 6;
              ref_resolved_shear_stress[3] = 3;
              ref_resolved_shear_stress[4] = 6;
              ref_resolved_shear_stress[5] = 4;
              ref_resolved_shear_stress[6] = 6;
              ref_resolved_shear_stress[7] = 6;
              ref_resolved_shear_stress[8] = 30;
              ref_resolved_shear_stress[9] = 25;
              ref_resolved_shear_stress[10] = 25;
              break;

            // Kaminski, Ribe and Browaeys, JGI, 2004 (same as in the matlab code)
      //      case DeformationType::OlivineEFabric :
      //        ref_resolved_shear_stress[0] = 2;
      //        ref_resolved_shear_stress[1] = 1;
     //         ref_resolved_shear_stress[2] = max_value;
     //         ref_resolved_shear_stress[3] = 3;
      //        break;

            // from Kaminski and Ribe, GJI 2004.
            // Todo: this one is not used in practice, since there is an optimalisation in
            // the code. So maybe remove it in the future.
           case DeformationType::Enstatite :
              ref_resolved_shear_stress[0] = 1;
              ref_resolved_shear_stress[1] = 4;
              ref_resolved_shear_stress[2] = 6;
              ref_resolved_shear_stress[3] = 3;
              ref_resolved_shear_stress[4] = 6;
              ref_resolved_shear_stress[5] = 4;
              ref_resolved_shear_stress[6] = 6;
              ref_resolved_shear_stress[7] = 6;
              ref_resolved_shear_stress[8] = 30;
              ref_resolved_shear_stress[9] = 25;
              ref_resolved_shear_stress[10] = 25;
              break;

            default:
              break;
          }
        return ref_resolved_shear_stress;
      }


      template<int dim>
      std::vector<std::vector<double> >
      LPO<dim>::volume_weighting(std::vector<double> fv, std::vector<std::vector<double>> angles) const
      {
        // Get volume weighted euler angles, using random draws to convert odf
        // to a discrete number of orientations, weighted by volume
        // 1a. Get index that would sort volume fractions AND
        //ix = np.argsort(fv[q,:]);
        // 1b. Get the sorted volume and angle arrays
        std::vector<double> fv_to_sort = fv;
        std::vector<double> fv_sorted = fv;
        std::vector<std::vector<double>> angles_sorted = angles;

        unsigned int n_grain = fv_to_sort.size();


        /**
         * ...
         */
        for (int i = n_grain-1; i >= 0; --i)
          {
            unsigned int index_max_fv = std::distance(fv_to_sort.begin(),max_element(fv_to_sort.begin(), fv_to_sort.end()));

            fv_sorted[i] = fv_to_sort[index_max_fv];
            angles_sorted[i] = angles[index_max_fv];
            Assert(angles[index_max_fv].size() == 3, ExcMessage("angles vector (size = " + std::to_string(angles[index_max_fv].size()) +
                                                                ") should have size 3."));
            Assert(angles_sorted[i].size() == 3, ExcMessage("angles_sorted vector (size = " + std::to_string(angles_sorted[i].size()) +
                                                            ") should have size 3."));
            fv_to_sort[index_max_fv] = -1;
          }


        // 2. Get cumulative weight for volume fraction
        std::vector<double> cum_weight(n_grains);
        std::partial_sum(fv_sorted.begin(),fv_sorted.end(),cum_weight.begin());
        // 3. Generate random indices
        boost::random::uniform_real_distribution<> dist(0, 1);
        std::vector<double> idxgrain(n_grains);
        for (unsigned int grain_i = 0; grain_i < n_grains; ++grain_i)
          {
            idxgrain[grain_i] = dist(this->random_number_generator); //random.rand(ngrains,1);
          }

        // 4. Find the maximum cum_weight that is less than the random value.
        // the euler angle index is +1. For example, if the idxGrain(g) < cumWeight(1),
        // the index should be 1 not zero)
        std::vector<std::vector<double>> angles_out(n_grains,std::vector<double>(3));
        for (unsigned int grain_i = 0; grain_i < n_grains; ++grain_i)
          {
            unsigned int counter = 0;
            for (unsigned int grain_j = 0; grain_j < n_grains-1; ++grain_j)
              {
                if (cum_weight[grain_j] < idxgrain[grain_i])
                  {
                    counter++;
                  }

                Assert(angles_sorted[counter].size() == 3, ExcMessage("angles_sorted vector (size = " + std::to_string(angles_sorted[counter].size()) +
                                                                      ") should have size 3."));
                angles_out[grain_i] = angles_sorted[counter];
                Assert(angles_out[counter].size() == 3, ExcMessage("angles_out vector (size = " + std::to_string(angles_out[counter].size()) +
                                                                   ") should have size 3."));
              }
          }
        return angles_out;
      }


      template <int dim>
      UpdateTimeFlags
      LPO<dim>::need_update() const
      {
        return update_time_step;
      }

      template <int dim>
      InitializationModeForLateParticles
      LPO<dim>::late_initialization_mode () const
      {
        return InitializationModeForLateParticles::initialize;
      }

      template <int dim>
      UpdateFlags
      LPO<dim>::get_needed_update_flags () const
      {
        return update_values | update_gradients;
      }

      template <int dim>
      std::vector<std::pair<std::string, unsigned int> >
      LPO<dim>::get_property_information() const
      {
        std::vector<std::pair<std::string,unsigned int> > property_information;

        for (size_t mineral_i = 0; mineral_i < n_minerals; mineral_i++)
          {
            property_information.push_back(std::make_pair("cpo mineral " + std::to_string(mineral_i) + " type",1));
            property_information.push_back(std::make_pair("cpo mineral " + std::to_string(mineral_i) + " volume fraction",1));
            for (unsigned int grain_i = 0; grain_i < n_grains; ++grain_i)
              {
                property_information.push_back(std::make_pair("cpo mineral " + std::to_string(mineral_i) + " grain " + std::to_string(grain_i) + " volume fraction",1));
                for (unsigned int index = 0; index < Tensor<2,3>::n_independent_components; index++)
                  {
                    property_information.push_back(std::make_pair("cpo mineral " + std::to_string(mineral_i) + " grain " + std::to_string(grain_i) + " a_cosine_matrix " + std::to_string(index),1));
                  }
              }
          }

        if (this->advection_method == AdvectionMethod::CrankNicolson)
          {
            for (size_t mineral_i = 0; mineral_i < n_minerals; mineral_i++)
              {
                for (unsigned int grain_i = 0; grain_i < n_grains; ++grain_i)
                  {
                    property_information.push_back(std::make_pair("cpo mineral " + std::to_string(mineral_i) + " grain " + std::to_string(grain_i) + " volume fraction derivative",1));
                    for (unsigned int index = 0; index < Tensor<2,3>::n_independent_components; index++)
                      {
                        property_information.push_back(std::make_pair("cpo mineral " + std::to_string(mineral_i) + " grain " + std::to_string(grain_i) + " a_cosine_matrix derivative" + std::to_string(index),1));
                      }
                  }
              }
          }
        return property_information;
      }

      template <int dim>
      double
      LPO<dim>::advect_forward_euler(std::vector<double> &volume_fractions,
                                     std::vector<Tensor<2,3> > &a_cosine_matrices,
                                     const std::pair<std::vector<double>, std::vector<Tensor<2,3> > > &derivatives,
                                     const double dt) const
      {
        double sum_volume_fractions = 0;
        for (unsigned int grain_i = 0; grain_i < a_cosine_matrices.size(); ++grain_i)
          {
            Assert(std::isfinite(volume_fractions[grain_i]),ExcMessage("volume_fractions[grain_i] is not finite before it is set."));
            volume_fractions[grain_i] = volume_fractions[grain_i] + dt * volume_fractions[grain_i] * derivatives.first[grain_i];
            Assert(std::isfinite(volume_fractions[grain_i]),ExcMessage("volume_fractions[grain_i] is not finite. grain_i = "
                                                                       + std::to_string(grain_i) + ", volume_fractions[grain_i] = " + std::to_string(volume_fractions[grain_i])
                                                                       + ", derivatives.first[grain_i] = " + std::to_string(derivatives.first[grain_i])));

            sum_volume_fractions += volume_fractions[grain_i];
          }


        for (unsigned int grain_i = 0; grain_i < n_grains; ++grain_i)
          {
            a_cosine_matrices[grain_i] = a_cosine_matrices[grain_i] + dt * a_cosine_matrices[grain_i] * derivatives.second[grain_i];
            // Should actually check them all, but this is the most important one.
            //Assert(a_cosine_matrices[grain_i][2][2] <= 1.0, ExcMessage("Internal error."));
          }

        Assert(sum_volume_fractions != 0, ExcMessage("Sum of volumes is equal to zero, which is not supporsed to happen."));
        return sum_volume_fractions;
      }


      template <int dim>
      double
      LPO<dim>::advect_backward_euler(std::vector<double> &volume_fractions,
                                      std::vector<Tensor<2,3> > &a_cosine_matrices,
                                      const std::pair<std::vector<double>, std::vector<Tensor<2,3> > > &derivatives,
                                      const double dt) const
      {
        double sum_volume_fractions = 0;
        for (unsigned int grain_i = 0; grain_i < n_grains; ++grain_i)
          {
            auto vf_old = volume_fractions[grain_i];
            auto vf_new = volume_fractions[grain_i];
            Assert(std::isfinite(vf_new),ExcMessage("vf_new is not finite before it is set."));
            for (size_t iteration = 0; iteration < property_advection_max_iterations; iteration++)
              {
                Assert(std::isfinite(vf_new),ExcMessage("vf_new is not finite before it is set in iteration " + std::to_string(iteration) + ": grain_i = "
                                                        + std::to_string(grain_i) + ", volume_fractions[grain_i] = " + std::to_string(volume_fractions[grain_i])
                                                        + ", derivatives.first[grain_i] = " + std::to_string(derivatives.first[grain_i]) + ", vf_new = " + std::to_string(vf_new)));

                vf_new = volume_fractions[grain_i] + dt * vf_new * derivatives.first[grain_i];

                Assert(std::isfinite(volume_fractions[grain_i]),ExcMessage("volume_fractions[grain_i] is not finite. grain_i = "
                                                                           + std::to_string(grain_i) + ", volume_fractions[grain_i] = " + std::to_string(volume_fractions[grain_i])
                                                                           + ", derivatives.first[grain_i] = " + std::to_string(derivatives.first[grain_i])+ ", vf_new = " + std::to_string(vf_new)));
                if (std::fabs(vf_new-vf_old) < property_advection_tolerance)
                  {
                    break;
                  }
                vf_old = vf_new;
              }

            volume_fractions[grain_i] = vf_new;
            Assert(std::isfinite(volume_fractions[grain_i]),ExcMessage("vf_new is not finite before it is set. grain_i = "
                                                                       + std::to_string(grain_i) + ", volume_fractions[grain_i] = " + std::to_string(volume_fractions[grain_i])
                                                                       + ", derivatives.first[grain_i] = " + std::to_string(derivatives.first[grain_i]) + ", vf_new = " + std::to_string(vf_new)));

            sum_volume_fractions += volume_fractions[grain_i];
          }

        for (unsigned int grain_i = 0; grain_i < n_grains; ++grain_i)
          {
            auto cosine_old = a_cosine_matrices[grain_i];
            auto cosine_new = a_cosine_matrices[grain_i];

            for (size_t iteration = 0; iteration < property_advection_max_iterations; iteration++)
              {
                cosine_new = a_cosine_matrices[grain_i] + dt * cosine_new * derivatives.second[grain_i];

                if ((cosine_new-cosine_old).norm() < property_advection_tolerance)
                  {
                    break;
                  }
                cosine_old = cosine_new;
              }

            a_cosine_matrices[grain_i] = cosine_new;
          }

        Assert(sum_volume_fractions != 0, ExcMessage("Sum of volumes is equal to zero, which is not supporsed to happen."));
        return sum_volume_fractions;
      }




      template <int dim>
      double
      LPO<dim>::advect_Crank_Nicolson(std::vector<double> &volume_fractions,
                                      std::vector<Tensor<2,3> > &a_cosine_matrices,
                                      const std::pair<std::vector<double>, std::vector<Tensor<2,3> > > &derivatives,
                                      std::vector<double> &previous_volume_fraction_derivatives,
                                      std::vector<Tensor<2,3> > &previous_a_cosine_matrices_derivatives,
                                      const double dt) const
      {
        double sum_volume_fractions = 0;
        for (unsigned int grain_i = 0; grain_i < n_grains; ++grain_i)
          {
            auto vf_old = volume_fractions[grain_i];
            auto vf_new = volume_fractions[grain_i];
            for (size_t iteration = 0; iteration < property_advection_max_iterations; iteration++)
              {
                vf_new = volume_fractions[grain_i] + dt * 0.5 * ((volume_fractions[grain_i] * previous_volume_fraction_derivatives[grain_i]) + (vf_new * derivatives.first[grain_i]));

                if (std::fabs(vf_new-vf_old) < property_advection_tolerance)
                  {
                    break;
                  }
                vf_old = vf_new;
              }

            previous_volume_fraction_derivatives[grain_i] = derivatives.first[grain_i];

            volume_fractions[grain_i] = vf_new;

            sum_volume_fractions += volume_fractions[grain_i];
          }

        for (unsigned int grain_i = 0; grain_i < n_grains; ++grain_i)
          {
            auto cosine_old = a_cosine_matrices[grain_i];
            auto cosine_new = a_cosine_matrices[grain_i];
            for (size_t iteration = 0; iteration < property_advection_max_iterations; iteration++)
              {
                cosine_new = a_cosine_matrices[grain_i] + dt * 0.5 * ((a_cosine_matrices[grain_i] * previous_a_cosine_matrices_derivatives[grain_i])
                                                                      + (cosine_new * derivatives.second[grain_i]));

                if ((cosine_new-cosine_old).norm() < property_advection_tolerance)
                  {
                    break;
                  }
                cosine_old = cosine_new;
              }

            previous_a_cosine_matrices_derivatives[grain_i] = derivatives.second[grain_i];

            a_cosine_matrices[grain_i] = cosine_new;
          }


        Assert(sum_volume_fractions != 0, ExcMessage("Sum of volumes is equal to zero, which is not supporsed to happen."));
        return sum_volume_fractions;
      }

      template <int dim>
      std::pair<std::vector<double>, std::vector<Tensor<2,3> > >
      LPO<dim>::compute_derivatives(const std::vector<double> &volume_fractions,
                                    const std::vector<Tensor<2,3> > &a_cosine_matrices,
                                    const SymmetricTensor<2,3> &strain_rate,
                                    const Tensor<2,3> &velocity_gradient_tensor,
                                    const double volume_fraction_mineral,
                                    const std::array<double,4> &ref_resolved_shear_stress) const
      {
        //
        std::vector<double> k_volume_fractions_zero = volume_fractions;
        std::vector<Tensor<2,3> > a_cosine_matrices_zero = a_cosine_matrices;
        std::pair<std::vector<double>, std::vector<Tensor<2,3> > > derivatives;
        switch (lpo_derivative_algorithm)
          {
            case LpoDerivativeAlgorithm::Zero:
            {
              return std::pair<std::vector<double>, std::vector<Tensor<2,3> > >(std::vector<double>(n_grains,0.0), std::vector<Tensor<2,3>>(n_grains, Tensor<2,3>()));
            }
            case LpoDerivativeAlgorithm::SpinTensor:
            {
              return this->compute_derivatives_spin_tensor(k_volume_fractions_zero,
                                                           a_cosine_matrices_zero,
                                                           strain_rate,
                                                           velocity_gradient_tensor,
                                                           volume_fraction_mineral,
                                                           ref_resolved_shear_stress);
              break;
            }
            case LpoDerivativeAlgorithm::DRex2004:
              return this->compute_derivatives_drex2004(k_volume_fractions_zero,
                                                        a_cosine_matrices_zero,
                                                        strain_rate,
                                                        velocity_gradient_tensor,
                                                        volume_fraction_mineral,
                                                        ref_resolved_shear_stress);
              break;
            default:
              AssertThrow(false, ExcMessage("Internal error."));
              break;
          }
        AssertThrow(false, ExcMessage("Internal error."));
        return derivatives;
      }

      template <int dim>
      std::pair<std::vector<double>, std::vector<Tensor<2,3> > >
      LPO<dim>::compute_derivatives_spin_tensor(const std::vector<double> &,
                                                const std::vector<Tensor<2,3> > &,
                                                const SymmetricTensor<2,3> &,
                                                const Tensor<2,3> &velocity_gradient_tensor,
                                                const double,
                                                const std::array<double,4> &) const
      {
        // dA/dt = W * A, where W is the spin tensor and A is the rotation matrix
        // The spin tensor is defined as W = 0.5 * ( L - L^T ), where L is the velocity gradient tensor.
        const Tensor<2,3> spin_tensor = -0.5 *(velocity_gradient_tensor - dealii::transpose(velocity_gradient_tensor));

        return std::pair<std::vector<double>, std::vector<Tensor<2,3> > >(std::vector<double>(n_grains,0.0), std::vector<Tensor<2,3>>(n_grains, spin_tensor));
      }

      template <int dim>
      std::pair<std::vector<double>, std::vector<Tensor<2,3> > >
      LPO<dim>::compute_derivatives_drex2004(const std::vector<double> &volume_fractions,
                                             const std::vector<Tensor<2,3> > &a_cosine_matrices,
                                             const SymmetricTensor<2,3> &strain_rate,
                                             const Tensor<2,3> &velocity_gradient_tensor,
                                             const double volume_fraction_mineral,
                                             const std::array<double,4> &ref_resolved_shear_stress,
                                             const bool prevent_nondimensionalization) const
      {
        SymmetricTensor<2,3> strain_rate_nondimensional = strain_rate;
        Tensor<2,3> velocity_gradient_tensor_nondimensional = velocity_gradient_tensor;
        // This if statement is only there for the unit test. In normal sitations it should always be set to false,
        // because the nondimensionalization should always be done (in this exact way), unless you really know what
        // you are doing.
        double nondimensionalization_value = 1.0;
        if (!prevent_nondimensionalization)
          {
            const std::array< double, 3 > eigenvalues = dealii::eigenvalues(strain_rate);
            nondimensionalization_value = std::max(std::abs(eigenvalues[0]),std::abs(eigenvalues[2]));

            Assert(!std::isnan(nondimensionalization_value), ExcMessage("The second invariant of the strain rate is not a number."));

            // Make the strain-rate and velocity gradient tensor non-dimensional
            // by dividing it through the second invariant
            if (nondimensionalization_value != 0)
              {
                strain_rate_nondimensional /= nondimensionalization_value;
                velocity_gradient_tensor_nondimensional /=  nondimensionalization_value;
              }
          }

        // create output variables
        std::vector<double> deriv_volume_fractions(n_grains);
        std::vector<Tensor<2,3> > deriv_a_cosine_matrices(n_grains);

        // create shorcuts
        const std::array<double, 4> &tau = ref_resolved_shear_stress;

        std::vector<double> strain_energy(n_grains);
        double mean_strain_energy = 0;

        // loop over grains
        for (unsigned int grain_i = 0; grain_i < n_grains; ++grain_i)
          {
            // Compute the Schmidt tensor for this grain (nu), s is the slip system.
            // We first compute beta_s,nu (equation 5, Kaminski & Ribe, 2001)
            // Then we use the beta to calculate the Schmidt tensor G_{ij} (Eq. 5, Kaminski & Ribe, 2001)
            Tensor<2,3> G;
            Tensor<1,3> w;
            Tensor<1,4> beta({1.0, 1.0, 1.0, 1.0});

            // these are variables we only need for olivine, but we need them for both
            // within this if bock and the next ones
            // todo: initialize to dealii uninitialized value
            unsigned int index_max_q = 0;
            unsigned int index_intermediate_q = 0;
            unsigned int index_min_q = 0;
            unsigned int index_inactive_q = 0;

            // compute G and beta
            Tensor<1,4> bigI;
            // this should be equal to a_cosine_matrices[grain_i]*a_cosine_matrices[grain_i]?
            // todo: check and maybe replace?
            for (unsigned int i = 0; i < 3; ++i)
              {
                for (unsigned int j = 0; j < 3; ++j)
                  {
                    bigI[0] = bigI[0] + strain_rate_nondimensional[i][j] * a_cosine_matrices[grain_i][0][i] * a_cosine_matrices[grain_i][1][j];

                    bigI[1] = bigI[1] + strain_rate_nondimensional[i][j] * a_cosine_matrices[grain_i][0][i] * a_cosine_matrices[grain_i][2][j];

                    bigI[2] = bigI[2] + strain_rate_nondimensional[i][j] * a_cosine_matrices[grain_i][2][i] * a_cosine_matrices[grain_i][1][j];

                    bigI[3] = bigI[3] + strain_rate_nondimensional[i][j] * a_cosine_matrices[grain_i][2][i] * a_cosine_matrices[grain_i][0][j];
                  }
              }

            if (bigI[0] == 0.0 && bigI[1] == 0.0 && bigI[2] == 0.0 && bigI[3] == 0.0)
              {
                // In this case there is no shear, only (posibily) a rotation. So \gamma_y and/or G should be zero.
                // Which is the default value, so do nothing.
              }
            else
              {
                // compute the element wise absolute value of the element wise
                // division of BigI by tau (tau = ref_resolved_shear_stress).
                std::vector<double> q_abs(4);
                for (unsigned int i = 0; i < 4; i++)
                  {
                    q_abs[i] = std::abs(bigI[i] / tau[i]);
                  }

                // here we find the indices starting at the largest value and ending at the smallest value
                // and assign them to special variables. Because all the variables are absolute values,
                // we can set them to a negative value to ignore them. This should be faster then deleting
                // the element, which would require allocation. (not tested)
                index_max_q = std::distance(q_abs.begin(),max_element(q_abs.begin(), q_abs.end()));

                q_abs[index_max_q] = -1;

                index_intermediate_q = std::distance(q_abs.begin(),max_element(q_abs.begin(), q_abs.end()));

                q_abs[index_intermediate_q] = -1;

                index_min_q = std::distance(q_abs.begin(),max_element(q_abs.begin(), q_abs.end()));

                q_abs[index_min_q] = -1;

                index_inactive_q = std::distance(q_abs.begin(),max_element(q_abs.begin(), q_abs.end()));

                // todo: explain
                Assert(bigI[index_max_q] != 0.0, ExcMessage("Internal error: bigI is zero."));
                double ratio = tau[index_max_q]/bigI[index_max_q];

                double q_intermediate = ratio * (bigI[index_intermediate_q]/tau[index_intermediate_q]);

                double q_min = ratio * (bigI[index_min_q]/tau[index_min_q]);

                beta[index_max_q] = 1.0; // max q_abs, weak system (most deformation) "s=1"
                beta[index_intermediate_q] = q_intermediate * std::pow(std::abs(q_intermediate), stress_exponent-1);
                beta[index_min_q] = q_min * std::pow(std::abs(q_min), stress_exponent-1);
                beta[index_inactive_q] = 0.0;

                for (unsigned int i = 0; i < 3; i++)
                  {
                    for (unsigned int j = 0; j < 3; j++)
                      {
                        G[i][j] = 2.0 * (beta[0] * a_cosine_matrices[grain_i][0][i] * a_cosine_matrices[grain_i][1][j]
                                         + beta[1] * a_cosine_matrices[grain_i][0][i] * a_cosine_matrices[grain_i][2][j]
                                         + beta[2] * a_cosine_matrices[grain_i][2][i] * a_cosine_matrices[grain_i][1][j]
                                         + beta[3] * a_cosine_matrices[grain_i][2][i] * a_cosine_matrices[grain_i][0][j]);
                      }
                  }
              }

            // Now calculate the analytic solution to the deformation minimization problem
            // compute gamma (equation 7, Kaminiski & Ribe, 2001)
            // todo: expand
            double top = 0;
            double bottom = 0;
            for (unsigned int i = 0; i < 3; ++i)
              {
                // Following the Drex code, which differs from EPSL paper,
                // which says gamma_nu depends on i+1: actually uses i+2
                unsigned int ip2 = i + 2;
                if (ip2 > 2)
                  ip2 = ip2-3;

                top = top - (velocity_gradient_tensor_nondimensional[i][ip2]-velocity_gradient_tensor_nondimensional[ip2][i])*(G[i][ip2]-G[ip2][i]);
                bottom = bottom - (G[i][ip2]-G[ip2][i])*(G[i][ip2]-G[ip2][i]);

                for (unsigned int j = 0; j < 3; ++j)
                  {
                    top = top + 2.0 * G[i][j]*velocity_gradient_tensor_nondimensional[i][j];
                    bottom = bottom + 2.0* G[i][j] * G[i][j];
                  }
              }
            // see comment on if all BigI are zero. In that case gamma should be zero.
            double gamma = bottom != 0.0 ? top/bottom : 0;

            // compute w (equation 8, Kaminiski & Ribe, 2001)
            // todo: explain what w is
            // todo: there was a loop around this in the phyton code, discuss/check
            w[0] = 0.5*(velocity_gradient_tensor_nondimensional[2][1]-velocity_gradient_tensor_nondimensional[1][2]) - 0.5*(G[2][1]-G[1][2])*gamma;
            w[1] = 0.5*(velocity_gradient_tensor_nondimensional[0][2]-velocity_gradient_tensor_nondimensional[2][0]) - 0.5*(G[0][2]-G[2][0])*gamma;
            w[2] = 0.5*(velocity_gradient_tensor_nondimensional[1][0]-velocity_gradient_tensor_nondimensional[0][1]) - 0.5*(G[1][0]-G[0][1])*gamma;

            // Compute strain energy for this grain (abrivated Estr)
            // For olivine: DREX only sums over 1-3. But Thissen's matlab code corrected
            // this and writes each term using the indices created when calculating bigI.
            // Note tau = RRSS = (tau_m^s/tau_o), this why we get tau^(p-n)
            const double rhos1 = std::pow(tau[index_max_q],exponent_p-stress_exponent) *
                                 std::pow(std::abs(gamma*beta[index_max_q]),exponent_p/stress_exponent);

            const double rhos2 = std::pow(tau[index_intermediate1_q],exponent_p-stress_exponent) *
                                 std::pow(std::abs(gamma*beta[index_intermediate1_q]),exponent_p/stress_exponent);

            const double rhos3 = std::pow(tau[index_intermediate2_q],exponent_p-stress_exponent) *
                                 std::pow(std::abs(gamma*beta[index_intermediate2_q]),exponent_p/stress_exponent);

            const double rhos4 = std::pow(tau[index_min_q],exponent_p-stress_exponent) *
                                 std::pow(std::abs(gamma*beta[index_intermediate3_q]),exponent_p/stress_exponent);


            const double rhos5 = std::pow(tau[index_inactive_q],exponent_p-stress_exponent) *
                                 std::pow(std::abs(gamma*beta[index_intermediate4_q]),exponent_p/stress_exponent);

            const double rhos6 = std::pow(tau[index_inactive_q],exponent_p-stress_exponent) *
                                 std::pow(std::abs(gamma*beta[index_intermediate5_q]),exponent_p/stress_exponent);

            const double rhos7 = std::pow(tau[index_inactive_q],exponent_p-stress_exponent) *
                                 std::pow(std::abs(gamma*beta[index_intermediate6_q]),exponent_p/stress_exponent);

            const double rhos8 = std::pow(tau[index_inactive_q],exponent_p-stress_exponent) *
                                 std::pow(std::abs(gamma*beta[index_intermediate7_q]),exponent_p/stress_exponent);

            const double rhos9 = std::pow(tau[index_inactive_q],exponent_p-stress_exponent) *
                                 std::pow(std::abs(gamma*beta[index_intermediate8_q]),exponent_p/stress_exponent);

            const double rhos10 = std::pow(tau[index_inactive_q],exponent_p-stress_exponent) *
                                 std::pow(std::abs(gamma*beta[index_min_q]),exponent_p/stress_exponent);

            const double rhos11 = std::pow(tau[index_inactive_q],exponent_p-stress_exponent) *
                                 std::pow(std::abs(gamma*beta[index_inactive_q]),exponent_p/stress_exponent);

            strain_energy[grain_i] = (rhos1 * exp(-nucleation_efficientcy * rhos1 * rhos1)
                                      + rhos2 * exp(-nucleation_efficientcy * rhos2 * rhos2)
                                      + rhos3 * exp(-nucleation_efficientcy * rhos3 * rhos3)
                                      + rhos4 * exp(-nucleation_efficientcy * rhos4 * rhos4)
                                      + rhos5 * exp(-nucleation_efficiency * rhos5 * rhos5)
                                      + rhos6 * exp(-nucleation_efficiency * rhos6 * rhos6)
                                      + rhos7 * exp(-nucleation_efficiency * rhos7 * rhos7)
                                      + rhos8 * exp(-nucleation_efficiency * rhos8 * rhos8)
                                      + rhos9 * exp(-nucleation_efficiency * rhos9 * rhos9)
                                      + rhos10 * exp(-nucleation_efficiency * rhos10 * rhos10)
                                      + rhos11 * exp(-nucleation_efficiency * rhos11 * rhos11));

            Assert(isfinite(strain_energy[grain_i]), ExcMessage("strain_energy[" + std::to_string(grain_i) + "] is not finite: " + std::to_string(strain_energy[grain_i])
                                                                + ", rhos1 = " + std::to_string(rhos1) + ", rhos2 = " + std::to_string(rhos2) + ", rhos3 = " + std::to_string(rhos3) + ", rhos4 = " + std::to_string(rhos4) 
                                                                + ",rhos5 = " + std::to_string(rhos5) + ", rhos6 = " + std::to_string(rhos6) + ", rhos7 = " + std::to_string(rhos7) + ", rhos8 = " + std::to_string(rhos8) 
                                                                + ", rhos9 = " + std::to_string(rhos9) + ", rhos10 = " + std::to_string(rhos10)
                                                                + ", rhos11= " + std::to_string(rhos11) + ", nucleation_efficientcy = " + std::to_string(nucleation_efficientcy) + "."));

            // compute the derivative of the cosine matrix a: \frac{\partial a_{ij}}{\partial t}
            // (Eq. 9, Kaminski & Ribe 2001)
            deriv_a_cosine_matrices[grain_i] = 0;
            if (volume_fractions[grain_i] >= threshold_GBS/n_grains)
              {
                deriv_a_cosine_matrices[grain_i] = permutation_operator_3d * w  * nondimensionalization_value;

                // volume averaged strain energy
                mean_strain_energy += volume_fractions[grain_i] * strain_energy[grain_i];

                Assert(isfinite(mean_strain_energy), ExcMessage("mean_strain_energy when adding grain " + std::to_string(grain_i) + " is not finite: " + std::to_string(mean_strain_energy)
                                                                + ", volume_fractions[grain_i] = " + std::to_string(volume_fractions[grain_i]) + ", strain_energy[grain_i] = " + std::to_string(strain_energy[grain_i]) + "."));
              }
            else
              {
                strain_energy[grain_i] = 0;
              }
          }

        // Change of volume fraction of grains by grain boundary migration
        for (unsigned int grain_i = 0; grain_i < n_grains; ++grain_i)
          {
            // Different than D-Rex. Here we actually only compute the derivative and do not multiply it with the volume_fractions. We do that when we advect.
            deriv_volume_fractions[grain_i] = volume_fraction_mineral * mobility * (mean_strain_energy - strain_energy[grain_i]) * nondimensionalization_value;

            Assert(isfinite(deriv_volume_fractions[grain_i]),
                   ExcMessage("deriv_volume_fractions[" + std::to_string(grain_i) + "] is not finite: "
                              + std::to_string(deriv_volume_fractions[grain_i])));
          }

        return std::pair<std::vector<double>, std::vector<Tensor<2,3> > >(deriv_volume_fractions, deriv_a_cosine_matrices);
      }

      template<int dim>
      unsigned int
      LPO<dim>::get_number_of_grains()
      {
        return n_grains;
      }

      template<int dim>
      unsigned int
      LPO<dim>::get_number_of_minerals()
      {
        return n_minerals;
      }


      template <int dim>
      void
      LPO<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Postprocess");
        {
          prm.enter_subsection("Particles");
          {
            prm.enter_subsection("LPO");
            {
              prm.declare_entry ("Random number seed", "1",
                                 Patterns::Integer (0),
                                 "The seed used to generate random numbers. This will make sure that "
                                 "results are reproducable as long as the problem is run with the "
                                 "same amount of MPI processes. It is implemented as final seed = "
                                 "user seed + MPI Rank. ");

              prm.declare_entry ("Number of grains per praticle", "50",
                                 Patterns::Integer (0),
                                 "The number of grains of olivine and the number of grain of enstatite "
                                 "each particle contains.");

              prm.declare_entry ("Property advection method", "Forward Euler",
                                 Patterns::Anything(),
                                 "Options: Forward Euler, Backward Euler, Crank-Nicolson");

              prm.declare_entry ("Property advection tolerance", "1e-10",
                                 Patterns::Double(0),
                                 "The Backward Euler and Crank-Nicolson property advection methods contain an iterations. "
                                 "This option allows for setting a tolerance. When the norm of tensor_new - tensor_old is "
                                 "smaller than this tolerance, the iteration is stopped.");

              prm.declare_entry ("Property advection max iterations", "100",
                                 Patterns::Integer(0),
                                 "The Backward Euler and Crank-Nicolson property advection methods contain an iterations. "
                                 "This option allows for setting the maximum ammount of iterations. Note that when the iteration "
                                 "is ended by the max iteration amount an assert is thrown.");

              prm.declare_entry ("LPO derivatives algorithm", "D-Rex 2004",
                                 Patterns::List(Patterns::Anything()),
                                 "Options: Spin tensor, D-Rex 2004");

              prm.enter_subsection("D-Rex 2004");
              {

                prm.declare_entry ("Minerals", "Olivine: Karato 2008, Enstatite",
                                   Patterns::List(Patterns::Anything()),
                                   "This determines what minerals and fabrics or fabric selectors are used used for the LPO calculation. "
                                   "The options are Olivine: A-fabric, Olivine: B-fabric, Olivine: C-fabric, Olivine: D-fabric, "
                                   "Olivine: E-fabric, Olivine: Karato 2008 or Enstatite. The "
                                   "Karato 2008 selector selects a fabric based on stress and water content as defined in "
                                   "figure 4 of the Karato 2008 review paper (doi: 10.1146/annurev.earth.36.031207.124120).");

                prm.declare_entry ("Mobility", "50",
                                   Patterns::Double(0),
                                   "The intrinsic grain boundary mobility for both olivine and enstatite. "
                                   "Todo: split for olivine and enstatite.");

                prm.declare_entry ("Volume fractions minerals", "0.5, 0.5",
                                   Patterns::List(Patterns::Double(0)),
                                   "The volume fraction for the different minerals. "
                                   "There need to be the same amount of values as there are minerals");

                prm.declare_entry ("Stress exponents", "3.5",
                                   Patterns::Double(0),
                                   "This is the power law exponent that characterizes the rheology of the "
                                   "slip systems. It is used in equation 11 of Kaminski et al., 2004. "
                                   "This is used for both olivine and enstatite. Todo: split?");

                prm.declare_entry ("Exponents p", "1.5",
                                   Patterns::Double(0),
                                   "This is exponent p as defined in equation 11 of Kaminski et al., 2004. ");

                prm.declare_entry ("Nucleation efficientcy", "5",
                                   Patterns::Double(0),
                                   "This is the dimensionless nucleation rate as defined in equation 8 of "
                                   "Kaminski et al., 2004. ");

                prm.declare_entry ("Threshold GBS", "0.3",
                                   Patterns::Double(0),
                                   "This is the grain-boundary sliding threshold. ");

                prm.declare_entry ("Use World Builder", "false",
                                   Patterns::Anything(),
                                   "Whether to use the world builder for setting the LPO.");
              }
              prm.leave_subsection();
            }
            prm.leave_subsection ();
          }
          prm.leave_subsection ();
        }
        prm.leave_subsection ();
      }


      template <int dim>
      void
      LPO<dim>::parse_parameters (ParameterHandler &prm)
      {

        prm.enter_subsection("Postprocess");
        {
          prm.enter_subsection("Particles");
          {
            prm.enter_subsection("LPO");
            {

              random_number_seed = prm.get_integer ("Random number seed");
              n_grains = prm.get_integer("Number of grains per praticle");

              property_advection_tolerance = prm.get_double("Property advection tolerance");
              property_advection_max_iterations = prm.get_integer ("Property advection max iterations");

              const std::string temp_lpo_derivative_algorithm = prm.get("LPO derivatives algorithm");

              if (temp_lpo_derivative_algorithm == "Spin tensor")
                {
                  lpo_derivative_algorithm = LpoDerivativeAlgorithm::SpinTensor;
                }
              else if (temp_lpo_derivative_algorithm ==  "D-Rex 2004")
                {
                  lpo_derivative_algorithm = LpoDerivativeAlgorithm::DRex2004;
                }
              else if (temp_lpo_derivative_algorithm ==  "Zero derivative")
                {
                  lpo_derivative_algorithm = LpoDerivativeAlgorithm::Zero;
                }
              else
                {
                  AssertThrow(false,
                              ExcMessage("The LPO derivatives algorithm needs to be one of the following: "
                                         "Spin tensor, D-Rex 2004."))
                }

              const std::string temp_advection_method = prm.get("Property advection method");
              if (temp_advection_method == "Forward Euler")
                {
                  advection_method = AdvectionMethod::ForwardEuler;
                }
              else if (temp_advection_method == "Backward Euler")
                {
                  advection_method = AdvectionMethod::BackwardEuler;
                }
              else if (temp_advection_method == "Crank-Nicolson")
                {
                  advection_method = AdvectionMethod::CrankNicolson;
                }
              else
                {
                  AssertThrow(false, ExcMessage("particle property advection method not found: \"" + temp_advection_method + "\""));
                }

              prm.enter_subsection("D-Rex 2004");
              {
                mobility = prm.get_double("Mobility");
                volume_fractions_minerals = Utilities::string_to_double(dealii::Utilities::split_string_list(prm.get("Volume fractions minerals")));
                stress_exponent = prm.get_double("Stress exponents");
                exponent_p = prm.get_double("Exponents p");
                nucleation_efficientcy = prm.get_double("Nucleation efficientcy");
                threshold_GBS = prm.get_double("Threshold GBS");
                use_world_builder = prm.get_bool("Use World Builder");

                const std::vector<std::string> temp_deformation_type_selector = dealii::Utilities::split_string_list(prm.get("Minerals"));
                n_minerals = temp_deformation_type_selector.size();
                deformation_type_selector.resize(n_minerals);

                for (size_t mineral_i = 0; mineral_i < n_minerals; mineral_i++)
                  {
                    if (temp_deformation_type_selector[mineral_i] == "Passive")
                      {
                        deformation_type_selector[mineral_i] = DeformationTypeSelector::Passive;
                      }
                    else if (temp_deformation_type_selector[mineral_i] == "Olivine: Karato 2008")
                      {
                        deformation_type_selector[mineral_i] = DeformationTypeSelector::OlivineKarato2008;
                      }
                    else if (temp_deformation_type_selector[mineral_i] ==  "Olivine: A-fabric")
                      {
                        deformation_type_selector[mineral_i] = DeformationTypeSelector::OlivineAFabric;
                      }
                    else if (temp_deformation_type_selector[mineral_i] ==  "Olivine: B-fabric")
                      {
                        deformation_type_selector[mineral_i] = DeformationTypeSelector::OlivineBFabric;
                      }
                    else if (temp_deformation_type_selector[mineral_i] ==  "Olivine: C-fabric")
                      {
                        deformation_type_selector[mineral_i] = DeformationTypeSelector::OlivineCFabric;
                      }
                    else if (temp_deformation_type_selector[mineral_i] ==  "Olivine: D-fabric")
                      {
                        deformation_type_selector[mineral_i] = DeformationTypeSelector::OlivineDFabric;
                      }
                   // else if (temp_deformation_type_selector[mineral_i] ==  "Enstatite: A-fabric")
                     // {
                     //   deformation_type_selector[mineral_i] = DeformationTypeSelector::EnstatiteBFabric;
                     // }
                    else if (temp_deformation_type_selector[mineral_i] ==  "Enstatite")
                      {
                        deformation_type_selector[mineral_i] = DeformationTypeSelector::Enstatite;
                      }
                    else
                      {
                        AssertThrow(false,
                                    ExcMessage("The  fabric needs to be one of the following: Olivine: Karato 2008, "
                                               "Olivine: A-fabric,Olivine: B-fabric,Olivine: C-fabric,Olivine: D-fabric,"
                                               "Olivine: E-fabric and Enstatite."))
                      }
                  }
              }
              prm.leave_subsection();


            }
            prm.leave_subsection ();
          }
          prm.leave_subsection ();
        }
        prm.leave_subsection ();


      }
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace Particle
  {
    namespace Property
    {
      ASPECT_REGISTER_PARTICLE_PROPERTY(LPO,
                                        "lpo",
                                        "A plugin in which the particle property tensor is "
                                        "defined as the deformation gradient tensor "
                                        "$\\mathbf F$ this particle has experienced. "
                                        "$\\mathbf F$ can be polar-decomposed into the left stretching tensor "
                                        "$\\mathbf L$ (the finite strain we are interested in), and the "
                                        "rotation tensor $\\mathbf Q$. See the corresponding cookbook in "
                                        "the manual for more detailed information.")
    }
  }
}
