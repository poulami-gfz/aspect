/*
  Copyright (C) 2020 by the authors of the ASPECT code.

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

#include <aspect/postprocess/visualization/principal_strain.h>

#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>

namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      template <int dim>
      PrincipalStrain<dim>::
      PrincipalStrain ()
        :
        DataPostprocessor<dim> ()
      {}


     
     template <int dim>
      std::vector<std::string>
      PrincipalStrain<dim>::
      get_names () const
      {
        std::vector<std::string> solution_names;

        // dim principal strain values
        for (unsigned int i=0; i<dim; ++i)
          solution_names.push_back("principal_strain_" + std::to_string(i+1));

        // dim principal strain directions
        for (unsigned int i=0; i<dim; ++i)
          for (unsigned int j=0; j<dim; ++j)
            solution_names.push_back("principal_strain_direction_" + std::to_string(i+1));

        return solution_names;
      }



      template <int dim>
      std::vector<DataComponentInterpretation::DataComponentInterpretation>
      PrincipalStrain<dim>::
      get_data_component_interpretation () const
      {
        std::vector<DataComponentInterpretation::DataComponentInterpretation> solution_components;

        // dim principal strain values
        for (unsigned int i=0; i<dim; ++i)
          solution_components.push_back(DataComponentInterpretation::component_is_scalar);

        // dim principal strain directions
        for (unsigned int i=0; i<dim; ++i)
          for (unsigned int j=0; j<dim; ++j)
            solution_components.push_back(DataComponentInterpretation::component_is_part_of_vector);

        return solution_components;
      }



      template <int dim>
      UpdateFlags
      PrincipalStrain<dim>::
      get_needed_update_flags () const
      {
        return update_values | update_gradients | update_quadrature_points;
      }



      template <int dim>
      void
      PrincipalStrain<dim>::
      evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &input_data,
                            std::vector<Vector<double> > &computed_quantities) const
      {
//         const std::vector<double>solution_values= input_data.solution_values;
        const unsigned int n_quadrature_points = input_data.solution_values.size();
        Assert (computed_quantities.size() == n_quadrature_points, ExcInternalError());
        Assert (computed_quantities[0].size() == dim*dim + dim, ExcInternalError());
        Assert (input_data.solution_values[0].size() == this->introspection().n_components,   ExcInternalError());
        Assert (input_data.solution_gradients[0].size() == this->introspection().n_components,  ExcInternalError());  
        Assert(this->introspection().compositional_name_exists("strain_xx"),
                          ExcMessage("strain_xx <" +
                                     field_name +
                                     "> exists for which you want to visualize the principal strain."));
        Assert(this->introspection().compositional_index_for_name("strain_xx") == 0, ExcMessage("The strain components should be the first compositional fields."))

//                 // Set use_strain_rates to true since the compaction viscosity might also depend on the strain rate.
//         MaterialModel::MaterialModelInputs<dim> in(input_data,
//                                                    this->introspection());
//         MaterialModel::MaterialModelOutputs<dim> out(n_quadrature_points,
//                                                      this->n_compositional_fields());
// 
//         this->get_material_model().evaluate(in, out);
        
          // Assign the strain components to the compositional fields reaction terms.
          // If there are too many fields, we simply fill only the first fields with the
          // existing strain rate tensor components.
          for (unsigned int q=0; q < n_quadrature_points; ++q)
            {
              // Convert the compositional fields into the tensor quantity they represent.
              Tensor<2,dim> strain;
              for (unsigned int i = 0; i < Tensor<2,dim>::n_independent_components ; ++i)
              {
                strain[Tensor<2,dim>::unrolled_to_component_indices(i)] = input_data.solution_values[q][this->introspection().component_indices.compositional_fields[i]];
              }

            const std::array<std::pair<double, Tensor<1,dim>>, dim> principal_straines_and_directions = eigenvectors(symmetrize(strain));

            // dim principal strain values
            for (unsigned int i=0; i<dim; ++i)
              {
              computed_quantities[q][i] = principal_straines_and_directions[i].first;
              }

            // dim principal strain directions
            for (unsigned int i=0; i<dim; ++i){
            {
              for (unsigned int j=0; j<dim; ++j)
              {
                computed_quantities[q][dim + i*dim + j] = principal_straines_and_directions[i].second[j];  
              }                           
            }
          }
        }    
      }
    

          
//       template <int dim>
//       void
//       PrincipalStrain<dim>::
//       declare_parameters (ParameterHandler &prm)
//       {
//       }

      template <int dim>
      void
      PrincipalStrain<dim>::parse_parameters (ParameterHandler &prm)
      {
            prm.enter_subsection("Compositional fields");
            {
              const std::string field_name = prm.get("Names of fields");

////             AssertThrow(this->introspection().compositional_name_exists("strain_xx"),
 ///                         ExcMessage("No compositional field with name <" +
//                                     "strain_xx" +
//                                     "> exists for which you want to visualize the principal strain."));

//               int compositional_field = this->introspection().compositional_index_for_name("strain_xx");
            }
            prm.leave_subsection();
      }          
        
        // average the values if requested
      // const auto &viz = this->get_postprocess_manager().template get_matching_postprocessor<Postprocess::Visualization<dim> >();
       //if (!viz.output_pointwise_stress_and_strain())
        // average_quantities(computed_quantities);   
      }


    }
    
  }



// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(PrincipalStrain,
                                                  "principal strain",
                                                  "A visualization output object that outputs the "
                                                  "principal strain values and directions, i.e., the "
                                                  "eigenvalues and eigenvectors of the strain tensor. "
                                                  "The postprocessor can either operate on the full "
                                                  "strain tensor or only on the deviatoric strain tensor.")
    }
  }
}
