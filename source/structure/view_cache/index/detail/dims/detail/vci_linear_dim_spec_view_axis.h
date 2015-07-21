/// \file
/// \brief The vci_linear_dim_spec_view_axis class header

/// \copyright
/// CATH Binaries - Protein structure comparison tools such as SSAP and SNAP
/// Copyright (C) 2011, Orengo Group, University College London
///
/// This program is free software: you can redistribute it and/or modify
/// it under the terms of the GNU General Public License as published by
/// the Free Software Foundation, either version 3 of the License, or
/// (at your option) any later version.
///
/// This program is distributed in the hope that it will be useful,
/// but WITHOUT ANY WARRANTY; without even the implied warranty of
/// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
/// GNU General Public License for more details.
///
/// You should have received a copy of the GNU General Public License
/// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#ifndef VCI_LINEAR_DIM_SPEC_VIEW_AXIS_H_INCLUDED
#define VCI_LINEAR_DIM_SPEC_VIEW_AXIS_H_INCLUDED

#include "structure/view_cache/index/detail/view_cache_index_type_aliases.h"

namespace cath {
	namespace index {
		namespace detail {
			namespace detail {

				/// \brief TODOCUMENT
				template <typename F>
				struct vci_linear_dim_spec_view_axis final {
					/// \brief TODOCUMENT
					using value_type      = view_base_type;

					/// \brief TODOCUMENT
					using value_pair      = std::pair<value_type, value_type>;

					/// \brief TODOCUMENT
					using linear_dim_type = view_cache_index_dim_linear<vci_linear_dim_spec_view_axis>;

					/// \brief TODOCUMENT
					std::string get_name() const {
						return F().get_name();
					}

					/// \brief TODOCUMENT
					void check_cell_width(const value_type &arg_cell_width ///< TODOCUMENT
					                      ) {
						if ( ! boost::math::isnormal( arg_cell_width ) || arg_cell_width <= 0.0 ) {
							BOOST_THROW_EXCEPTION(cath::common::invalid_argument_exception("Cannot create a view_cache_index_dim_linear<vci_linear_dim_spec_view_axis<>> with a cell_width that isn't a sensible, strictly-positive value"));
						}
					}

					/// \brief TODOCUMENT
					inline value_pair prepare_search_begin_and_end(const value_type &arg_search_min, ///< TODOCUMENT
					                                               const value_type &arg_search_max  ///< TODOCUMENT
					                                               ) {
						// std::cerr << "In vci_linear_dim_spec_view_axis::prepare_search_begin_and_end(), arg_search_min is " << arg_search_min << std::endl;
						// std::cerr << "In vci_linear_dim_spec_view_axis::prepare_search_begin_and_end(), arg_search_max is " << arg_search_max << std::endl;
						// return std::make_tuple(
						// 	true,
						// 	arg_linear_dim.cell_index_of_value_in_current( arg_search_min ),
						// 	arg_linear_dim.cell_index_of_value_in_current( arg_search_max ) + 1
						// );
						assert( arg_search_min < arg_search_max );
						return std::make_pair( arg_search_min, arg_search_max );
					}

					/// \brief TODOCUMENT
					value_type get_search_radius(const vcie_match_criteria &arg_criteria ///< TODOCUMENT
					                             ) {
						return std::sqrt( arg_criteria.get_maximum_squared_distance() );
					}

					/// \brief TODOCUMENT
					const value_type & get_index_value(const view_cache_index_entry &arg_entry ///< TODOCUMENT
					                                   ) {
						return F()( arg_entry );
					}
				};

			}
		}
	}
}

#endif

