/// \file
/// \brief The view_cache_index_layer class header

/// \copyright
/// CATH Tools - Protein structure comparison tools such as SSAP and SNAP
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

#ifndef _CATH_TOOLS_SOURCE_UNI_STRUCTURE_VIEW_CACHE_INDEX_DETAIL_SCAFFOLD_VIEW_CACHE_INDEX_LAYER_H
#define _CATH_TOOLS_SOURCE_UNI_STRUCTURE_VIEW_CACHE_INDEX_DETAIL_SCAFFOLD_VIEW_CACHE_INDEX_LAYER_H

#include <boost/core/demangle.hpp>
#include <boost/log/trivial.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/range/begin.hpp>

#include "cath/common/debug_numeric_cast.hpp"
#include "cath/common/exception/invalid_argument_exception.hpp"
#include "cath/structure/view_cache/index/view_cache_index_entry.hpp"

#include <cstddef>
#include <vector>

namespace cath {
	namespace index {
		namespace detail {

			/// \brief Manage one layer of the view_cache_index that indexes on the dimension DIM (eg view_cache_index_dim_linear_x)
			///
			/// T is the type of the next layer down (or view_cache_index_tail if this is the last layer).
			///
			/// The layer uses DIM to store view_cache_index_entry objects in bins and then when later searching
			/// for similar entries to some query view_cache_index_entry, the DIM can identify the bins which
			/// might possible contain matches so that all others can safely be ignored.
			///
			/// For example, view_cache_index_dim_linear_x puts view_cache_index_entry objects in bins according to the
			/// x component of their view. Then a search with a given query view_cache_index_entry only needs to look in
			/// bins covering ranges of x values that might be close enough to the query's to be interesting; all other
			/// bins can be skipped.
			template <typename DIM, typename T>
			class view_cache_index_layer final {
			private:
				/// \brief An instance of the dimension with which this layer is organising Ts
				DIM the_dimension;

				/// \brief The cells into which the Ts are indexed by DIM
				std::vector<T> cells;

				int cell_index_of_value_in_current(const double &) const;

			public:
				explicit view_cache_index_layer(const DIM &);

				// const double & get_cell_width() const;
				bool empty() const;
				size_t get_num_cells() const;
	
				// const T & get_cell_entry(const size_t &) const;
				// bool has_cell_at_value(const double &) const;
				// T & cell_at_value(const double &);
				// const T & cell_at_value(const double &) const;

				template <typename DEFAULTS>
				void store(const view_cache_index_entry &,
				           const DEFAULTS &);

				template <typename ACTN>
				void perform_action_on_matches(const view_cache_index_entry &,
				                               const vcie_match_criteria &,
				                               ACTN &) const;

				template <typename ACTN>
				void perform_action_on_all_match_at_leaves(ACTN &) const;


				template <typename ACTN>
				void perform_action_on_all_match_at_nodes(const view_cache_index_layer<DIM, T> &,
				                                          const vcie_match_criteria &,
				                                          ACTN &) const;

				/// \brief TODOCUMENT
				constexpr static size_t num_dims_remaining = T::num_dims_remaining + 1;
			};

			/// \brief Ctor to build a layer from a DIM 
			template <typename DIM, typename T>
			view_cache_index_layer<DIM, T>::view_cache_index_layer(const DIM &prm_dim ///< The DIM to use to index values in this layer
			                                                       ) : the_dimension( prm_dim ) {
			}

			// /// \brief TODOCUMENT
			// template <typename DIM, typename T>
			// const double & view_cache_index_layer<DIM, T>::get_cell_width() const {
			// 	return the_dimension.get_cell_width();
			// }

			/// \brief TODOCUMENT
			template <typename DIM, typename T>
			bool view_cache_index_layer<DIM, T>::empty() const {
				return cells.empty();
			}
	
			/// \brief TODOCUMENT
			template <typename DIM, typename T>
			size_t view_cache_index_layer<DIM, T>::get_num_cells() const {
				return cells.size();
			}
	
			// /// \brief TODOCUMENT
			// template <typename DIM, typename T>
			// const T & view_cache_index_layer<DIM, T>::get_cell_entry(const size_t &prm_index ///< TODOCUMENT
			//                                                        ) const {
			// 	return cells[ prm_index ];
			// }
	
			// /// \brief TODOCUMENT
			// template <typename DIM, typename T>
			// bool view_cache_index_layer<DIM, T>::has_cell_at_value(const double &prm_value ///< TODOCUMENT
			//                                                      ) const {
			// 	return the_dimension.has_cell_at_value( prm_value );
			// }
	
			// /// \brief TODOCUMENT
			// template <typename DIM, typename T>
			// T & view_cache_index_layer<DIM, T>::cell_at_value(const double &prm_value ///< TODOCUMENT
			//                                                 ) {
			// 	return the_dimension.cell_at_value( cells, prm_value );
			// }
	
			// /// \brief TODOCUMENT
			// template <typename DIM, typename T>
			// const T & view_cache_index_layer<DIM, T>::cell_at_value(const double &prm_value ///< TODOCUMENT
			//                                                       ) const {
			// 	return the_dimension.cell_at_value( cells, prm_value );
			// }


			/// \brief Store a view_cache_index_entry in this layer, using 
			template <typename DIM, typename T>
			template <typename DEFAULTS>
			void view_cache_index_layer<DIM, T>::store(const view_cache_index_entry &prm_entry,   ///< TODOCUMENT
			                                           const DEFAULTS               &prm_defaults ///< TODOCUMENT
			                                           ) {
				the_dimension.store( prm_entry, cells, prm_defaults );
			}

			/// \brief TODOCUMENT
			template <typename DIM, typename T>
			template <typename ACTN>
			void view_cache_index_layer<DIM, T>::perform_action_on_matches(const view_cache_index_entry &prm_entry,    ///< TODOCUMENT
			                                                               const vcie_match_criteria    &prm_criteria, ///< TODOCUMENT
			                                                               ACTN                         &prm_action    ///< TODOCUMENT
			                                                               ) const {
				the_dimension.perform_action_on_matches( prm_entry, cells, prm_criteria, prm_action );
			}

			/// \brief TODOCUMENT
			template <typename DIM, typename T>
			template <typename ACTN>
			void view_cache_index_layer<DIM, T>::perform_action_on_all_match_at_leaves(ACTN &prm_action ///< TODOCUMENT
			                                                                           ) const {
				BOOST_LOG_TRIVIAL( trace ) << ::boost::core::demangle( typeid( *this ).name() ) << " : " << cells.size();
				for (const T &cell : cells) {
					cell.perform_action_on_all_match_at_leaves( prm_action );
				}
			}


			/// \brief TODOCUMENT
			template <typename DIM, typename T>
			template <typename ACTN>
			void view_cache_index_layer<DIM, T>::perform_action_on_all_match_at_nodes(const view_cache_index_layer<DIM, T> &prm_match_layer, ///< TODOCUMENT
			                                                                          const vcie_match_criteria            &prm_criteria,    ///< TODOCUMENT
			                                                                          ACTN                                 &prm_action       ///< TODOCUMENT
			                                                                          ) const {
				the_dimension.perform_action_on_all_match_at_nodes(
					cells,
					prm_match_layer.the_dimension,
					prm_match_layer.cells,
					prm_criteria,
					prm_action
				);
			}

		} // namespace detail
	} // namespace index
} // namespace cath

#endif

