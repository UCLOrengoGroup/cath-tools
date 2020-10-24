/// \file
/// \brief The view_cache_index_dim_linear class header

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

#ifndef _CATH_TOOLS_SOURCE_CT_UNI_CATH_STRUCTURE_VIEW_CACHE_INDEX_DETAIL_DIMS_DETAIL_VIEW_CACHE_INDEX_DIM_LINEAR_HPP
#define _CATH_TOOLS_SOURCE_CT_UNI_CATH_STRUCTURE_VIEW_CACHE_INDEX_DETAIL_DIMS_DETAIL_VIEW_CACHE_INDEX_DIM_LINEAR_HPP

#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/range/begin.hpp>
#include <boost/range/irange.hpp>
#include <boost/range/join.hpp>

#include "cath/common/algorithm/constexpr_clamp.hpp"
#include "cath/common/algorithm/constexpr_floor.hpp"
#include "cath/common/boost_addenda/range/indices.hpp"
#include "cath/common/cpp14/cbegin_cend.hpp"
#include "cath/common/debug_numeric_cast.hpp"
#include "cath/common/exception/invalid_argument_exception.hpp"
#include "cath/common/exception/not_implemented_exception.hpp"
#include "cath/common/size_t_literal.hpp"
#include "cath/structure/view_cache/index/detail/vcie_match_criteria.hpp"
#include "cath/structure/view_cache/index/view_cache_index_entry.hpp"

#include <cstddef>
#include <vector>

using namespace ::cath::common::literals;

namespace cath {
	namespace index {
		namespace detail {
			namespace detail {
				namespace detail {

					/// \brief TODOCUMENT
					template <typename T>
					inline constexpr int cell_index_of_value_in_current(const T   &prm_cell_width,   ///< TODOCUMENT
					                                                    const int &prm_start_offset, ///< TODOCUMENT
					                                                    const T   &prm_value         ///< TODOCUMENT
					                                                    ) {
						/// \todo Consider changing this to a numeric_cast (or debug_numeric_cast) if/when numeric_cast becomes constexpr
						///
						/// \todo Handle: libstdc++'s floor is constexpr but the standard doesn't define it as constexpr and
						///       libc++'s isn't. Find some constexpr workaround here. Is it adequate to just static_cast the fraction
						///       to int?
						return static_cast<int>( common::constexpr_floor( prm_value / prm_cell_width ) ) - prm_start_offset;
					}

					/// \brief TODOCUMENT
					template <typename T>
					inline constexpr size_t clamped_cell_index_of_value_in_current(const T      &prm_cell_width,   ///< TODOCUMENT
					                                                               const int    &prm_start_offset, ///< TODOCUMENT
					                                                               const size_t &prm_num_cells,    ///< TODOCUMENT
					                                                               const T      &prm_value         ///< TODOCUMENT
					                                                               ) {
						/// \todo Consider changing this to a numeric_cast (or debug_numeric_cast) if/when numeric_cast becomes constexpr
						return static_cast<size_t>( common::constexpr_clamp(
							cell_index_of_value_in_current( prm_cell_width, prm_start_offset, prm_value ),
							0,
							static_cast<int>( prm_num_cells )
						) );
					}

					/// \brief TODOCUMENT
					template <typename T>
					inline constexpr T min_value_in_cell_of_index_in_current(const T      &prm_cell_width,   ///< TODOCUMENT
					                                                         const int    &prm_start_offset, ///< TODOCUMENT
					                                                         const size_t &prm_cell_index    ///< TODOCUMENT
					                                                         ) {
						return cath::debug_numeric_cast<multiplier_type>(
							prm_start_offset + cath::debug_numeric_cast<int>( prm_cell_index )
						) * prm_cell_width;
					}

					/// \brief TODOCUMENT
					template <typename T>
					inline constexpr size_size_size_size_tpl search_cell_ranges(const T      &prm_cell_width,   ///< TODOCUMENT
					                                                            const int    &prm_start_offset, ///< TODOCUMENT
					                                                            const size_t &prm_num_cells,    ///< TODOCUMENT
					                                                            const T      &prm_search_begin, ///< TODOCUMENT
					                                                            const T      &prm_search_end    ///< TODOCUMENT
					                                                            ) {
						// The expressions for calculating the end values adds prm_cell_width to prm_search_end
						// so that the generate a one-past-the-end (in the standard C++/STL manner)
						return std::make_tuple(
							                                         clamped_cell_index_of_value_in_current( prm_cell_width, prm_start_offset, prm_num_cells, prm_search_begin                ),
							( prm_search_begin <  prm_search_end ) ? clamped_cell_index_of_value_in_current( prm_cell_width, prm_start_offset, prm_num_cells, prm_search_end + prm_cell_width ) : prm_num_cells,
							0_z,
							( prm_search_begin >= prm_search_end ) ? clamped_cell_index_of_value_in_current( prm_cell_width, prm_start_offset, prm_num_cells, prm_search_end + prm_cell_width ) : 0_z
						);
					}
				} // namespace detail

				/// \brief Perform indexing for a single layer based on a single linear, continuous dimension as specified in T
				///
				/// These classes (view_cache_index_dim_linear<>, view_cache_index_dim_dirn) don't hold the next levels down to the data being indexed.
				/// They simply specify how to manage that data, being passed a suitable vector as necessary.
				///
				/// The basic idea of this indexing on a continuous, linear dimension is to:
				///  * index things into one of a contiguous array of cells based on its value, and then
				///  * facilitate searching in the cells that may contain a match for a given value.
				template <typename T>
				class view_cache_index_dim_linear final {
				private:
					/// \brief TODOCUMENT
					using value_type = typename T::value_type;

					/// \brief The width of the cells in which the values will be stored
					value_type cell_width;

					/// \brief The offset (from 0) of the first cell,
					///        so the first cell is: [ cell_width * start_offset,  cell_width * ( start_offset + 1) )
					int start_offset;

					const value_type & get_cell_width() const;
					const int & get_start_offset() const;

					size_size_size_size_tpl search_cell_ranges(const size_t &,
					                                           const value_type &,
					                                           const value_type &) const;

					template <typename CELLS>
					bool has_cell_at_value(const CELLS &,
					                       const value_type &) const;

					template <typename CELLS>
					typename CELLS::value_type & cell_at_value(CELLS &,
					                                           const typename CELLS::value_type &,
					                                           const value_type &);

					template <typename CELLS>
					const typename CELLS::value_type & cell_at_value(const CELLS &,
					                                                 const value_type &) const;

				public:
					explicit view_cache_index_dim_linear(const value_type &);

					int cell_index_of_value_in_current(const value_type &) const;
					value_type min_value_in_cell_of_index_in_current(const size_t &) const;

					template <typename CELLS, typename DEFAULTS>
					void store(const view_cache_index_entry &,
					           CELLS &,
					           const DEFAULTS &);

					template <typename CELLS, typename ACTN>
					void perform_action_on_matches(const view_cache_index_entry &,
					                               const CELLS &,
					                               const vcie_match_criteria &,
					                               ACTN &) const;

					template <typename CELLS, typename ACTN>
					inline void perform_action_on_all_match_at_nodes(const CELLS &,
					                                                 const view_cache_index_dim_linear<T> &,
					                                                 const CELLS &,
					                                                 const vcie_match_criteria &,
					                                                 ACTN &) const;
				};

				/// \brief Getter for the cell width
				template <typename T>
				const typename view_cache_index_dim_linear<T>::value_type & view_cache_index_dim_linear<T>::get_cell_width() const {
					return cell_width;
				}

				/// \brief Getter for the start offset
				template <typename T>
				const int & view_cache_index_dim_linear<T>::get_start_offset() const {
					return start_offset;
				}

				/// \brief TODOCUMENT
				template <typename T>
				size_size_size_size_tpl view_cache_index_dim_linear<T>::search_cell_ranges(const size_t     &prm_num_cells,     ///< TODOCUMENT
				                                                                           const value_type &prm_search_begin, ///< TODOCUMENT
				                                                                           const value_type &prm_search_end    ///< TODOCUMENT
				                                                                           ) const {
					return detail::search_cell_ranges(
						cell_width,
						start_offset,
						prm_num_cells,
						prm_search_begin,
						prm_search_end
					);
				}

				/// \brief TODOCUMENT
				template <typename T>
				template <typename CELLS>
				bool view_cache_index_dim_linear<T>::has_cell_at_value(const CELLS      &prm_cells, ///< TODOCUMENT
				                                                       const value_type &prm_value  ///< TODOCUMENT
				                                                       ) const {
					const int cell_index_in_current = cell_index_of_value_in_current( prm_value );
					return ( cell_index_in_current > 0 && static_cast<size_t>( cell_index_in_current ) < prm_cells.size() );
				}

				/// \brief TODOCUMENT
				template <typename T>
				template <typename CELLS>
				typename CELLS::value_type & view_cache_index_dim_linear<T>::cell_at_value(CELLS                            &prm_cells,        ///< TODOCUMENT
				                                                                           const typename CELLS::value_type &prm_default_cell, ///< TODOCUMENT
				                                                                           const value_type                 &prm_value         ///< TODOCUMENT
				                                                                           ) {
					const int cell_index_in_current = cell_index_of_value_in_current( prm_value );
					if ( prm_cells.empty() ) {
						prm_cells.assign( 1, prm_default_cell );
						start_offset = cell_index_in_current;
						return prm_cells.front();
					}
					if ( cell_index_in_current < 0 ) {
						const size_t num_to_prepend = cath::debug_numeric_cast<size_t>( -cell_index_in_current );
						const int new_start_offest = get_start_offset() + cell_index_in_current;
						// Come GCC v4.9 above, replace this begin() with cbegin()
						prm_cells.insert( std::begin( prm_cells ), num_to_prepend, prm_default_cell );
						start_offset = new_start_offest;
						return prm_cells.front();
					}
					if ( static_cast<size_t>( cell_index_in_current ) >= prm_cells.size() ) {
						const size_t new_size    = 1 + cath::debug_numeric_cast<size_t>( cell_index_in_current );
						prm_cells.resize( new_size, prm_default_cell );
						return prm_cells.back();
					}
					return prm_cells[ cath::debug_numeric_cast<size_t>( cell_index_in_current ) ];
				}

				/// \brief TODOCUMENT
				template <typename T>
				template <typename CELLS>
				const typename CELLS::value_type & view_cache_index_dim_linear<T>::cell_at_value(const CELLS      &prm_cells, ///< TODOCUMENT
				                                                                                 const value_type &prm_value  ///< TODOCUMENT
				                                                                                 ) const {
	#ifndef NDEBUG
					if ( prm_cells.empty() ) {
						BOOST_THROW_EXCEPTION(cath::common::invalid_argument_exception("Cannot get entry at_value() with no populated cells"));
					}
					if ( ! has_cell_at_value( prm_cells, prm_value ) ) {
						BOOST_THROW_EXCEPTION(cath::common::invalid_argument_exception("Value is not found in any of the cells"));
					}
	#endif
					return prm_cells[ cell_index_of_value_in_current( prm_value ) ];
				}

				/// \brief TODOCUMENT
				template <typename T>
				view_cache_index_dim_linear<T>::view_cache_index_dim_linear(const value_type &prm_cell_width ///< TODOCUMENT
				                                                            ) : cell_width   ( prm_cell_width ),
				                                                                start_offset ( 0              ) {
					T().check_cell_width( get_cell_width() );
				}

				/// \brief The cell index of the specified value in the current state of the index
				template <typename T>
				int view_cache_index_dim_linear<T>::cell_index_of_value_in_current(const value_type &prm_value ///< TODOCUMENT
				                                                                   ) const {
					return detail::cell_index_of_value_in_current( cell_width, start_offset, prm_value );
				}

				/// \brief TODOCUMENT
				template <typename T>
				typename view_cache_index_dim_linear<T>::value_type view_cache_index_dim_linear<T>::min_value_in_cell_of_index_in_current(const size_t &prm_cell_index ///< TODOCUMENT
				                                                                                                                          ) const {
					return detail::min_value_in_cell_of_index_in_current( cell_width, start_offset, prm_cell_index );
				}

				/// \brief TODOCUMENT
				template <typename T>
				template <typename CELLS, typename DEFAULTS>
				void view_cache_index_dim_linear<T>::store(const view_cache_index_entry &prm_entry,   ///< TODOCUMENT
				                                           CELLS                        &prm_cells,   ///< TODOCUMENT
				                                           const DEFAULTS               &prm_defaults ///< TODOCUMENT
				                                           ) {
					// Grab the value in question and TODOCUMENT
					const value_type &value = T().get_index_value( prm_entry );
					cell_at_value( prm_cells, typename CELLS::value_type{ prm_defaults.get_head() }, value ).store(
						prm_entry,
						prm_defaults.get_tail()
					);
				}

				/// \brief TODOCUMENT
				template <typename T>
				template <typename CELLS, typename ACTN>
				void view_cache_index_dim_linear<T>::perform_action_on_matches(const view_cache_index_entry &prm_entry,    ///< TODOCUMENT
				                                                               const CELLS                  &prm_cells,    ///< TODOCUMENT
				                                                               const vcie_match_criteria    &prm_criteria, ///< TODOCUMENT
				                                                               ACTN                         &prm_action    ///< TODOCUMENT
				                                                               ) const {
					const value_type &value         = T().get_index_value  ( prm_entry    );
					const value_type &search_radius = T().get_search_radius( prm_criteria );
					const size_t      num_cells     = prm_cells.size();

					if ( ! prm_cells.empty() ) {
						const auto prepped_begin_and_end = T().prepare_search_begin_and_end(
							value - search_radius,
							value + search_radius
						);
						const auto search_ranges = search_cell_ranges( num_cells, prepped_begin_and_end.first, prepped_begin_and_end.second );
						// std::cerr << "Search ["         << value - search_radius
						//           << ", "               << value + search_radius
						//           << "] in "            << num_cells
						//           << " cells of width " << cell_width
						//           << " from "           << cath::debug_numeric_cast<multiplier_type>( start_offset             ) * cell_width
						//           << " to "             << cath::debug_numeric_cast<multiplier_type>( start_offset + num_cells ) * cell_width
						//           << " - [ "            << std::get<0>( search_ranges )
						//           << " ("               << min_value_in_cell_of_index_in_current( std::get<0>( search_ranges ) )
						//           << "), "              << std::get<1>( search_ranges )
						//           << " ("               << min_value_in_cell_of_index_in_current( std::get<1>( search_ranges ) )
						//           << ");  "              << std::get<2>( search_ranges )
						//           << " ("               << min_value_in_cell_of_index_in_current( std::get<2>( search_ranges ) )
						//           << "), "               << std::get<3>( search_ranges )
						//           << " ("               << min_value_in_cell_of_index_in_current( std::get<3>( search_ranges ) )
						//           << ") ]"               << std::endl;

						const auto joined_iranges = boost::range::join(
							boost::irange( std::get<0>( search_ranges ), std::get<1>( search_ranges ) ),
							boost::irange( std::get<2>( search_ranges ), std::get<3>( search_ranges ) )
						);
						for (const size_t &search_cell_ctr : joined_iranges) {
							prm_cells[ search_cell_ctr ].perform_action_on_matches( prm_entry, prm_criteria, prm_action );
						}
					}
				}

				/// \brief TODOCUMENT
				template <typename T>
				template <typename CELLS, typename ACTN>
				inline void view_cache_index_dim_linear<T>::perform_action_on_all_match_at_nodes(const CELLS                          &prm_query_cells,     ///< TODOCUMENT
				                                                                                 const view_cache_index_dim_linear<T> &prm_match_dimension, ///< TODOCUMENT
				                                                                                 const CELLS                          &prm_match_cells,     ///< TODOCUMENT
				                                                                                 const vcie_match_criteria            &prm_criteria,        ///< TODOCUMENT
				                                                                                 ACTN                                 &prm_action           ///< TODOCUMENT
				                                                                                 ) const {
					const value_type &search_radius   = T().get_search_radius( prm_criteria );


					if ( ! prm_query_cells.empty() && ! prm_match_cells.empty() ) {
						const size_t &num_query_cells = prm_query_cells.size();

						for (const size_t &cell_ctr : common::indices( num_query_cells ) ) {
							const auto &query_cell = prm_query_cells[ cell_ctr ];
							if ( ! query_cell.empty() ) {

								const auto prepped_begin_and_end = T().prepare_search_begin_and_end(
									min_value_in_cell_of_index_in_current( cell_ctr     ) - search_radius,
									min_value_in_cell_of_index_in_current( cell_ctr + 1 ) + search_radius
								);

								const auto search_ranges = prm_match_dimension.search_cell_ranges(
									prm_match_cells.size(),
									prepped_begin_and_end.first,
									prepped_begin_and_end.second
								);
								const auto joined_iranges = boost::range::join(
									boost::irange( std::get<0>( search_ranges ), std::get<1>( search_ranges ) ),
									boost::irange( std::get<2>( search_ranges ), std::get<3>( search_ranges ) )
								);
								for (const size_t &match_cell_ctr : joined_iranges) {
									const auto &match_cell = prm_match_cells[ match_cell_ctr ];
									if ( ! match_cell.empty() ) {
										// std::cerr << std::string( 10 - CELLS::value_type::num_dims_remaining, ' ' )
										//           << "In dimension "                          << T().get_name()
										//           << ", scanning for matches to query cell [" << min_value_in_cell_of_index_in_current( cell_ctr     )
										//           << ", "                                     << min_value_in_cell_of_index_in_current( cell_ctr + 1 )
										//           << "].size_"                                << query_cell.get_num_cells()
										//           << ", trying match cell ["                  << prm_match_dimension.min_value_in_cell_of_index_in_current( match_cell_ctr     )
										//           << ", "                                     << prm_match_dimension.min_value_in_cell_of_index_in_current( match_cell_ctr + 1 )
										//           << "].size_"                                << match_cell.get_num_cells()
										//           << std::endl;

										query_cell.perform_action_on_all_match_at_nodes( match_cell, prm_criteria, prm_action );
									}
								}
							}
						}
					}
				}

			} // namespace detail
		} // namespace detail
	} // namespace index
} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_UNI_CATH_STRUCTURE_VIEW_CACHE_INDEX_DETAIL_DIMS_DETAIL_VIEW_CACHE_INDEX_DIM_LINEAR_HPP
