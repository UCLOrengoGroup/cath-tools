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

#ifndef VIEW_CACHE_INDEX_DIM_LINEAR_H_INCLUDED
#define VIEW_CACHE_INDEX_DIM_LINEAR_H_INCLUDED

#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/range/begin.hpp>
#include <boost/range/irange.hpp>
#include <boost/range/join.hpp>

#include "common/algorithm/constexpr_clamp.h"
#include "common/algorithm/constexpr_floor.h"
#include "common/cpp14/cbegin_cend.h"
#include "common/debug_numeric_cast.h"
#include "common/size_t_literal.h"
#include "exception/invalid_argument_exception.h"
#include "exception/not_implemented_exception.h"
#include "structure/view_cache/index/detail/vcie_match_criteria.h"
#include "structure/view_cache/index/view_cache_index_entry.h"

#include <cstddef>
#include <vector>

using namespace cath::common::literals;

namespace cath {
	namespace index {
		namespace detail {
			namespace detail {
				namespace detail {

					/// \brief TODOCUMENT
					template <typename T>
					constexpr inline int cell_index_of_value_in_current(const T   &arg_cell_width,   ///< TODOCUMENT
					                                                    const int &arg_start_offset, ///< TODOCUMENT
					                                                    const T   &arg_value         ///< TODOCUMENT
					                                                    ) {
						/// \todo Consider changing this to a numeric_cast (or debug_numeric_cast) if/when numeric_cast becomes constexpr
						///
						/// \todo Handle: libstdc++'s floor is constexpr but the standard doesn't define it as constexpr and
						///       libc++'s isn't. Find some constexpr workaround here. Is it adequate to just static_cast the fraction
						///       to int?
						return static_cast<int>( common::constexpr_floor( arg_value / arg_cell_width ) ) - arg_start_offset;
					}

					/// \brief TODOCUMENT
					template <typename T>
					constexpr inline size_t clamped_cell_index_of_value_in_current(const T      &arg_cell_width,   ///< TODOCUMENT
					                                                               const int    &arg_start_offset, ///< TODOCUMENT
					                                                               const size_t &arg_num_cells,    ///< TODOCUMENT
					                                                               const T      &arg_value         ///< TODOCUMENT
					                                                               ) {
						/// \todo Consider changing this to a numeric_cast (or debug_numeric_cast) if/when numeric_cast becomes constexpr
						return static_cast<size_t>( common::constexpr_clamp(
							cell_index_of_value_in_current( arg_cell_width, arg_start_offset, arg_value ),
							0,
							static_cast<int>( arg_num_cells )
						) );
					}

					/// \brief TODOCUMENT
					template <typename T>
					constexpr inline T min_value_in_cell_of_index_in_current(const T      &arg_cell_width,   ///< TODOCUMENT
					                                                         const int    &arg_start_offset, ///< TODOCUMENT
					                                                         const size_t &arg_cell_index    ///< TODOCUMENT
					                                                         ) {
						return cath::debug_numeric_cast<multiplier_type>(
							arg_start_offset + cath::debug_numeric_cast<int>( arg_cell_index )
						) * arg_cell_width;
					}

					/// \brief TODOCUMENT
					template <typename T>
					constexpr inline size_size_size_size_tpl search_cell_ranges(const T      &arg_cell_width,   ///< TODOCUMENT
					                                                            const int    &arg_start_offset, ///< TODOCUMENT
					                                                            const size_t &arg_num_cells,    ///< TODOCUMENT
					                                                            const T      &arg_search_begin, ///< TODOCUMENT
					                                                            const T      &arg_search_end    ///< TODOCUMENT
					                                                            ) {
						// The expressions for calculating the end values adds arg_cell_width to arg_search_end
						// so that the generate a one-past-the-end (in the standard C++/STL manner)
						return std::make_tuple(
							                                         clamped_cell_index_of_value_in_current( arg_cell_width, arg_start_offset, arg_num_cells, arg_search_begin                ),
							( arg_search_begin <  arg_search_end ) ? clamped_cell_index_of_value_in_current( arg_cell_width, arg_start_offset, arg_num_cells, arg_search_end + arg_cell_width ) : arg_num_cells,
							0_z,
							( arg_search_begin >= arg_search_end ) ? clamped_cell_index_of_value_in_current( arg_cell_width, arg_start_offset, arg_num_cells, arg_search_end + arg_cell_width ) : 0_z
						);
					}
				}

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
					view_cache_index_dim_linear(const value_type &);

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
				size_size_size_size_tpl view_cache_index_dim_linear<T>::search_cell_ranges(const size_t     &arg_num_cells,     ///< TODOCUMENT
				                                                                           const value_type &arg_search_begin, ///< TODOCUMENT
				                                                                           const value_type &arg_search_end    ///< TODOCUMENT
				                                                                           ) const {
					return detail::search_cell_ranges(
						cell_width,
						start_offset,
						arg_num_cells,
						arg_search_begin,
						arg_search_end
					);
				}

				/// \brief TODOCUMENT
				template <typename T>
				template <typename CELLS>
				bool view_cache_index_dim_linear<T>::has_cell_at_value(const CELLS      &arg_cells, ///< TODOCUMENT
				                                                       const value_type &arg_value  ///< TODOCUMENT
				                                                       ) const {
					const int cell_index_in_current = cell_index_of_value_in_current( arg_value );
					return ( cell_index_in_current > 0 && static_cast<size_t>( cell_index_in_current ) < arg_cells.size() );
				}

				/// \brief TODOCUMENT
				template <typename T>
				template <typename CELLS>
				typename CELLS::value_type & view_cache_index_dim_linear<T>::cell_at_value(CELLS                            &arg_cells,        ///< TODOCUMENT
				                                                                           const typename CELLS::value_type &arg_default_cell, ///< TODOCUMENT
				                                                                           const value_type                 &arg_value         ///< TODOCUMENT
				                                                                           ) {
					const int cell_index_in_current = cell_index_of_value_in_current( arg_value );
					if ( arg_cells.empty() ) {
						arg_cells.assign( 1, arg_default_cell );
						start_offset = cell_index_in_current;
						return arg_cells.front();
					}
					else if ( cell_index_in_current < 0 ) {
						const size_t num_to_prepend = cath::debug_numeric_cast<size_t>( -cell_index_in_current );
						const int new_start_offest = get_start_offset() + cell_index_in_current;
						// Come GCC v4.9 above, replace this begin() with cbegin()
						arg_cells.insert( std::begin( arg_cells ), num_to_prepend, arg_default_cell );
						start_offset = new_start_offest;
						return arg_cells.front();
					}
					else if ( static_cast<size_t>( cell_index_in_current ) >= arg_cells.size() ) {
						const size_t new_size    = 1 + cath::debug_numeric_cast<size_t>( cell_index_in_current );
						arg_cells.resize( new_size, arg_default_cell );
						return arg_cells.back();
					}
					else {
						return arg_cells[ cath::debug_numeric_cast<size_t>( cell_index_in_current ) ];
					}
				}

				/// \brief TODOCUMENT
				template <typename T>
				template <typename CELLS>
				const typename CELLS::value_type & view_cache_index_dim_linear<T>::cell_at_value(const CELLS      &arg_cells, ///< TODOCUMENT
				                                                                                 const value_type &arg_value  ///< TODOCUMENT
				                                                                                 ) const {
	#ifndef NDEBUG
					if ( arg_cells.empty() ) {
						BOOST_THROW_EXCEPTION(cath::common::invalid_argument_exception("Cannot get entry at_value() with no populated cells"));
					}
					if ( ! has_cell_at_value( arg_cells, arg_value ) ) {
						BOOST_THROW_EXCEPTION(cath::common::invalid_argument_exception("Value is not found in any of the cells"));
					}
	#endif
					return arg_cells[ cell_index_of_value_in_current( arg_value ) ];
				}

				/// \brief TODOCUMENT
				template <typename T>
				view_cache_index_dim_linear<T>::view_cache_index_dim_linear(const value_type &arg_cell_width ///< TODOCUMENT
				                                                            ) : cell_width   ( arg_cell_width ),
				                                                                start_offset ( 0              ) {
					T().check_cell_width( get_cell_width() );
				}

				/// \brief The cell index of the specified value in the current state of the index
				template <typename T>
				int view_cache_index_dim_linear<T>::cell_index_of_value_in_current(const value_type &arg_value ///< TODOCUMENT
				                                                                   ) const {
					return detail::cell_index_of_value_in_current( cell_width, start_offset, arg_value );
				}

				/// \brief TODOCUMENT
				template <typename T>
				typename view_cache_index_dim_linear<T>::value_type view_cache_index_dim_linear<T>::min_value_in_cell_of_index_in_current(const size_t &arg_cell_index ///< TODOCUMENT
				                                                                                                                          ) const {
					return detail::min_value_in_cell_of_index_in_current( cell_width, start_offset, arg_cell_index );
				}

				/// \brief TODOCUMENT
				template <typename T>
				template <typename CELLS, typename DEFAULTS>
				void view_cache_index_dim_linear<T>::store(const view_cache_index_entry &arg_entry,   ///< TODOCUMENT
				                                           CELLS                        &arg_cells,   ///< TODOCUMENT
				                                           const DEFAULTS               &arg_defaults ///< TODOCUMENT
				                                           ) {
					// Grab the value in question and TODOCUMENT
					const value_type &value = T().get_index_value( arg_entry );
					cell_at_value( arg_cells, arg_defaults.get_head(), value ).store(
						arg_entry,
						arg_defaults.get_tail()
					);
				}

				/// \brief TODOCUMENT
				template <typename T>
				template <typename CELLS, typename ACTN>
				void view_cache_index_dim_linear<T>::perform_action_on_matches(const view_cache_index_entry &arg_entry,    ///< TODOCUMENT
				                                                               const CELLS                  &arg_cells,    ///< TODOCUMENT
				                                                               const vcie_match_criteria    &arg_criteria, ///< TODOCUMENT
				                                                               ACTN                         &arg_action    ///< TODOCUMENT
				                                                               ) const {
					const value_type &value         = T().get_index_value  ( arg_entry    );
					const value_type &search_radius = T().get_search_radius( arg_criteria );
					const size_t      num_cells     = arg_cells.size();

					if ( ! arg_cells.empty() ) {
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
							arg_cells[ search_cell_ctr ].perform_action_on_matches( arg_entry, arg_criteria, arg_action );
						}
					}
				}

				/// \brief TODOCUMENT
				template <typename T>
				template <typename CELLS, typename ACTN>
				inline void view_cache_index_dim_linear<T>::perform_action_on_all_match_at_nodes(const CELLS                          &arg_query_cells,     ///< TODOCUMENT
				                                                                                 const view_cache_index_dim_linear<T> &arg_match_dimension, ///< TODOCUMENT
				                                                                                 const CELLS                          &arg_match_cells,     ///< TODOCUMENT
				                                                                                 const vcie_match_criteria            &arg_criteria,        ///< TODOCUMENT
				                                                                                 ACTN                                 &arg_action           ///< TODOCUMENT
				                                                                                 ) const {
					const value_type &search_radius   = T().get_search_radius( arg_criteria );


					if ( ! arg_query_cells.empty() && ! arg_match_cells.empty() ) {
						const size_t &num_query_cells = arg_query_cells.size();

						for (const size_t &cell_ctr : boost::irange( 0_z, num_query_cells ) ) {
							const auto &query_cell = arg_query_cells[ cell_ctr ];
							if ( ! query_cell.empty() ) {

								const auto prepped_begin_and_end = T().prepare_search_begin_and_end(
									min_value_in_cell_of_index_in_current( cell_ctr     ) - search_radius,
									min_value_in_cell_of_index_in_current( cell_ctr + 1 ) + search_radius
								);

								const auto search_ranges = arg_match_dimension.search_cell_ranges(
									arg_match_cells.size(),
									prepped_begin_and_end.first,
									prepped_begin_and_end.second
								);
								const auto joined_iranges = boost::range::join(
									boost::irange( std::get<0>( search_ranges ), std::get<1>( search_ranges ) ),
									boost::irange( std::get<2>( search_ranges ), std::get<3>( search_ranges ) )
								);
								for (const size_t &match_cell_ctr : joined_iranges) {
									const auto &match_cell = arg_match_cells[ match_cell_ctr ];
									if ( ! match_cell.empty() ) {
										// std::cerr << std::string( 10 - CELLS::value_type::num_dims_remaining, ' ' )
										//           << "In dimension "                          << T().get_name()
										//           << ", scanning for matches to query cell [" << min_value_in_cell_of_index_in_current( cell_ctr     )
										//           << ", "                                     << min_value_in_cell_of_index_in_current( cell_ctr + 1 )
										//           << "].size_"                                << query_cell.get_num_cells()
										//           << ", trying match cell ["                  << arg_match_dimension.min_value_in_cell_of_index_in_current( match_cell_ctr     )
										//           << ", "                                     << arg_match_dimension.min_value_in_cell_of_index_in_current( match_cell_ctr + 1 )
										//           << "].size_"                                << match_cell.get_num_cells()
										//           << std::endl;

										query_cell.perform_action_on_all_match_at_nodes( match_cell, arg_criteria, arg_action );
									}
								}
							}
						}
					}
				}

			}
		}
	}
}

#endif
