/// \file
/// \brief The view_cache_index_dim_orient class header

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

#ifndef CATH_TOOLS_SOURCE_CT_UNI_CATH_STRUCTURE_VIEW_CACHE_INDEX_DETAIL_DIMS_VIEW_CACHE_INDEX_DIM_ORIENT_HPP
#define CATH_TOOLS_SOURCE_CT_UNI_CATH_STRUCTURE_VIEW_CACHE_INDEX_DETAIL_DIMS_VIEW_CACHE_INDEX_DIM_ORIENT_HPP

// ********************************************************
// ********************************************************
// ****                                                ****
// **** !! AT PRESENT, THIS HEADER IS COMMENTED OUT !! ****
// **** !! BECAUSE IT DOES NOT COMPILE              !! ****
// ****                                                ****
// ********************************************************
// ********************************************************

/*

#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/range/begin.hpp>

#include "cath/common/algorithm/sort_uniq_copy.hpp"
#include "cath/common/debug_numeric_cast.hpp"
#include "cath/common/exception/invalid_argument_exception.hpp"
#include "cath/common/exception/not_implemented_exception.hpp"
#include "cath/common/type_aliases.hpp"
#include "cath/structure/view_cache/index/detail/vcie_match_criteria.hpp"
#include "cath/structure/view_cache/index/detail/view_cache_index_type_aliases.hpp"
#include "cath/structure/view_cache/index/view_cache_index_entry.hpp"

#include <cstddef>
#include <vector>

namespace cath::index::detail {

	/// \brief TODOCUMENT
	///
	/// These classes (view_cache_index_dim_orient<>, view_cache_index_dim_dirn) don't hold the next levels down to the data being indexed.
	/// They simply specify how to manage that data, being passed a suitable vector as necessary.
	class view_cache_index_dim_orient final {
	private:
		/// \brief TODOCUMENT
		angle_type search_radius;

		/// \brief TODOCUMENT
		frame_quat_rot_vec cell_width;

		/// \brief TODOCUMENT
		size_vec_vec neighbours;

		static size_vec_vec calc_neighbours(const frame_quat_rot_vec &,
		                                    const angle_type &);

		// const value_type & get_cell_width() const;
		// const int & get_start_offset() const;

		size_t cell_index_of_value_in_current(const frame_quat_rot &) const;

		// template <typename CELLS>
		// bool has_cell_at_value(const CELLS &,
		//                        const value_type &) const;

		template <typename CELLS>
		typename CELLS::value_type & cell_at_value(CELLS &,
		                                           const typename CELLS::value_type &,
		                                           const frame_quat_rot &);

		template <typename CELLS>
		const typename CELLS::value_type & cell_at_value(const CELLS &,
		                                                 const frame_quat_rot &) const;

	public:
		view_cache_index_dim_orient(const detail::vcie_match_criteria &,
		                            const frame_quat_rot_vec &);

		template <typename CELLS, typename DEFAULTS>
		void store(const view_cache_index_entry &,
		           CELLS &,
		           const DEFAULTS &);

		template <typename CELLS, typename ACTN>
		void perform_action_on_matches(const view_cache_index_entry &,
		                               const CELLS &,
		                               const detail::vcie_match_criteria &,
		                               ACTN &) const;
	};

	/// \brief TODOCUMENT
	///
	/// Note: this makes a few mathematical assumptions that I haven't completely checked
	///       so it's important to test that this calculates correct results...
	///
	/// \todo Test this by scanning the all-against-all anchor pairs and then for each, check:
	///        * the half-way quaternion is the same angle to both of the pair
	///        * the pair's in neighbours iff the half-way quaternion is within prm_search_radius of both
	///
	/// \todo Once the above tests are in place, try changing the code to calculate the
	///       distance_1_of_angle(2.0 * prm_search_radius) at the start and then compare
	///       the distance_1_of_quat_rot( anchor, neighbour ) to that value
	size_vec_vec view_cache_index_dim_orient::calc_neighbours(const frame_quat_rot_vec &prm_anchor_quat_rots, ///< TODOCUMENT
	                                                          const angle_type         &prm_search_radius     ///< TODOCUMENT
	                                                          ) {
		const size_t num_anchor_quat_rots = prm_anchor_quat_rots.size();

		size_vec_vec new_neighbours;
		new_neighbours.reserve( num_anchor_quat_rots );

		for (const frame_quat_rot &anchor : prm_anchor_quat_rots) {
			size_vec local_neighbours;
			for (const size_t &neighbour_ctr : indices( num_anchor_quat_rots ) ) {
				const frame_quat_rot &neighbour = prm_anchor_quat_rots[ neighbour_ctr ];
				if ( angle_between_quat_rots( anchor, neighbour ) / 2.0 < prm_search_radius ) {
					local_neighbours.push_back( neighbour_ctr );
				}
			}
			new_neighbours.push_back( common::sort_copy( local_neighbours ) );
		}
		return new_neighbours;
	}

	/// \brief Getter for the cell width
	const typename view_cache_index_dim_orient::value_type & view_cache_index_dim_orient::get_cell_width() const {
		return cell_width;
	}

	/// \brief Getter for the start offset
	const int & view_cache_index_dim_orient::get_start_offset() const {
		return start_offset;
	}

	/// \brief The cell index of the specified value in the current state of the index
	int view_cache_index_dim_orient::cell_index_of_value_in_current(const value_type &prm_value ///< TODOCUMENT
	                                                                ) const {
		return cath::debug_numeric_cast<int>( floor( prm_value / cell_width) ) - start_offset;
	}

	/// \brief TODOCUMENT
	template <typename CELLS>
	bool view_cache_index_dim_orient::has_cell_at_value(const CELLS      &prm_cells, ///< TODOCUMENT
	                                                    const value_type &prm_value  ///< TODOCUMENT
	                                                    ) const {
		const int cell_index_in_current = cell_index_of_value_in_current( prm_value );
		return ( cell_index_in_current > 0 && static_cast<size_t>( cell_index_in_current ) < prm_cells.size() );
	}

	/// \brief TODOCUMENT
	template <typename CELLS>
	typename CELLS::value_type & view_cache_index_dim_orient::cell_at_value(CELLS                            &prm_cells,        ///< TODOCUMENT
	                                                                        const typename CELLS::value_type &prm_default_cell, ///< TODOCUMENT
	                                                                        const value_type                 &prm_value         ///< TODOCUMENT
	                                                                        ) {
		const int cell_index_in_current = cell_index_of_value_in_current( prm_value );
		if ( prm_cells.empty() ) {
			prm_cells.assign( 1, prm_default_cell );
			start_offset = cell_index_in_current;
			return prm_cells.front();
		}
		else if ( cell_index_in_current < 0 ) {
			const size_t num_to_prepend = cath::debug_numeric_cast<size_t>( -cell_index_in_current );
			const int new_start_offest = get_start_offset() + cell_index_in_current;
			prm_cells.insert( ::std::cbegin( prm_cells ), num_to_prepend, prm_default_cell );
			start_offset = new_start_offest;
			return prm_cells.front();
		}
		else if ( static_cast<size_t>( cell_index_in_current ) >= prm_cells.size() ) {
			const size_t new_size    = 1 + cath::debug_numeric_cast<size_t>( cell_index_in_current );
			prm_cells.resize( new_size, prm_default_cell );
			return prm_cells.back();
		}
		else {
			return prm_cells[ cath::debug_numeric_cast<size_t>( cell_index_in_current ) ];
		}
	}

	/// \brief TODOCUMENT
	template <typename CELLS>
	const typename CELLS::value_type & view_cache_index_dim_orient::cell_at_value(const CELLS      &prm_cells, ///< TODOCUMENT
	                                                                              const value_type &prm_value  ///< TODOCUMENT
	                                                                              ) const {
		if constexpr ( IS_IN_DEBUG_MODE ) {
			if ( prm_cells.empty() ) {
				BOOST_THROW_EXCEPTION(cath::common::invalid_argument_exception("Cannot get entry at_value() with no populated cells"));
			}
			if ( ! has_cell_at_value( prm_cells, prm_value ) ) {
				BOOST_THROW_EXCEPTION(cath::common::invalid_argument_exception("Value is not found in any of the cells"));
			}
		}
		return prm_cells[ cell_index_of_value_in_current( prm_value ) ];
	}

	/// \brief TODOCUMENT
	view_cache_index_dim_orient::view_cache_index_dim_orient(const value_type &prm_cell_width ///< TODOCUMENT
	                                                         ) : cell_width   ( prm_cell_width ),
	                                                             start_offset ( 0              ) {
		T().check_cell_width( get_cell_width() );
	}

	/// \brief TODOCUMENT
	template <typename CELLS, typename DEFAULTS>
	void view_cache_index_dim_orient::store(const view_cache_index_entry &prm_entry,   ///< TODOCUMENT
	                                        CELLS                        &prm_cells,   ///< TODOCUMENT
	                                        const DEFAULTS               &prm_defaults ///< TODOCUMENT
	                                        ) {
		// Grab the value in question and TODOCUMENT
		const value_type &value = T().get_index_value( prm_entry );
		cell_at_value( prm_cells, prm_defaults.get_head(), value ).store(
			prm_entry,
			prm_defaults.get_tail()
		);
	}

	/// \brief TODOCUMENT
	template <typename CELLS, typename ACTN>
	void view_cache_index_dim_orient::perform_action_on_matches(const view_cache_index_entry      &prm_entry,    ///< TODOCUMENT
	                                                            const CELLS                       &prm_cells,    ///< TODOCUMENT
	                                                            const detail::vcie_match_criteria &prm_criteria, ///< TODOCUMENT
	                                                            ACTN                              &prm_action    ///< TODOCUMENT
	                                                            ) const {
		const value_type &value         = T().get_index_value  ( prm_entry    );
		const value_type &search_radius = T().get_search_radius( prm_criteria );

		const int num_cells      = debug_numeric_cast<int>( prm_cells.size() - 1 );
		const int min_cell_index = std::max( 0,         cell_index_of_value_in_current( value - search_radius ) );
		const int max_cell_index = std::min( num_cells, cell_index_of_value_in_current( value + search_radius ) );

		for (const int &cell_index : irange( min_cell_index, max_cell_index + 1 ) ) {
			prm_cells[ cell_index ].perform_action_on_matches( prm_entry, prm_criteria, prm_action );
		}
	}

	/// \brief TODOCUMENT
	size_t index_of_closest_quat_rep(const spanning_quads &prm_spanning_quats, ///< TODOCUMENT
	                                 const frame_quat_rot &prm_quat_rot        ///< TODOCUMENT
	                                 ) {
		const auto closest_itr = boost::range::min_element(
			prm_spanning_quats,
			[] (const frame_quat_rot &x, const frame_quat_rot &y) {
				return (
					distance_1_between_quat_rots( prm_quat_rot, x )
					<
					distance_1_between_quat_rots( prm_quat_rot, y )
				);
			}
		);
		if ( closest_itr == ::std::cend( prm_spanning_quats ) ) {
			BOOST_THROW_EXCEPTION(out_of_range_exception(""));
		}
		return boost::numeric_cast<size_t>( std::distance(
			::std::cbegin( prm_spanning_quats ),
			closest_itr
		) );
	}

	/// \brief TODOCUMENT
	size_t index_of_closest_quat_rep(const spanning_quads &prm_spanning_quats, ///< TODOCUMENT
	                                 const frame_quat_rot &prm_quat_rot        ///< TODOCUMENT
	                                 ) {
		const auto closest_itr = boost::range::min_element(
			prm_spanning_quats,
			[] (const frame_quat_rot &x, const frame_quat_rot &y) {
				return (
					distance_1_between_quat_rots( prm_quat_rot, x )
					<
					distance_1_between_quat_rots( prm_quat_rot, y )
				);
			}
		);
		if ( closest_itr == ::std::cend( prm_spanning_quats ) ) {
			BOOST_THROW_EXCEPTION(out_of_range_exception(""));
		}
		return boost::numeric_cast<size_t>( std::distance(
			::std::cbegin( prm_spanning_quats ),
			closest_itr
		) );
	}

} // namespace cath::index::detail

*/

#endif // CATH_TOOLS_SOURCE_CT_UNI_CATH_STRUCTURE_VIEW_CACHE_INDEX_DETAIL_DIMS_VIEW_CACHE_INDEX_DIM_ORIENT_HPP
