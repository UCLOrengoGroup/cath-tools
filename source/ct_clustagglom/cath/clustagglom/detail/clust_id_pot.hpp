/// \file
/// \brief The clust_id_pot class header

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

#ifndef _CATH_TOOLS_SOURCE_CT_CLUSTAGGLOM_CATH_CLUSTAGGLOM_DETAIL_CLUST_ID_POT_HPP
#define _CATH_TOOLS_SOURCE_CT_CLUSTAGGLOM_CATH_CLUSTAGGLOM_DETAIL_CLUST_ID_POT_HPP

#include <optional>

#include <boost/range/adaptor/filtered.hpp>
#include <boost/range/algorithm/min_element.hpp>

#include "cath/clustagglom/clustagglom_type_aliases.hpp"
#include "cath/common/algorithm/copy_build.hpp"
#include "cath/common/boost_addenda/range/indices.hpp"
#include "cath/common/debug_numeric_cast.hpp"
#include "cath/common/size_t_literal.hpp"
#include "cath/common/type_aliases.hpp"

using namespace ::cath::common::literals;

namespace cath::clust::detail {

	/// \brief Implements a pot of indices, ascending from 0, that supports O(1) versions of operations:
	///         * add a new index that's one larger than any previous index
	///         * remove some specified index
	///         * query whether some index remains present
	///         * get the n-th index for some ordering, which is arbitrary
	///           but consistent across calls to const methods
	class clust_id_pot final {
	private:
		/// \brief The present indices, possibly out of order
		item_vec jumbled_values;

		/// \brief Correspond to the full list of all indices that have been in this pot and, in each case,
		///        use a size_opt to indicate whether that index is still present and, if so, the index
		///        of jumbled_values in which it's currently stored.
		size_opt_vec indices;

	public:
		clust_id_pot() = delete;
		explicit clust_id_pot(const size_t &);

		[[nodiscard]] const item_idx &get_jumbled_nth_index( const size_t & ) const;
		[[nodiscard]] bool            has_index( const item_idx & ) const;
		clust_id_pot & remove_index(const item_idx &);
		const item_idx & add_new_index();

		[[nodiscard]] item_idx get_min_value_excluding_spec( const item_idx & ) const;
	};

	/// \brief Ctor from the number of items
	inline clust_id_pot::clust_id_pot(const size_t &prm_num_indices ///< The number of items with which to initialise
	                                  ) : jumbled_values( common::copy_build<item_vec    >( common::indices( prm_num_indices ) ) ),
	                                      indices       ( common::copy_build<size_opt_vec>( common::indices( prm_num_indices ) ) ) {
	}

	/// \brief get the n-th index for some ordering, which is arbitrary
	///           but consistent across calls to const methods
	inline const item_idx & clust_id_pot::get_jumbled_nth_index(const size_t &prm_index ///< The index of the index to get under some arbitrary ordering
	                                                            ) const {
		return jumbled_values[ prm_index ];
	}

	/// \brief Return whether the specified index is still present in the pot
	inline bool clust_id_pot::has_index(const item_idx &prm_index ///< The index to query
	                                    ) const {
		return static_cast<bool>( indices[ prm_index ] );
	}

	/// \brief Remove the specified index
	inline clust_id_pot & clust_id_pot::remove_index(const item_idx &prm_index ///< The index to remove
	                                                 ) {
		const auto x_index = *indices[ prm_index ];
		if ( x_index + 1 < jumbled_values.size() ) {
			std::swap( jumbled_values[ x_index ], jumbled_values.back() );
			indices[ jumbled_values[ x_index ] ] = x_index;
		}
		indices[ jumbled_values.back() ] = ::std::nullopt;
		jumbled_values.pop_back();
		return *this;
	}

	/// \brief Add a new index that's one larger than any previous index and return it
	inline const item_idx & clust_id_pot::add_new_index() {
		jumbled_values.push_back( debug_numeric_cast<item_idx>( indices.size() ) );
		indices.push_back( debug_numeric_cast<item_idx>( jumbled_values.size() - 1_z ) );
		return jumbled_values.back();
	}

	/// \brief Get the minimum index, excluding the specified one
	inline item_idx clust_id_pot::get_min_value_excluding_spec(const item_idx &prm_item ///< The index to exclude
	                                                           ) const {
		return *::boost::range::min_element(
			jumbled_values
				| ::boost::adaptors::filtered(
					[ & ]( const item_idx &x ) { return ( x != prm_item ); }
				)
		);
	}

} // namespace cath::clust::detail

#endif // _CATH_TOOLS_SOURCE_CT_CLUSTAGGLOM_CATH_CLUSTAGGLOM_DETAIL_CLUST_ID_POT_HPP
