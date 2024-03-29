/// \file
/// \brief The multi_struc_res_rep_pair_list class header

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

#ifndef CATH_TOOLS_SOURCE_CT_UNI_CATH_SCAN_DETAIL_RES_PAIR_MULTI_STRUC_RES_REP_PAIR_LIST_HPP
#define CATH_TOOLS_SOURCE_CT_UNI_CATH_SCAN_DETAIL_RES_PAIR_MULTI_STRUC_RES_REP_PAIR_LIST_HPP

#include "cath/scan/detail/res_pair/multi_struc_res_rep_pair.hpp"
#include "cath/scan/detail/scan_type_aliases.hpp"

#include <cstddef>
#include <type_traits>

// clang-format off
namespace cath::scan::detail { class multi_struc_res_rep_pair; }
// clang-format on

namespace cath::scan::detail {

	/// \brief Store an ordered list of multi_struc_res_rep_pair objects
	///
	/// This is useful for implementing a cell of reps for all-vs-all scanning
	class multi_struc_res_rep_pair_list final {
	private:
		/// \brief The multi_struc_res_rep_pairs, stored in a vector
		multi_struc_res_rep_pair_vec multi_struc_res_rep_pairs;

	public:
		/// \brief const_iterator type alias at part of making this a range over the multi_struc_res_rep_pairs
		using const_iterator = multi_struc_res_rep_pair_vec_citr;
		using iterator       = multi_struc_res_rep_pair_vec_citr;

		using value_type = typename multi_struc_res_rep_pair_vec::value_type;

		/// \brief Default ctor
		multi_struc_res_rep_pair_list() = default;
		explicit multi_struc_res_rep_pair_list(multi_struc_res_rep_pair_vec);

		[[nodiscard]] bool               empty() const;
		[[nodiscard]] size_t             size() const;
		const multi_struc_res_rep_pair & operator[](const size_t &) const;

		template <class... Ts>
		void emplace_back(Ts&& ...);
		void push_back(const multi_struc_res_rep_pair &);

		[[nodiscard]] const_iterator begin() const;
		[[nodiscard]] const_iterator end() const;
	};

	/// \brief Ctor from a vector of multi_struc_res_rep_pair objects
	inline multi_struc_res_rep_pair_list::multi_struc_res_rep_pair_list(multi_struc_res_rep_pair_vec prm_multi_struc_res_rep_pairs ///< The vector of multi_struc_res_rep_pairs from which to construct the multi_struc_res_rep_pair_list
	                                                                    ) : multi_struc_res_rep_pairs { std::move( prm_multi_struc_res_rep_pairs ) } {
	}

	/// \brief Return whether this multi_struc_res_rep_pair_list is empty
	inline bool multi_struc_res_rep_pair_list::empty() const {
		return multi_struc_res_rep_pairs.empty();
	}

	/// \brief Return the number of multi_struc_res_rep_pair entries
	inline size_t multi_struc_res_rep_pair_list::size() const {
		return multi_struc_res_rep_pairs.size();
	}

	/// \brief Standard subscript operator
	inline const multi_struc_res_rep_pair & multi_struc_res_rep_pair_list::operator[](const size_t &prm_index ///< The index of the entry to be queried
	                                                                                  ) const {
		return multi_struc_res_rep_pairs[ prm_index ];
	}

	/// \brief TODOCUMENT
	template <class... Ts>
	void multi_struc_res_rep_pair_list::emplace_back(Ts&& ... prm_values ///< TODOCUMENT
	                                                 ) {
		multi_struc_res_rep_pairs.emplace_back( std::forward<Ts>( prm_values ) ... );
	}

	/// \brief TODOCUMENT
	inline void multi_struc_res_rep_pair_list::push_back(const multi_struc_res_rep_pair &prm_res_pair ///< TODOCUMENT
	                                                     ) {
		multi_struc_res_rep_pairs.push_back( prm_res_pair );
	}

	/// \brief Standard const begin() method
	inline auto multi_struc_res_rep_pair_list::begin() const -> const_iterator {
		return ::std::cbegin( multi_struc_res_rep_pairs );
	}

	/// \brief Standard const end() method
	inline auto multi_struc_res_rep_pair_list::end() const -> const_iterator {
		return ::std::cend  ( multi_struc_res_rep_pairs );
	}

} // namespace cath::scan::detail

#endif // CATH_TOOLS_SOURCE_CT_UNI_CATH_SCAN_DETAIL_RES_PAIR_MULTI_STRUC_RES_REP_PAIR_LIST_HPP
