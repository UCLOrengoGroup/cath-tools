/// \file
/// \brief The single_struc_res_pair_list class header

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

#ifndef CATH_TOOLS_SOURCE_CT_UNI_CATH_SCAN_DETAIL_RES_PAIR_SINGLE_STRUC_RES_PAIR_LIST_HPP
#define CATH_TOOLS_SOURCE_CT_UNI_CATH_SCAN_DETAIL_RES_PAIR_SINGLE_STRUC_RES_PAIR_LIST_HPP

#include <boost/range/algorithm_ext/for_each.hpp>

#include "cath/scan/detail/quad_criteria_are_met_by.hpp"
#include "cath/scan/detail/res_pair/single_struc_res_pair.hpp"

namespace cath::scan::detail {

	/// \brief Store an ordered list of single_struc_res_pair objects
	///
	/// This is useful for implementing a cell of neighbours around a rep
	/// that can be compared one-to-one with the equivalent neighbours of another rep
	class single_struc_res_pair_list final {
	private:
		/// \brief The single_struc_res_rep_pairs, stored in a vector
		single_struc_res_pair_vec single_struc_res_pairs;

	public:
		/// \brief const_iterator type alias at part of making this a range over the single_struc_res_pairs
		using const_iterator = single_struc_res_pair_vec_citr;

		/// \brief Default ctor
		single_struc_res_pair_list() = default;
		explicit single_struc_res_pair_list(single_struc_res_pair_vec);

		void reserve(const size_t &);
		[[nodiscard]] bool   empty() const;
		[[nodiscard]] size_t size() const;

		const single_struc_res_pair & operator[](const size_t &) const;

		template <class... Ts>
		void emplace_back(Ts&& ...);
		void push_back(const single_struc_res_pair &);

		[[nodiscard]] const_iterator begin() const;
		[[nodiscard]] const_iterator end() const;
	};

	/// \brief Ctor from a vector of single_struc_res_pair objects
	inline single_struc_res_pair_list::single_struc_res_pair_list(single_struc_res_pair_vec prm_single_struc_res_pairs ///< TODOCUMENT
	                                                              ) : single_struc_res_pairs { std::move( prm_single_struc_res_pairs ) } {
	}

	/// \brief Return whether this single_struc_res_pair_list is empty
	inline void single_struc_res_pair_list::reserve(const size_t &prm_size ///< TODOCUMENT
	                                                ) {
		single_struc_res_pairs.reserve( prm_size );
	}

	/// \brief Return whether this single_struc_res_pair_list is empty
	inline bool single_struc_res_pair_list::empty() const {
		return single_struc_res_pairs.empty();
	}

	/// \brief Return the number of single_struc_res_pair entries
	inline size_t single_struc_res_pair_list::size() const {
		return single_struc_res_pairs.size();
	}

	/// \brief Standard subscript operator
	inline const single_struc_res_pair & single_struc_res_pair_list::operator[](const size_t &prm_index ///< The index of the entry to be queried
	                                                                            ) const {
		return single_struc_res_pairs[ prm_index ];
	}

	/// \brief TODOCUMENT
	template <class... Ts>
	inline void single_struc_res_pair_list::emplace_back(Ts&& ...prm_values ///< TODOCUMENT
	                                                     ) {
		single_struc_res_pairs.emplace_back( std::forward<Ts>( prm_values ) ... );
	}

	/// \brief TODOCUMENT
	inline void single_struc_res_pair_list::push_back(const single_struc_res_pair &prm_res_pair ///< TODOCUMENT
	                                                  ) {
		single_struc_res_pairs.push_back( prm_res_pair );
	}

	/// \brief Standard const begin() method
	inline auto single_struc_res_pair_list::begin() const -> const_iterator {
		return ::std::cbegin( single_struc_res_pairs );
	}

	/// \brief Standard const end() method
	inline auto single_struc_res_pair_list::end() const -> const_iterator {
		return ::std::cend  ( single_struc_res_pairs );
	}

	/// \brief Perform action on the corresponding pairwise entries in two single_struc_res_pair_lists that meet prm_criteria
	///
	/// \pre prm_list_a and prm_list_b should be of matching lengths (checked by an assert statement in debug build)
	///
	/// \pre All single_struc_res_pairs in prm_list_a and prm_list_b should pass are_not_violated_by( prm_criteria, x )
	///      for the specified quad_criteria (because this isn't checked)
	///
	/// \pre (implicit) prm_list_a and prm_list_b should be sorted to have equivalent entries in equivalent positions
	///
	/// Note that single_struc_res_pair_lists may contain dummy entries so this checks each entry is not
	/// a dummy before proceeding further
	///
	/// \relates single_struc_res_pair_list
	template <typename FN>
	inline void act_on_single_matches(const single_struc_res_pair_list &prm_list_a,      ///< TODOCUMENT
	                                  const single_struc_res_pair_list &prm_list_b,      ///< TODOCUMENT
	                                  const index_type                 &prm_structure_a, ///< TODOCUMENT
	                                  const index_type                 &prm_structure_b, ///< TODOCUMENT
	                                  const quad_criteria              &prm_criteria,    ///< TODOCUMENT
	                                  FN                               &prm_function     ///< TODOCUMENT
	                                  ) {
		assert( prm_list_a.size() == prm_list_b.size() );
		boost::range::for_each(
			prm_list_a,
			prm_list_b,
			[&] (const single_struc_res_pair &x, const single_struc_res_pair &y) {
				if ( ! is_dummy( x ) && ! is_dummy( y ) ) {
					/// \todo Consider moving these are_not_violated_by(const single_struc_res_pair &) tests
					///       into the building so that single_struc_res_pair entries that fail are never built int
					///       the scan_structure_data. This has the disadvantage of preventing
					if ( are_not_violated_by( prm_criteria, x ) && are_not_violated_by( prm_criteria, y ) ) {
						if ( are_met_by( prm_criteria, x, y ) ) {
//					std::cerr << "**** MATCHED : " << x
//					          << "(dummy:"         << std::right << std::setw( 5 ) << std::boolalpha << is_dummy( x ) << std::noboolalpha
//					          << ") vs "           << y
//					          << "(dummy:"         << std::right << std::setw( 5 ) << std::boolalpha << is_dummy( y ) << std::noboolalpha
//					          << ")\n";
							prm_function( x, y, prm_structure_a, prm_structure_b );
						}
					}
				}
//				else {
//					std::cerr << "**** SKIPPED : " << x
//					          << "(dummy:"         << std::right << std::setw( 5 ) << std::boolalpha << is_dummy( x ) << std::noboolalpha
//					          << ") vs "           << y
//					          << "(dummy:"         << std::right << std::setw( 5 ) << std::boolalpha << is_dummy( y ) << std::noboolalpha
//					          << ")\n";
//				}
			}
		);
	}

} // namespace cath::scan::detail

#endif // CATH_TOOLS_SOURCE_CT_UNI_CATH_SCAN_DETAIL_RES_PAIR_SINGLE_STRUC_RES_PAIR_LIST_HPP
