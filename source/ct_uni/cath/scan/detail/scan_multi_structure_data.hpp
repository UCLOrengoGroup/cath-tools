/// \file
/// \brief The scan_multi_structure_data class header

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

#ifndef _CATH_TOOLS_SOURCE_CT_UNI_CATH_SCAN_DETAIL_SCAN_MULTI_STRUCTURE_DATA_HPP
#define _CATH_TOOLS_SOURCE_CT_UNI_CATH_SCAN_DETAIL_SCAN_MULTI_STRUCTURE_DATA_HPP

#include "cath/common/boost_addenda/range/front.hpp"
#include "cath/scan/detail/res_pair/multi_struc_res_rep_pair_list.hpp"
#include "cath/scan/detail/scan_structure_data.hpp"

namespace cath::scan::detail {

	/// \brief For a list of structures, store the data needed for following up a multi_struc_res_rep_pair match
	///
	/// For more information on the information held, see scan_structure_data
	///
	/// This is roughly a wrapper for a vector of scan_structure_data objects.
	/// It won't do anything to identify and prevent duplicates.
	///
	/// It might be worth making this store a reference (or reference_wrapper) back to the original protein?
	class scan_multi_structure_data final {
	private:
		/// \brief The vector of scan_structure_data objects
		scan_structure_data_vec structures_data;

	public:
		/// \brief const_iterator type alias to structures_data's const_iterator type
		///        as part of making scan_multi_structure_data into a random access range
		using const_iterator = scan_structure_data_vec_citr;

		using iterator = scan_structure_data_vec_citr;

		const scan_structure_data & operator[](const index_type &) const;

		[[nodiscard]] bool   empty() const;
		[[nodiscard]] size_t size() const;

		[[nodiscard]] const_iterator begin() const;
		[[nodiscard]] const_iterator end() const;

		template <class... Ts>
		void emplace_back(Ts&& ...);
		void push_back(const scan_structure_data &);
	};

	/// \brief Subscript operator to get a const reference to the scan_structure_data at the specified index
	inline const scan_structure_data & scan_multi_structure_data::operator[](const index_type &prm_index ///< The index of the scan_structure_data to return
	                                                                         ) const {
		return structures_data[ prm_index ];
	}

	/// \brief Whether this is empty (ie doesn't contain any data on any structures)
	inline bool scan_multi_structure_data::empty() const {
		return structures_data.empty();
	}

	/// \brief The number of structures for which this contains data
	inline size_t scan_multi_structure_data::size() const {
		return structures_data.size();
	}

	/// \brief Standard begin() operator returning a const_iterator over the scan_structure_data objects
	///
	/// This is part of making scan_multi_structure_data into a random access range
	inline auto scan_multi_structure_data::begin() const -> const_iterator {
		return ::std::cbegin( structures_data );
	}

	/// \brief Standard end() operator returning a const_iterator over the scan_structure_data objects
	///
	/// This is part of making scan_multi_structure_data into a random access range
	inline auto scan_multi_structure_data::end() const -> const_iterator {
		return ::std::cend  ( structures_data );
	}

	/// \brief TODOCUMENT
	template <class... Ts>
	void scan_multi_structure_data::emplace_back(Ts&& ... prm_values ///< TODOCUMENT
	                                             ) {
		structures_data.emplace_back( std::forward<Ts>( prm_values ) ... );
	}


	/// \brief Add the specified scan_structure_data to the end of the list
	///
	/// \pre TODOCUMENT (roled_scan_strides must match)
	inline void scan_multi_structure_data::push_back(const scan_structure_data &prm_scan_structure_data ///< The scan_structure_data to append to the scan_multi_structure_data
	                                                 ) {
		if ( ! empty() ) {
			if ( common::front( *this ).get_roled_scan_stride() != prm_scan_structure_data.get_roled_scan_stride() ) {
				BOOST_THROW_EXCEPTION(common::invalid_argument_exception("Cannot add a scan_data_structure to an existing scan_multi_structure_data which already has scan_data_structures with a different roled_scan_stride"));
			}
		}
		structures_data.push_back( prm_scan_structure_data );
	}

	/// \brief Get the list of single_struc_res_pair neighbours associated with
	///        the specified rep specified multi_struc_res_rep_pair
	///
	/// \relates scan_multi_structure_data
	inline const single_struc_res_pair_list & get_neighbours_of_rep_pair(const scan_multi_structure_data &prm_multi_structure_data, ///< The scan_multi_structure_data to query
	                                                                     const multi_struc_res_rep_pair  &prm_res_pair              ///< The rep multi_struc_res_rep_pair for which the neighbours should be extracted
	                                                                     ) {
		return prm_multi_structure_data[ prm_res_pair.get_structure_index() ].get_res_pairs_of_rep_indices(
			prm_res_pair.get_from_res_rep_index(),
			prm_res_pair.get_to_res_rep_index()
		);
	}

	/// \brief Build a new scan_structure_data from a protein and roled_scan_stride and
	///        append it to the specified scan_multi_structure_data
	///
	/// The roled_scan_stride is used to specify whether a query/index scan_structure_data object should be appended.
	/// (They're different because the scan_stride can specify different from/to strides for each.)
	///
	/// \relates scan_multi_structure_data
	inline void add_structure_data(scan_multi_structure_data &prm_scan_multi_structure_data, ///< TODOCUMENT
	                               const protein             &prm_protein,                   ///< TODOCUMENT
	                               const roled_scan_stride   &prm_roled_scan_stride          ///< TODOCUMENT
	                               ) {
		prm_scan_multi_structure_data.emplace_back(
			prm_protein,
			prm_roled_scan_stride
		);
	}

	/// \brief Returns whether the specified criteria are met by the specified multi_struc_res_rep_pair
	///
	/// NOTE: As indicated by the name, this can be used to rule some multi_struc_res_rep_pairs out of
	///       consideration but the criteria assess the similarity of two res_pairs so it doesn't make sense
	///       to ask whether the criteria are met by a single res_pair.
	///
	/// Violation can occur if the from_index matches the to_index or their absolute difference is
	/// otherwise less than the criteria's minimum_index_distance
	///
	/// \relates quad_criteria
	///
	/// \relatesalso scan_multi_structure_data
	inline bool are_not_violated_by(const quad_criteria             &prm_criteria,           ///< The criteria to apply
	                                const multi_struc_res_rep_pair  &prm_res_pair,           ///< The res_pair to be tested for whether it has any chance of matching another
	                                const scan_multi_structure_data &prm_scan_structure_data ///< The corresponding scan_multi_structure_data containing the data for this res_pair's structure
	                                ) {
		return are_not_violated_by(
			prm_criteria,
			prm_res_pair,
			prm_scan_structure_data[ prm_res_pair.get_structure_index() ]
		);
	}


	/// \brief TODOCUMENT
	///
	/// \relates scan_multi_structure_data
	template <typename FN>
	inline void act_on_multi_matches(const multi_struc_res_rep_pair_list &prm_list_a,               ///< TODOCUMENT
	                                 const multi_struc_res_rep_pair_list &prm_list_b,               ///< TODOCUMENT
	                                 const scan_multi_structure_data     &prm_query_structures_data, ///< TODOCUMENT
	                                 const scan_multi_structure_data     &prm_index_structures_data, ///< TODOCUMENT
	                                 const quad_criteria                 &prm_criteria,             ///< TODOCUMENT
//	                                 const scan_stride                   &prm_stride,               ///< TODOCUMENT
	                                 FN                                  &prm_fn                    ///< TODOCUMENT
	                                 ) {
//		const auto query_from_strider = prm_stride.get_query_from_strider();
//		const auto query_to_strider   = prm_stride.get_query_to_strider();
//		const auto index_from_strider = prm_stride.get_index_from_strider();
//		const auto index_to_strider   = prm_stride.get_index_to_strider();

		for (const multi_struc_res_rep_pair &res_pair_a : prm_list_a) {
//			const auto query_from_idx = get_index_of_rep_index( query_from_strider, res_pair_a.get_from_res_rep_index() );
//			const auto query_to_idx   = get_index_of_rep_index( query_to_strider,   res_pair_a.get_to_res_rep_index()   );
//			std::cerr << "Query[ "
//			          << query_from_idx
//			          << "(rep:"
//			          << res_pair_a.get_from_res_rep_index()
//			          << ")"
//			          << " - "
//			          << query_to_idx
//			          << "(rep:"
//			          << res_pair_a.get_to_res_rep_index()
//			          << ")"
//			          << " in "
//			          << res_pair_a.get_structure_index()
//			          << " ]; - "
//			          << res_pair_a
//			          << "\n";
			if ( are_not_violated_by( prm_criteria, res_pair_a ) ) {
				for (const multi_struc_res_rep_pair &res_pair_b : prm_list_b) {
//					const auto index_from_idx = get_index_of_rep_index( index_from_strider, res_pair_b.get_from_res_rep_index() );
//					const auto index_to_idx   = get_index_of_rep_index( index_to_strider,   res_pair_b.get_to_res_rep_index()   );
					if ( are_not_violated_by( prm_criteria, res_pair_a ) && are_met_by( prm_criteria, res_pair_a, res_pair_b ) ) {
//						std::cerr << "Investigating hit: query[ "
//						          << query_from_idx << " - " << query_to_idx
//						          << " in "
//						          << res_pair_a.get_structure_index()
//						          << " ]; index [ "
//						          << index_from_idx << " - " << index_to_idx
//						          << " in "
//						          << res_pair_b.get_structure_index()
//						          << " ]\n";
						act_on_single_matches(
							get_neighbours_of_rep_pair( prm_query_structures_data, res_pair_a ),
							get_neighbours_of_rep_pair( prm_index_structures_data, res_pair_b ),
							res_pair_a.get_structure_index(),
							res_pair_b.get_structure_index(),
							prm_criteria,
							prm_fn
						);
					}
//					else {
//						if ( are_not_violated_by( prm_criteria, res_pair_a ) ) {
//							std::cerr << "Skipping index " << index_from_idx << " - " << index_to_idx << "\n";
//						}
//						else if ( query_from_idx == index_from_idx && query_to_idx == index_to_idx) {
//							std::cerr << "**** SKIPPED A SELF-QUAD\n";
//						}
//					}
				}
			}
//			else {
//				std::cerr << "Skipping query " << query_from_idx << " - " << query_to_idx << "\n";
//			}
		}
	}

} // namespace cath::scan::detail

#endif // _CATH_TOOLS_SOURCE_CT_UNI_CATH_SCAN_DETAIL_SCAN_MULTI_STRUCTURE_DATA_HPP
