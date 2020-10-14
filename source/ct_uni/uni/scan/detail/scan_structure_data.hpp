/// \file
/// \brief The scan_structure_data class header

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

#ifndef _CATH_TOOLS_SOURCE_UNI_SCAN_DETAIL_SCAN_STRUCTURE_DATA_H
#define _CATH_TOOLS_SOURCE_UNI_SCAN_DETAIL_SCAN_STRUCTURE_DATA_H

#include <boost/throw_exception.hpp>
#include <boost/units/quantity.hpp>
#include <boost/units/systems/information/byte.hpp>

#include "common/boost_addenda/range/front.hpp"
#include "common/debug_numeric_cast.hpp"
#include "common/exception/out_of_range_exception.hpp"
#include "scan/detail/detail/scan_structure_data_helper.hpp"
#include "scan/detail/quad_criteria_are_met_by.hpp"
#include "scan/detail/res_pair/multi_struc_res_rep_pair.hpp"
#include "scan/detail/res_pair/single_struc_res_pair.hpp"
#include "scan/detail/res_pair/single_struc_res_pair_list.hpp"
#include "scan/detail/scan_type_aliases.hpp"
#include "scan/detail/stride/rep_strider.hpp"
#include "scan/detail/stride/roled_scan_stride.hpp"
#include "structure/protein/protein.hpp"

#include <cstddef>

namespace cath { namespace scan { class quad_criteria; } }

namespace cath {
	namespace scan {
		namespace detail {
			class scan_structure_data;
			rep_strider get_this_from_strider(const scan_structure_data &);
			rep_strider get_this_to_strider(const scan_structure_data &);

			/// \brief For a given structure, store the data needed for following up a multi_struc_res_rep_pair match
			///
			/// The interesting data is a list of the single_struc_res_pair neighbours for each
			/// of the multi_struc_res_rep_pair representatives.
			///
			/// For a given scan_stride, the pattern of neighbours is consistent making it easy to
			/// perform pairwise of two reps' neighbours. Conversely, neighbours should not be compared
			/// from scan_structure_data objects generated from different scan_strides.
			///
			/// This also stores data about the strides and number of residues to make it possible
			/// to go back from the rep indices stored in the multi_struc_res_pair to the full
			/// indices of those reps in the original protein.
			///
			/// If there's reason, it might make sense to store a reference (or reference_wrapper) back to the original protein
			class scan_structure_data final {
			private:
				/// \brief A vector of single_struc_res_pair lists, one for each representative multi_struct_res_rep_pair
				///
				/// Each list stores the same pattern of neighbours around the representative res_pairs
				/// (which involves storing dummy single_struc_res_pairs for neighbours that overrun the ends).
				///
				/// The ordering is: to_res_rep_indices are minor; from_res_rep_indices are major
				/// (ie `foreach (from_rep_idx) { foreach (to_rep_idx) { add_entry( from_rep_idx, to_rep_idx ); } }`)
				single_struc_res_pair_list_vec rep_sets;

				/// \brief The total number of residues in the source structure
				index_type num_residues;

				/// \brief TODOCUMENT
				roled_scan_stride the_roled_scan_stride;

				void sanity_check() const;
				size_t get_rep_sets_index_of_res_rep_indices(const res_rep_index_type &,
				                                             const res_rep_index_type &) const;

				scan_structure_data(single_struc_res_pair_list_vec,
				                    const index_type &,
				                    roled_scan_stride);

			public:
				scan_structure_data(const protein &,
				                    const roled_scan_stride &);

				const single_struc_res_pair_list & get_res_pairs_of_rep_indices(const res_rep_index_type &,
				                                                                const res_rep_index_type &) const;
				info_quantity get_info_size() const;
				const index_type & get_num_residues() const;
				const roled_scan_stride & get_roled_scan_stride() const;
			};

			/// \brief Perform sanity checks for consistency in the scan_structure_data
			///
			/// \pre The scan_structure_data should be consistent else an out_of_range_exception will be thrown
			inline void scan_structure_data::sanity_check() const {
				const auto num_from_res_reps = debug_numeric_cast<size_t>(get_num_reps_of_num_residues( get_this_from_strider( *this ), get_num_residues() ) );
				const auto num_to_res_reps   = debug_numeric_cast<size_t>(get_num_reps_of_num_residues( get_this_to_strider  ( *this ), get_num_residues() ) );
				if ( rep_sets.size() != num_from_res_reps * num_to_res_reps ) {
					BOOST_THROW_EXCEPTION(common::out_of_range_exception("The number of rep_sets in scan_structure_data does not match the expected number"));
				}
			}

			/// \brief Calculate the index of the entry in rep_sets that corresponding to the two specified rep indices
			///
			/// NOTE: the inputs are indices in the lists of reps, not indices in the whole protein
			///
			/// Even in a debug build, this doesn't check the result is in range of rep_sets - that can be done
			/// by the range-checked std::vector on the attempt to use the index on rep_sets
			inline size_t scan_structure_data::get_rep_sets_index_of_res_rep_indices(const res_rep_index_type &prm_from_rep_index, ///< The rep res_pair index of the from-residue being queried
			                                                                         const res_rep_index_type &prm_to_rep_index    ///< The rep res_pair index of the to-residue   being queried
			                                                                         ) const {
				const auto num_to_res_reps = debug_numeric_cast<size_t>( get_num_reps_of_num_residues(
					get_this_to_strider( *this ),
					get_num_residues()
				) );
				return ( prm_from_rep_index * num_to_res_reps ) + prm_to_rep_index;
			}

			/// \brief Ctor for building from required data and sanity checking
			///
			/// This is private; a public ctor that ensures consistency is provided and that ctor delegates to this one
			inline scan_structure_data::scan_structure_data(single_struc_res_pair_list_vec  prm_rep_sets,         ///< The lists of neighbours for each of the rep res_pairs
			                                                const index_type               &prm_num_residues,     ///< The total number of residues in the source protein
			                                                roled_scan_stride               prm_roled_scan_stride ///< TODOCUMENT
			                                                ) : rep_sets              { std::move( prm_rep_sets          ) },
			                                                    num_residues          { prm_num_residues                   },
			                                                    the_roled_scan_stride { std::move( prm_roled_scan_stride ) } {
				sanity_check();
			}

			/// \brief Ctor for building from a protein and the relevant roled_scan_stride
			inline scan_structure_data::scan_structure_data(const protein           &prm_protein,          ///< The protein that this scan_structure_data should represent
			                                                const roled_scan_stride &prm_roled_scan_stride ///< TODOCUMENT
			                                                ) : scan_structure_data(
			                                                    	detail::build_rep_sets( prm_protein, prm_roled_scan_stride ),
			                                                    	debug_unwarned_numeric_cast<index_type>( prm_protein.get_length() ),
			                                                    	prm_roled_scan_stride
			                                                    ) {
			}

			/// \brief Get the list of neighbour single_struc_res_pairs relating to the specified from/to rep indices
			inline const single_struc_res_pair_list & scan_structure_data::get_res_pairs_of_rep_indices(const res_rep_index_type &prm_from_rep_index, ///< The rep index of the from residue of interest
			                                                                                            const res_rep_index_type &prm_to_rep_index    ///< The rep index of the to   residue of interest
			                                                                                            ) const {
				return rep_sets[
					get_rep_sets_index_of_res_rep_indices(
						prm_from_rep_index,
						prm_to_rep_index
					)
				];
			}

			/// \brief TODOCUMENT
			inline info_quantity scan_structure_data::get_info_size() const {
				const auto num_bytes = rep_sets.empty()
					? 0
					: rep_sets.size() * common::front( rep_sets ).size() * sizeof( single_struc_res_pair );
				return num_bytes * boost::units::information::bytes;
			}

			/// \brief Getter for the number of residues in the full source protein
			inline const index_type & scan_structure_data::get_num_residues() const {
				return num_residues;
			}

			/// \brief TODOCUMENT
			inline const roled_scan_stride & scan_structure_data::get_roled_scan_stride() const {
				return the_roled_scan_stride;
			}

			/// \brief TODOCUMENT
			inline rep_strider get_this_from_strider(const scan_structure_data &prm_scan_structure_data ///< TODOCUMENT
			                                         ) {
				return get_this_from_strider( prm_scan_structure_data.get_roled_scan_stride() );
			}

			/// \brief TODOCUMENT
			inline rep_strider get_this_to_strider(const scan_structure_data &prm_scan_structure_data ///< TODOCUMENT
			                                       ) {
				return get_this_to_strider  ( prm_scan_structure_data.get_roled_scan_stride() );
			}

			/// \brief Get the list of single_struc_res_pair neighbours associated with
			///        the specified rep specified multi_struc_res_rep_pair
			///
			/// \relates scan_structure_data
			inline const single_struc_res_pair_list & get_res_pairs_of_rep_res_pair(const scan_structure_data      &prm_scan_structure_data, ///< The scan_structure_data to query
			                                                                        const multi_struc_res_rep_pair &prm_res_pair             ///< The rep multi_struc_res_rep_pair for which the neighbours should be extracted
			                                                                        ) {
				return prm_scan_structure_data.get_res_pairs_of_rep_indices(
					prm_res_pair.get_from_res_rep_index(),
					prm_res_pair.get_to_res_rep_index()
				);
			}

			/// \brief Get the full index in the source protein represented by the from-residue in the specified res_pair
			///
			/// \relates scan_structure_data
			inline index_type get_from_res_index(const scan_structure_data      &prm_scan_structure_data, ///< The scan_structure_data associated with the specified res_pair
			                                     const multi_struc_res_rep_pair &prm_res_pair             ///< The res_pair with the from-index for which the full index (in the source protein) should be extracted
			                                     ) {
				return get_index_of_rep_index(
					get_this_from_strider( prm_scan_structure_data ),
					prm_res_pair.get_from_res_rep_index()
				);
			}

			/// \brief Get the full index in the source protein represented by the to-residue in the specified res_pair
			///
			/// \relates scan_structure_data
			inline index_type get_to_res_index(const scan_structure_data      &prm_scan_structure_data, ///< The scan_structure_data associated with the specified res_pair
			                                   const multi_struc_res_rep_pair &prm_res_pair             ///< The res_pair with the to-index for which the full index (in the source protein) should be extracted
			                                   ) {
				return get_index_of_rep_index(
					get_this_to_strider( prm_scan_structure_data ),
					prm_res_pair.get_to_res_rep_index()
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
			///
			/// \todo Provide an overload with the scan_multi_structure_data (correct name?) that uses
			///       prm_res_pair.get_structure_index() to get the correct scan_structure_data and
			///       then calls this
			inline bool are_not_violated_by(const quad_criteria            &prm_criteria,           ///< The criteria to apply
			                                const multi_struc_res_rep_pair &prm_res_pair,           ///< The res_pair to be tested for whether it has any chance of matching another
			                                const scan_structure_data      &prm_scan_structure_data ///< The corresponding scan_structure_data for this res_pair's structure
			                                ) {
				return are_not_violated_by(
					prm_criteria,
					get_from_res_index( prm_scan_structure_data, prm_res_pair ),
					get_to_res_index  ( prm_scan_structure_data, prm_res_pair )
				);
			}

		} // namespace detail
	} // namespace scan
} // namespace cath

#endif
