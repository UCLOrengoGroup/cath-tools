/// \file
/// \brief The alignment_split_mapping class header

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

#ifndef _CATH_TOOLS_SOURCE_ALIGNMENT_REFINER_DETAIL_ALIGNMENT_SPLIT_MAPPING_H
#define _CATH_TOOLS_SOURCE_ALIGNMENT_REFINER_DETAIL_ALIGNMENT_SPLIT_MAPPING_H

#include "alignment/alignment.hpp"
#include "alignment/refiner/detail/alignment_split_half.hpp"
#include "exception/invalid_argument_exception.hpp"

#include <iostream>

namespace cath {
	namespace align {
		namespace detail {
			class alignment_split;
		} // namespace detail
	} // namespace align
} // namespace cath

namespace cath {
	namespace align {
		namespace detail {

			/// \brief TODOCUMENT
			///
			/// Bear in mind that split that can be aligned the fastest is the one that puts one entry in one group
			/// and all the others. This only has \f$ (n-1) $ inter-residue comparisons to perform per index.
			///
			/// Much worse is equally splitting, which requires \f$ \frac{n^2}{4} $ comparisons. That's \f$\ frac{n}{4} . \frac{n}{n-1} \ge \frac{n}{4} $
			/// as many comparisons as for the single case above
			class alignment_split_mapping final {
			private:
				/// \brief TODOCUMENT
				size_t orig_length;

				/// \brief TODOCUMENT
				size_t orig_num_entries;

				/// \brief TODOCUMENT
				size_t did_insert_entries;


				/// \brief TODOCUMENT
				alignment local_aln;


				/// \brief TODOCUMENT
				///
				/// Since some rows may have been inserted to fill in missing residues, not all indices
				/// will map back to original alignment positions
				size_opt_vec idx_of_orig_aln_idx;

				/// \brief TODOCUMENT
				size_vec orig_aln_entries;

				/// \brief TODOCUMENT
				size_vec_vec index_of_pdb_res_index;

			public:
				alignment_split_mapping(const alignment &,
				                        const size_set &,
				                        const size_vec &);

				size_t orig_aln_length() const;
				size_t orig_aln_num_entries() const;
				bool inserted_entries() const;

				size_t length() const;
				size_t num_entries() const;

				size_opt index_of_orig_aln_index(const size_t &) const;
				size_opt entry_of_orig_aln_entry(const size_t &) const;
				size_t orig_aln_entry_of_entry(const size_t &) const;

				aln_posn_opt position_of_entry_of_index(const size_t &,
				                                        const size_t &) const;

				size_t index_of_protein_index(const size_t &,
				                              const size_t &) const;
			};

			/// \brief TODOCUMENT
			inline size_t alignment_split_mapping::index_of_protein_index(const size_t &arg_entry, ///< TODOCUMENT
			                                                              const size_t &arg_index  ///< TODOCUMENT
			                                                              ) const {
#ifndef NDEBUG
				if ( arg_entry >= index_of_pdb_res_index.size() ) {
					BOOST_THROW_EXCEPTION(cath::common::invalid_argument_exception("Entry is out of range"));
				}
#endif
				const size_vec &pdb_res_indices = index_of_pdb_res_index[ arg_entry ];
#ifndef NDEBUG
				if ( arg_index >= pdb_res_indices.size() ) {
			//#warning This bit of code should be removed
					std::cerr << "About to throw exception \"Index is out of range\"" << std::endl;
					std::cerr << "arg_entry is : " << arg_entry << std::endl;
					std::cerr << "arg_index is : " << arg_index << std::endl;
					std::cerr << "There are " << index_of_pdb_res_index.size() << " entries in index_of_pdb_res_index : " << std::endl;
					for (const size_vec &pdb_res_indices_local : index_of_pdb_res_index) {
						std::cerr << "\t" << pdb_res_indices_local.size() << "\tentries, from\t" << pdb_res_indices_local.front() << "\tto\t" << pdb_res_indices_local.back() << std::endl;
					}
			//		sleep(999999);
					BOOST_THROW_EXCEPTION(cath::common::invalid_argument_exception("Index is out of range"));
				}
#endif
				return pdb_res_indices[ arg_index ];
			}

			std::ostream & operator<<(std::ostream &,
			                          const alignment_split_mapping &);

			size_vec get_orig_entries(const alignment_split_mapping &);

			alignment_row get_row_of_alignment(const alignment_split_mapping &,
			                                   const size_t &);

			bool has_position_of_entry_of_index(const alignment_split_mapping &,
			                                    const size_t &,
			                                    const size_t &);

			size_t get_position_of_entry_of_index(const alignment_split_mapping &,
			                                      const size_t &,
			                                      const size_t &);

//			size_vec present_orig_aln_entries_of_orig_aln_index(const alignment_split_mapping &,
//			                                                    const size_t &);

			size_vec present_orig_aln_entries_of_index(const alignment_split_mapping &,
			                                           const size_t &);

			alignment_split_mapping make_alignment_split_mapping(const alignment &,
			                                                     const alignment_split &,
			                                                     const alignment_split_half &,
			                                                     const size_vec &);

			alignment build_alignment(const alignment &,
			                          const alignment_split_mapping &,
			                          const alignment_split_mapping &);

		} // namespace detail
	} // namespace align
} // namespace cath

#endif
