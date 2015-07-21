/// \file
/// \brief The alignment class header

/// \copyright
/// CATH Binaries - Protein structure comparison tools such as SSAP and SNAP
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

#ifndef ALIGNMENT_H_INCLUDED
#define ALIGNMENT_H_INCLUDED

#include <boost/operators.hpp>
#include <boost/optional.hpp>

#include "alignment/align_type_aliases.h"
#include "alignment/residue_score/alignment_residue_scores.h"
#include "common/type_aliases.h"

#include <cstddef>
#include <iosfwd>

namespace cath { namespace align { class alignment_row; } }

namespace cath {
	namespace align {

		/// \brief TODOCUMENT
		///
		/// There is one score for each position (which is how both SSAP and CORA score alignments)
		///
		/// \todo Consider making this class enforce that all positions within an entry are strictly increasing.
		///
		/// At present, the alignment demands that all scores are added together at the end and that the alignment is not
		/// modified after that. This makes it simple to check that the whole alignment is_scored() but it may turn out to be
		/// a touch too strict and hence require loosening.
		///
		/// Nomenclature:
		///  - index    : the index into the list of equivalences, starting from 0
		///               (ie the row index in a SSAP/CORA alignment file
		///                or the column index in a FASTA alignment)
		///  - entry    : the index of the entry being aligned, starting from 0
		///               (ie the column index in a SSAP/CORA alignment file
		///                or the row index in a FASTA alignment)
		///  - position : the contents of a value in the alignment, ie the position in the alignee
		///
		/// It's useful to allow an alignment with one entry so that cath-superpose doesn't
		/// have to treat that as a special case.
		class alignment final : private boost::equality_comparable<alignment> {
		public:
			/// \brief TODOCUMENT
			using size_type = size_vec_vec::size_type;

		private:
//			bool_deq_vec     position_presences;
//			aln_posn_vec_vec old_positions;

			/// \brief TODOCUMENT
			opt_aln_posn_vec_vec positions;

			/// \brief TODOCUMENT
			size_type logical_length;

			/// \brief TODOCUMENT
			boost::optional<alignment_residue_scores> new_scores;

			void check_scored() const;
			void check_not_scored() const;

			void check_entry_in_range(const size_type &) const;
			void check_index_in_range(const size_type &) const;

			size_type reserved_length() const;

		public:
			explicit alignment(const size_type &);
			explicit alignment(const opt_aln_posn_vec_vec &);

			void reserve(const size_type &);
			size_type num_entries() const;
			size_type length() const;

			bool is_scored() const;
			const alignment_residue_scores & get_alignment_residue_scores() const;
//			double get_score_of_index(const size_type &) const;

			opt_aln_posn position_of_entry_of_index(const size_type &,
			                                        const size_type &) const;

			void set_position_value(const size_type &,
			                        const size_type &,
			                        const aln_posn_type &);

			void set_scores(const alignment_residue_scores &);
	//		void set_non_raw_scores(const doub_vec &);

			/// \brief A constant for the number of entries in a pair alignment
			static constexpr size_type NUM_ENTRIES_IN_PAIR_ALIGNMENT = 2;
			/// \brief A constant for the index of the A entry (0) in a pair alignment
			static constexpr size_type PAIR_A_IDX = 0;
			/// \brief A constant for the index of the B entry (1) in a pair alignment
			static constexpr size_type PAIR_B_IDX = 1;
		};

		opt_size_size_pair first_non_consecutive_entry_positions(const alignment &);

		void check_entry_positions_are_consecutive(const alignment &);

		bool has_position_of_entry_of_index(const alignment &,
		                                    const size_t &,
		                                    const size_t &);

		bool has_position_of_both_entries_of_index(const alignment &,
		                                           const size_t &,
		                                           const size_t &,
		                                           const size_t &);

		bool has_position_of_all_entries_of_index(const alignment &,
		                                          const size_t &);

		aln_posn_type get_position_of_entry_of_index(const alignment &,
		                                             const size_t &,
		                                             const size_t &);

		/// \brief TODOCUMENT
		using opt_aln_size = boost::optional<alignment::size_type>;

		bool operator==(const alignment &,
		                const alignment &);

		std::ostream & operator<<(std::ostream &,
		                          const alignment &);

		size_vec present_positions_of_index(const alignment &,
		                                    const alignment::size_type &);

		size_t num_present_positions_of_index(const alignment &,
		                                      const alignment::size_type &);

		size_vec num_present_positions_by_index(const alignment &);

		size_vec present_positions_of_entry(const alignment &,
		                                    const alignment::size_type &);

		size_t num_present_positions_of_entry(const alignment &,
		                                      const alignment::size_type &);

		size_vec num_present_positions_by_entry(const alignment &);

		opt_aln_size get_index_of_first_present_position_of_entry(const alignment &,
		                                                          const alignment::size_type &);

		opt_aln_size get_index_of_last_present_position_of_entry(const alignment &,
		                                                         const alignment::size_type &);

		size_t num_present_positions_of_both_entries(const alignment &,
		                                             const alignment::size_type &,
		                                             const alignment::size_type &);

		size_t num_present_positions_of_all_entries(const alignment &);

		opt_aln_size get_index_of_first_present_position_of_both_entries(const alignment &,
		                                                                 const alignment::size_type &,
		                                                                 const alignment::size_type &);

		opt_aln_size get_index_of_last_present_position_of_both_entries(const alignment &,
		                                                                const alignment::size_type &,
		                                                                const alignment::size_type &);

		opt_aln_posn get_first_present_position_of_entry(const alignment &,
		                                                 const alignment::size_type &);

		opt_aln_posn get_last_present_position_of_entry(const alignment &,
		                                                const alignment::size_type &);

		opt_aln_posn get_max_last_present_position(const alignment &);

		float_score_type get_score_of_entry_and_index(const alignment &,
		                                              const size_t &,
		                                              const size_t &);

		float_score_type get_mean_score_of_index(const alignment &,
		                                         const size_t &);

		float_score_type get_total_score_of_index(const alignment &,
		                                          const size_t &);

		float_score_vec get_mean_score_by_index(const alignment &);

		float_score_vec get_total_score_by_index(const alignment &);

		float_score_vec get_total_score_or_num_positions_by_index(const alignment &);

		void set_scores(alignment &,
		                const opt_score_vec_vec &);

		alignment set_scores_copy(alignment,
		                          const opt_score_vec_vec &);

		void set_empty_scores(alignment &);

		alignment set_empty_scores_copy(alignment);

		alignment make_single_alignment(const size_t &);

		alignment alignment_offset_1_factory(opt_aln_posn_vec_vec);

		alignment make_permuted_alignment(const alignment &,
		                                  const size_vec &);

		void append_row(alignment &,
		                const alignment_row &);

		void append_row(alignment &,
		                const opt_aln_posn_vec &);

		void set_position_value(alignment &,
		                        const alignment::size_type &,
		                        const alignment::size_type &,
		                        const size_t &);

	}
}
#endif
