/// \file
/// \brief The alignment_row class header

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

#ifndef CATH_TOOLS_SOURCE_CT_UNI_CATH_ALIGNMENT_ALIGNMENT_ROW_HPP
#define CATH_TOOLS_SOURCE_CT_UNI_CATH_ALIGNMENT_ALIGNMENT_ROW_HPP

#include "cath/alignment/align_type_aliases.hpp"
#include "cath/common/type_aliases.hpp"

// clang-format off
namespace cath::align { class alignment; }
// clang-format on

namespace cath::align {

	/// \brief Represent a row in an alignment (ie one index across several entries being aligned)
	///
	/// This is a useful for manipulating alignment rows (eg appending them to the end of alignments)
	///
	/// The value for each entry is stored as an optional<aln_posn_type> and this is revealed in some of the interface.
	/// An absent position is represented by a value of nullopt; a present position is indicated by a value of the
	/// index in the thing being aligned (offset 0).
	///
	/// NOTE: there isn't any suggestion that this is how rows are stored in the alignment class (it isn't).
	class alignment_row final {
	private:
		/// \brief The positions in the row (with absent positions represented as nullopt)
		aln_posn_opt_vec positions;

		void sanity_check_entry(const size_t &) const;
		void check_has_position_of_entry(const size_t &) const;

	public:
		explicit alignment_row(aln_posn_opt_vec);

		[[nodiscard]] size_t       num_entries() const;
		[[nodiscard]] aln_posn_opt position_of_entry( const size_t & ) const;
	};

	bool has_position_of_entry(const alignment_row &,
	                           const size_t &);

	aln_posn_type get_position_of_entry(const alignment_row &,
	                                    const size_t &);

	bool any_entries_present(const alignment_row &);

	bool any_entries_present(const alignment &,
	                         const size_vec &,
	                         const size_t &);

	alignment_row get_row_of_entries_of_alignment(const alignment &,
	                                              const size_vec &,
	                                              const size_t &);

	alignment_row get_row_of_alignment(const alignment &,
	                                   const size_t &);

	alignment_row make_empty_aln_row(const size_t &);

	alignment_row make_row_with_single_value(const size_t &,
	                                         const size_t &,
	                                         const size_t &);

	void append_row_with_single_value(alignment &,
	                                  const size_t &,
	                                  const size_t &);

	alignment_row copy_aln_row_entry(const alignment_row &,
	                                 const size_t &,
	                                 const alignment_row &,
	                                 const size_t &);

	alignment_row glue_aln_row_onto_empties(const size_t &,
	                                        const size_t &,
	                                        const alignment_row &,
	                                        const size_t &);

	alignment_row glue_empties_onto_row(const alignment_row &,
	                                    const size_t &);

	alignment_row glue_aln_rows_together(const alignment_row &,
	                                     const size_t &,
	                                     const alignment_row &,
	                                     const size_t &);

	alignment_row join(const alignment_row &,
	                   const alignment_row &);

	void set_positions_of_entries_from_row(aln_posn_opt_vec &,
	                                       const alignment_row &,
	                                       const size_vec &);

	alignment_row weave(const alignment_row &,
	                    const size_vec &,
	                    const alignment_row &,
	                    const size_vec &);

	alignment_row remove_entry_from_row( const alignment_row &, const size_t & );

	aln_posn_opt_vec get_has_posns_and_posns(const alignment_row &);

} // namespace cath::align

#endif // CATH_TOOLS_SOURCE_CT_UNI_CATH_ALIGNMENT_ALIGNMENT_ROW_HPP
