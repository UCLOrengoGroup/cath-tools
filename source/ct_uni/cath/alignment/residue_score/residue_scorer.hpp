/// \file
/// \brief The residue_scorer class header

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

#ifndef CATH_TOOLS_SOURCE_CT_UNI_CATH_ALIGNMENT_RESIDUE_SCORE_RESIDUE_SCORER_HPP
#define CATH_TOOLS_SOURCE_CT_UNI_CATH_ALIGNMENT_RESIDUE_SCORE_RESIDUE_SCORER_HPP

#include <filesystem>

// clang-format off
namespace cath { class protein_list; }
namespace cath::align { class alignment; }
namespace cath::align { class alignment_residue_scores; }
// clang-format on

namespace cath::align {

	/// \brief TODOCUMENT
	class residue_scorer final {
	  private:
		// bool renormalise;
		// bool from_to_or_both;
		// bool best_twenty;
		// bool atom_policy_ca_or_cb;
		// bool score_function;

		// alignment_residue_scores should store the number of entries and the number of present entries in the position;

	  public:
		[[nodiscard]] alignment_residue_scores get_alignment_residue_scores( const alignment &, const protein_list & ) const;
	};

	void score_alignment(const residue_scorer &,
	                     alignment &,
	                     const protein_list &);

	alignment score_alignment_copy(const residue_scorer &,
	                               alignment,
	                               const protein_list &);

	alignment read_and_rescore_fasta_alignment(const ::std::filesystem::path &,
	                                           const protein_list &,
	                                           const residue_scorer &,
	                                           std::ostream &);

} // namespace cath::align

#endif // CATH_TOOLS_SOURCE_CT_UNI_CATH_ALIGNMENT_RESIDUE_SCORE_RESIDUE_SCORER_HPP
