/// \file
/// \brief The ssap_scores_file_alignment_acquirer class header

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

#ifndef CATH_TOOLS_SOURCE_CT_UNI_CATH_ACQUIRER_ALIGNMENT_ACQUIRER_SSAP_SCORES_FILE_ALIGNMENT_ACQUIRER_HPP
#define CATH_TOOLS_SOURCE_CT_UNI_CATH_ACQUIRER_ALIGNMENT_ACQUIRER_SSAP_SCORES_FILE_ALIGNMENT_ACQUIRER_HPP

#include <filesystem>
#include <optional>

#include "cath/acquirer/alignment_acquirer/alignment_acquirer.hpp"
#include "cath/alignment/align_type_aliases.hpp"
#include "cath/alignment/aln_glue_style.hpp"

// clang-format off
namespace cath::align { class alignment; }
namespace cath::file { class pdb_list; }
// clang-format on

namespace cath::align {

	/// \brief TODOCUMENT
	class ssap_scores_file_alignment_acquirer final : public alignment_acquirer {
	  private:
		using super = alignment_acquirer;
		::std::filesystem::path ssap_scores_filename;

		[[nodiscard]] std::unique_ptr<alignment_acquirer>      do_clone() const final;
		[[nodiscard]] bool                                     do_requires_backbone_complete_input() const final;
		[[nodiscard]] std::pair<alignment, size_size_pair_vec> do_get_alignment_and_spanning_tree(
		  const file::strucs_context &,
		  const align_refining &,
		  const ostream_ref_opt & = ::std::nullopt ) const final;

	  public:
		explicit ssap_scores_file_alignment_acquirer( ::std::filesystem::path );

		[[nodiscard]] ::std::filesystem::path get_ssap_scores_file() const;
	};

	std::pair<alignment, size_size_pair_vec> build_multi_alignment(const file::pdb_list &,
	                                                               const str_vec &,
	                                                               const size_size_doub_tpl_vec &,
	                                                               const ::std::filesystem::path &,
	                                                               const aln_glue_style &,
	                                                               const ostream_ref_opt & = ::std::nullopt);

} // namespace cath::align

#endif // CATH_TOOLS_SOURCE_CT_UNI_CATH_ACQUIRER_ALIGNMENT_ACQUIRER_SSAP_SCORES_FILE_ALIGNMENT_ACQUIRER_HPP
