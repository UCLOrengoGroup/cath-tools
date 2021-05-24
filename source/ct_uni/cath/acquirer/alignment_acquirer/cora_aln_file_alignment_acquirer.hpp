/// \file
/// \brief The cora_aln_file_alignment_acquirer class header

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

#ifndef _CATH_TOOLS_SOURCE_CT_UNI_CATH_ACQUIRER_ALIGNMENT_ACQUIRER_CORA_ALN_FILE_ALIGNMENT_ACQUIRER_HPP
#define _CATH_TOOLS_SOURCE_CT_UNI_CATH_ACQUIRER_ALIGNMENT_ACQUIRER_CORA_ALN_FILE_ALIGNMENT_ACQUIRER_HPP

#include <filesystem>

#include "cath/acquirer/alignment_acquirer/post_refine_alignment_acquirer.hpp"

namespace cath::align {

	/// \brief TODOCUMENT
	class cora_aln_file_alignment_acquirer final : public post_refine_alignment_acquirer {
	private:
		using super = post_refine_alignment_acquirer;

		::std::filesystem::path cora_alignment_file;

		[[nodiscard]] std::unique_ptr<alignment_acquirer>      do_clone() const final;
		[[nodiscard]] bool                                     do_requires_backbone_complete_input() const final;
		[[nodiscard]] std::pair<alignment, size_size_pair_vec> do_get_alignment_and_spanning_tree(
		  const file::strucs_context &,
		  const ostream_ref_opt & = ::std::nullopt ) const final;

	  public:
		explicit cora_aln_file_alignment_acquirer( ::std::filesystem::path );

		[[nodiscard]] ::std::filesystem::path get_cora_alignment_file() const;
	};

} // namespace cath::align

#endif // _CATH_TOOLS_SOURCE_CT_UNI_CATH_ACQUIRER_ALIGNMENT_ACQUIRER_CORA_ALN_FILE_ALIGNMENT_ACQUIRER_HPP
