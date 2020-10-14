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

#ifndef _CATH_TOOLS_SOURCE_UNI_ACQUIRER_ALIGNMENT_ACQUIRER_SSAP_SCORES_FILE_ALIGNMENT_ACQUIRER_HPP
#define _CATH_TOOLS_SOURCE_UNI_ACQUIRER_ALIGNMENT_ACQUIRER_SSAP_SCORES_FILE_ALIGNMENT_ACQUIRER_HPP

#include <boost/filesystem.hpp>
#include <boost/optional.hpp>

#include "cath/acquirer/alignment_acquirer/alignment_acquirer.hpp"
#include "cath/alignment/align_type_aliases.hpp"
#include "cath/alignment/aln_glue_style.hpp"

namespace cath { namespace align { class alignment; } }
namespace cath { namespace file { class pdb_list; } }

namespace cath {
	namespace align {

		/// \brief TODOCUMENT
		class ssap_scores_file_alignment_acquirer final : public alignment_acquirer {
		private:
			using super = alignment_acquirer;
			boost::filesystem::path ssap_scores_filename;

			std::unique_ptr<alignment_acquirer> do_clone() const final;
			bool do_requires_backbone_complete_input() const final;
			std::pair<alignment, size_size_pair_vec> do_get_alignment_and_spanning_tree(const file::strucs_context &,
			                                                                            const align_refining &) const final;

		public:
			explicit ssap_scores_file_alignment_acquirer(const boost::filesystem::path &);

			boost::filesystem::path get_ssap_scores_file() const;
		};

		std::pair<alignment, size_size_pair_vec> build_multi_alignment(const file::pdb_list &,
		                                                               const str_vec &,
		                                                               const size_size_doub_tpl_vec &,
		                                                               const boost::filesystem::path &,
		                                                               const aln_glue_style &,
		                                                               const ostream_ref_opt & = boost::none);

	} // namespace align
} // namespace cath

#endif
