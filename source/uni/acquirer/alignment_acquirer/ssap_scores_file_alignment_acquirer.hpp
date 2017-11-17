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

#ifndef _CATH_TOOLS_SOURCE_ACQUIRER_ALIGNMENT_ACQUIRER_SSAP_SCORES_FILE_ALIGNMENT_ACQUIRER_H
#define _CATH_TOOLS_SOURCE_ACQUIRER_ALIGNMENT_ACQUIRER_SSAP_SCORES_FILE_ALIGNMENT_ACQUIRER_H

#include <boost/filesystem.hpp>
#include <boost/optional.hpp>

#include "acquirer/alignment_acquirer/alignment_acquirer.hpp"
#include "alignment/align_type_aliases.hpp"
#include "uni/alignment/aln_glue_style.hpp"

namespace cath { namespace align { class alignment; } }

namespace cath {
	namespace align {

		/// \brief TODOCUMENT
		class ssap_scores_file_alignment_acquirer final : public alignment_acquirer {
		private:
			using super = alignment_acquirer;
			boost::filesystem::path ssap_scores_filename;

			/// \brief The approach that should be used for glueing alignments together
			aln_glue_style glue_style = aln_glue_style::SIMPLY;

			std::unique_ptr<alignment_acquirer> do_clone() const final;
			std::pair<alignment, size_size_pair_vec> do_get_alignment_and_spanning_tree(const file::pdb_list &) const final;

		public:
			/// \brief The default for the approach that should be used for glueing alignments together
			static constexpr aln_glue_style DEFAULT_ALN_GLUE_STYLE = aln_glue_style::SIMPLY;

			explicit ssap_scores_file_alignment_acquirer(const boost::filesystem::path &,
			                                             const aln_glue_style & = DEFAULT_ALN_GLUE_STYLE);

			boost::filesystem::path get_ssap_scores_file() const;
			const aln_glue_style & get_glue_style() const;
		};

		size_size_alignment_tuple_vec get_spanning_alignments(const file::pdb_list &,
		                                                      const str_vec &,
		                                                      const size_size_doub_tpl_vec &,
		                                                      const boost::filesystem::path &,
		                                                      const ostream_ref_opt & = boost::none);

		std::pair<alignment, size_size_pair_vec> build_multi_alignment(const file::pdb_list &,
		                                                               const str_vec &,
		                                                               const size_size_doub_tpl_vec &,
		                                                               const boost::filesystem::path &,
		                                                               const aln_glue_style &,
		                                                               const ostream_ref_opt & = boost::none);

	} // namespace align
} // namespace cath

#endif
