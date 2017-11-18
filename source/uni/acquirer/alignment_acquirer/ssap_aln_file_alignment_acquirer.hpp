/// \file
/// \brief The ssap_aln_file_alignment_acquirer class header

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

#ifndef _CATH_TOOLS_SOURCE_ACQUIRER_ALIGNMENT_ACQUIRER_SSAP_ALN_FILE_ALIGNMENT_ACQUIRER_H
#define _CATH_TOOLS_SOURCE_ACQUIRER_ALIGNMENT_ACQUIRER_SSAP_ALN_FILE_ALIGNMENT_ACQUIRER_H

#include <boost/filesystem.hpp>

#include "acquirer/alignment_acquirer/alignment_acquirer.hpp"

namespace cath {
	namespace align {

		/// \brief TODOCUMENT
		class ssap_aln_file_alignment_acquirer final : public alignment_acquirer {
		private:
			using super = alignment_acquirer;

			boost::filesystem::path ssap_alignment_file;

			std::unique_ptr<alignment_acquirer> do_clone() const final;
			std::pair<alignment, size_size_pair_vec> do_get_alignment_and_spanning_tree(const file::strucs_context &) const final;

		public:
			explicit ssap_aln_file_alignment_acquirer(const boost::filesystem::path &);

			boost::filesystem::path get_ssap_alignment_file() const;
		};

	} // namespace align
} // namespace cath

#endif
