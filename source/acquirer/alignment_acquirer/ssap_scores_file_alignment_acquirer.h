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

#ifndef SSAP_SCORES_FILE_ALIGNMENT_ACQUIRER_H_INCLUDED
#define SSAP_SCORES_FILE_ALIGNMENT_ACQUIRER_H_INCLUDED

#include <boost/filesystem.hpp>

#include "acquirer/alignment_acquirer/alignment_acquirer.h"
#include "alignment/align_type_aliases.h"

namespace cath { namespace align { class alignment; } }

namespace cath {
	namespace opts {

		/// \brief TODOCUMENT
		class ssap_scores_file_alignment_acquirer final : public alignment_acquirer {
		private:
			using super = alignment_acquirer;
			boost::filesystem::path ssap_scores_file;

			virtual std::unique_ptr<alignment_acquirer> do_clone() const override final;
			virtual std::pair<align::alignment, sup::superpose_orderer> do_get_alignment_and_orderer(const file::pdb_list &) const override final;

			align::size_size_alignment_tuple_vec get_spanning_alignments(const boost::filesystem::path &,
			                                                             const str_vec &,
			                                                             const file::pdb_list &,
			                                                             const size_size_pair_vec &,
			                                                             std::ostream &) const;

		public:
			ssap_scores_file_alignment_acquirer(const boost::filesystem::path &);
			virtual ~ssap_scores_file_alignment_acquirer() noexcept = default;

			boost::filesystem::path get_ssap_scores_file() const;
		};

	}
}

#endif
