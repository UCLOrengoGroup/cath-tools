/// \file
/// \brief The do_the_ssaps_alignment_acquirer class header

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

#ifndef _CATH_TOOLS_SOURCE_UNI_ACQUIRER_ALIGNMENT_ACQUIRER_DO_THE_SSAPS_ALIGNMENT_ACQUIRER_HPP
#define _CATH_TOOLS_SOURCE_UNI_ACQUIRER_ALIGNMENT_ACQUIRER_DO_THE_SSAPS_ALIGNMENT_ACQUIRER_HPP

#include <boost/filesystem/path.hpp>
#include <boost/optional.hpp>

#include "cath/acquirer/alignment_acquirer/alignment_acquirer.hpp"
#include "cath/common/path_type_aliases.hpp"

namespace cath { namespace align { class alignment; } }

namespace cath { namespace align {

	/// \brief Acquire the alignment by performing any necessary cath-ssaps in
	///        a temp directory and then using ssap_scores_file_alignment_acquirer
	class do_the_ssaps_alignment_acquirer final : public alignment_acquirer {
	private:
		using super = alignment_acquirer;

		/// \brief Where the magic shall happen
		path_opt directory_of_joy;

		std::unique_ptr<alignment_acquirer> do_clone() const final;
		bool do_requires_backbone_complete_input() const final;
		std::pair<alignment, size_size_pair_vec> do_get_alignment_and_spanning_tree(const file::strucs_context &,
		                                                                            const align_refining &) const final;

	public:
		explicit do_the_ssaps_alignment_acquirer(const path_opt & = boost::none);

		const path_opt & get_directory_of_joy() const;

		static boost::filesystem::path make_temp_dir_for_doing_ssaps(const file::strucs_context &);
	};

} } // namespace cath::align

#endif
