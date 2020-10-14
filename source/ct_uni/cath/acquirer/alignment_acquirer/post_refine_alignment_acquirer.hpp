/// \file
/// \brief The post_refine_alignment_acquirer class header

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

#ifndef _CATH_TOOLS_SOURCE_UNI_ACQUIRER_ALIGNMENT_ACQUIRER_POST_REFINE_ALIGNMENT_ACQUIRER_HPP
#define _CATH_TOOLS_SOURCE_UNI_ACQUIRER_ALIGNMENT_ACQUIRER_POST_REFINE_ALIGNMENT_ACQUIRER_HPP

#include <boost/filesystem.hpp>

#include "cath/acquirer/alignment_acquirer/alignment_acquirer.hpp"

namespace cath {
	namespace align {

		/// \brief Abstract alignment_acquirer that refines (and rescores) the alignment
		///        after the concrete alignment_acquirer has done its work
		class post_refine_alignment_acquirer : public alignment_acquirer {
		private:
			using super = alignment_acquirer;

			virtual std::unique_ptr<alignment_acquirer> do_clone() const = 0;

			virtual std::pair<alignment, size_size_pair_vec> do_get_alignment_and_spanning_tree(const file::strucs_context &) const = 0;

			std::pair<alignment, size_size_pair_vec> do_get_alignment_and_spanning_tree(const file::strucs_context &,
			                                                                            const align_refining &) const final;
		};

	} // namespace align
} // namespace cath

#endif
