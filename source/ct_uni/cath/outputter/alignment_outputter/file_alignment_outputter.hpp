/// \file
/// \brief The file_alignment_outputter class header

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

#ifndef CATH_TOOLS_SOURCE_CT_UNI_CATH_OUTPUTTER_ALIGNMENT_OUTPUTTER_FILE_ALIGNMENT_OUTPUTTER_HPP
#define CATH_TOOLS_SOURCE_CT_UNI_CATH_OUTPUTTER_ALIGNMENT_OUTPUTTER_FILE_ALIGNMENT_OUTPUTTER_HPP

#include <filesystem>

#include "cath/common/clone/clone_ptr.hpp"
#include "cath/outputter/alignment_outputter/alignment_outputter.hpp"

namespace cath::opts {

	/// \brief TODOCUMENT
	class file_alignment_outputter final : public alignment_outputter {
	private:
		/// \brief TODOCUMENT
		const ::std::filesystem::path output_file;

		/// \brief TODOCUMENT
		common::clone_ptr<alignment_outputter> ostream_alignment_outputter_ptr;

		[[nodiscard]] std::unique_ptr<alignment_outputter> do_clone() const final;
		void do_output_alignment(const align::alignment_context &,
		                         std::ostream &) const final;
		[[nodiscard]] bool                                 do_involves_display_spec() const final;
		[[nodiscard]] std::string                          do_get_name() const final;

	  public:
		file_alignment_outputter( ::std::filesystem::path, const alignment_outputter & );
	};

} // namespace cath::opts

#endif // CATH_TOOLS_SOURCE_CT_UNI_CATH_OUTPUTTER_ALIGNMENT_OUTPUTTER_FILE_ALIGNMENT_OUTPUTTER_HPP
