/// \file
/// \brief The hmmer_scores_file class header

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

#ifndef _CATH_TOOLS_SOURCE_CT_UNI_CATH_FILE_HMMER_SCORES_FILE_HMMER_SCORES_FILE_HPP
#define _CATH_TOOLS_SOURCE_CT_UNI_CATH_FILE_HMMER_SCORES_FILE_HMMER_SCORES_FILE_HPP

#include <filesystem>

#include "cath/file/file_type_aliases.hpp"
#include "cath/file/hmmer_scores_file/hmmer_name_handling.hpp"

namespace cath::file {

	/// \brief TODOCUMENT
	class hmmer_scores_file final {
	private:
		hmmer_scores_file() = delete;

	public:
		static hmmer_scores_entry_vec remove_duplicates(const hmmer_scores_entry_vec &);

		static hmmer_scores_entry_vec parse_hmmer_scores_file(std::istream &,
		                                                      const hmmer_name_handling & = hmmer_name_handling::STRIP);

		static hmmer_scores_entry_vec parse_hmmer_scores_file(const ::std::filesystem::path &,
		                                                      const hmmer_name_handling & = hmmer_name_handling::STRIP);

//		static hmmer_scores_entry_vec parse_hmmer_scores_file_fancy(std::istream &);

//		static hmmer_scores_entry_vec parse_hmmer_scores_file_fancy(const ::std::filesystem::path &);
	};

} // namespace cath::file

#endif // _CATH_TOOLS_SOURCE_CT_UNI_CATH_FILE_HMMER_SCORES_FILE_HMMER_SCORES_FILE_HPP
