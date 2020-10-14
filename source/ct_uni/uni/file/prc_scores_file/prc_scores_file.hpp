/// \file
/// \brief The prc_scores_file class header

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

#ifndef _CATH_TOOLS_SOURCE_UNI_FILE_PRC_SCORES_FILE_PRC_SCORES_FILE_HPP
#define _CATH_TOOLS_SOURCE_UNI_FILE_PRC_SCORES_FILE_PRC_SCORES_FILE_HPP

#include <boost/filesystem/path.hpp>

#include "file/file_type_aliases.hpp"

namespace cath {
	namespace file {

		/// \brief TODOCUMENT
		class prc_scores_file final {
		private:
			prc_scores_file() = delete;

		public:
			static prc_scores_entry_vec remove_duplicates(const prc_scores_entry_vec &);

			static prc_scores_entry_vec parse_prc_scores_file(std::istream &);

			static prc_scores_entry_vec parse_prc_scores_file(const std::string &);

			static prc_scores_entry_vec parse_prc_scores_file(const boost::filesystem::path &);

			static prc_scores_entry_vec parse_prc_scores_file_fancy(std::istream &);

			static prc_scores_entry_vec parse_prc_scores_file_fancy(const std::string &);

			static prc_scores_entry_vec parse_prc_scores_file_fancy(const boost::filesystem::path &);
		};

	} // namespace file
} // namespace cath

#endif
