/// \file
/// \brief The ssap_scores_file class header

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

#ifndef _CATH_TOOLS_SOURCE_FILE_SSAP_SCORES_FILE_SSAP_SCORES_FILE_H
#define _CATH_TOOLS_SOURCE_FILE_SSAP_SCORES_FILE_SSAP_SCORES_FILE_H

#include <boost/filesystem/path.hpp>

#include "common/type_aliases.hpp"
#include "file/file_type_aliases.hpp"

#include <iosfwd>

namespace cath {
	namespace file {

		/// \brief TODOCUMENT
		class ssap_scores_file final {
		private:
			ssap_scores_file() = delete;

		public:
			static ssap_scores_entry_vec parse_ssap_scores_file_simple(std::istream &);

			static ssap_scores_entry_vec parse_ssap_scores_file_simple(const std::string &);

			static ssap_scores_entry_vec parse_ssap_scores_file_simple(const boost::filesystem::path &);

			static std::pair<str_vec, size_size_pair_doub_map> parse_ssap_scores_file(std::istream &);

			static std::pair<str_vec, size_size_pair_doub_map> parse_ssap_scores_file(const boost::filesystem::path &);
		};

		str_str_pair_bool_map make_arbitrary_is_positive_data(const ssap_scores_entry_vec &);

	} // namespace file
} // namespace cath

#endif
