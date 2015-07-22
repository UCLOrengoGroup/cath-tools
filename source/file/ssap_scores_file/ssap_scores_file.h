/// \file
/// \brief The ssap_scores_file class header

/// \copyright
/// CATH Binaries - Protein structure comparison tools such as SSAP and SNAP
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

#ifndef SSAP_SCORES_FILE_H_INCLUDED
#define SSAP_SCORES_FILE_H_INCLUDED

#include <boost/filesystem/path.hpp>

#include "common/type_aliases.h"

#include <iosfwd>

namespace cath {
	namespace file {

		/// \brief TODOCUMENT
		class ssap_scores_file final {
		private:
			ssap_scores_file() = delete;

		public:
			static std::pair<str_vec, size_size_pair_doub_map> parse_ssap_scores_file(std::istream &);

			static std::pair<str_vec, size_size_pair_doub_map> parse_ssap_scores_file(const boost::filesystem::path &);
		};

	}
}

#endif