/// \file
/// \brief The path_or_istream header

/// \copyright
/// Tony Lewis's Common C++ Library Code (here imported into the CATH Tools project and then tweaked, eg namespaced in cath)
/// Copyright (C) 2007, Tony Lewis
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

#ifndef _CATH_TOOLS_SOURCE_SRC_COMMON_COMMON_FILE_PATH_OR_ISTREAM_HPP
#define _CATH_TOOLS_SOURCE_SRC_COMMON_COMMON_FILE_PATH_OR_ISTREAM_HPP

#include <boost/filesystem/path.hpp>
#include <boost/optional.hpp>

#include "common/file/open_fstream.hpp"

#include <fstream>
#include <functional>

namespace cath {
	namespace common {

		/// \brief Type alias for a reference_wrapper of istream
		using istream_ref = std::reference_wrapper<std::istream>;

		/// \brief Type alias for an optional reference_wrapper of istream
		using ostream_ref_opt = boost::optional<istream_ref>;

		/// \brief Handle opening a file and providing an istream, with support for
		///        a special flag which indicates input from the istream optionally specified on construction
		///
		/// Note: this has substantial overlap with ofstream_list and could perhaps share a common implementation
		class path_or_istream final {
		private:
			/// \brief An optional special istream from which input can be read (usually stdin)
			ostream_ref_opt standard_instream;

			/// \brief A flag that can be used when passing a path to indicate input should be read from the standard_instream
			boost::filesystem::path standard_instream_flag = "-";

			/// \brief The ifstream from which file should be read
			boost::optional<std::ifstream> input_file_stream;

		public:
			path_or_istream() = default;
			explicit path_or_istream(std::istream &,
			                         const boost::filesystem::path & = "-");

			path_or_istream & set_path(const boost::filesystem::path &);
			const boost::filesystem::path & get_flag() const;
			path_or_istream & close();
			std::istream & get_istream();
		};


	} // namespace common
} // namespace cath

#endif
