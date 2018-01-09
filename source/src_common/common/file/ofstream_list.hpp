/// \file
/// \brief The open_ifstream / open_ofstream header

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

#ifndef _CATH_TOOLS_SOURCE_SRC_COMMON_COMMON_FILE_OFSTREAM_LIST_HPP
#define _CATH_TOOLS_SOURCE_SRC_COMMON_COMMON_FILE_OFSTREAM_LIST_HPP

#include <boost/filesystem.hpp>
#include <boost/optional.hpp>

#include "common/path_type_aliases.hpp"
#include "common/type_aliases.hpp"

#include <fstream>
#include <functional>

namespace cath {
	namespace common {

		/// \brief Type alias for a vector of reference_wrapper of ostream
		using ostream_ref_vec = std::vector<ostream_ref>;

		/// \brief A list of ostreams, with support for populating from paths, which get automatically opened,
		///        and a special flag which indicates output to the ostream optionally specified on construction
		///
		/// Note: this has substantial overlap with path_or_istream and could perhaps share a common implementation
		class ofstream_list {
		private:
			/// \brief An optional special ostream to which output can be sent (usually stdout)
			ostream_ref_opt standard_outstream;

			/// \brief A flag that can be used when passing a path to indicate output should be sent to the standard_outstream
			boost::filesystem::path standard_outstream_flag = "-";

			/// \brief The standard list of ofstreams to which output should be sent
			std::deque<std::ofstream> ofstreams;

		public:
			ofstream_list() = default;

			explicit ofstream_list(std::ostream &,
			                       const boost::filesystem::path & = "-");

			ostream_ref_vec open_ofstreams(const path_vec &);
			const boost::filesystem::path & get_flag() const;
			void close_all();
		};

		ostream_ref open_ofstream(ofstream_list &,
		                          const boost::filesystem::path &);

	} // namespace common
} // namespace cath

#endif
