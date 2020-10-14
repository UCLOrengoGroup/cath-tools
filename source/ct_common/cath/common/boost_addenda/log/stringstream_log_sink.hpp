/// \file
/// \brief The stringstream_log_sink class header

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

#ifndef _CATH_TOOLS_SOURCE_SRC_COMMON_COMMON_BOOST_ADDENDA_LOG_STRINGSTREAM_LOG_SINK_HPP
#define _CATH_TOOLS_SOURCE_SRC_COMMON_COMMON_BOOST_ADDENDA_LOG_STRINGSTREAM_LOG_SINK_HPP

#include "cath/common/boost_addenda/log/log_to_ostream_guard.hpp"

#include <sstream>

namespace cath {

	/// \brief Grab all boost logging to an internal stringstream
	///
	/// This is useful for checking for the presence or absence of log
	/// messages in tests. Or for simply silencing warnings during tests.
	class stringstream_log_sink final {
	private:
		/// \brief The stringstream to which all Boost logging should be redirected
		::std::ostringstream out_ss;

		/// \brief The log_to_ostream_guard to grab the Boost logging
		log_to_ostream_guard out_ss_guard{ out_ss };

	public:
		stringstream_log_sink() = default;

		/// \brief Keep things simple by having a single instance of the guard
		stringstream_log_sink(const stringstream_log_sink &) = delete;
		/// \brief Keep things simple by having a single instance of the guard
		stringstream_log_sink & operator=(const stringstream_log_sink &) = delete;

		::std::ostringstream & stringstream();
		const ::std::ostringstream & stringstream() const;
		::std::string str() const;
		bool str_is_empty() const;
	};

} // namespace cath

#endif
