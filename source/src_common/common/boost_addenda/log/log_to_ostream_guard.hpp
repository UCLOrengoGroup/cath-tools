/// \file
/// \brief The log_to_ostream_guard class header

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

#ifndef _CATH_TOOLS_SOURCE_SRC_COMMON_COMMON_BOOST_ADDENDA_LOG_LOG_TO_OSTREAM_GUARD_H
#define _CATH_TOOLS_SOURCE_SRC_COMMON_COMMON_BOOST_ADDENDA_LOG_LOG_TO_OSTREAM_GUARD_H

#include <boost/log/utility/setup/console.hpp>
using sink_t     = boost::log::sinks::synchronous_sink<boost::log::sinks::basic_text_ostream_backend<char>>;
using sink_bsptr = boost::shared_ptr<sink_t>;

#include <iosfwd>

namespace cath {

	/// \brief TODOCUMENT
	class log_to_ostream_guard final {
	private:
		/// \brief A shared_ptr to a Boost Log sink that is added to the core for the lifetime of the log_to_ostream_guard
		sink_bsptr boost_log_sink_bsptr;

	public:
		explicit log_to_ostream_guard(std::ostream &);
		~log_to_ostream_guard() noexcept;

		/// \brief Specify that the copy-ctor shouldn't be used
		log_to_ostream_guard(const log_to_ostream_guard &) = delete;
		/// \brief Specify that the copy-assign shouldn't be used
		log_to_ostream_guard & operator=(const log_to_ostream_guard &) = delete;

		void remove_log_sink();
	};

} // namespace cath

#endif
