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

#ifndef _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_BOOST_ADDENDA_LOG_LOG_TO_OSTREAM_GUARD_HPP
#define _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_BOOST_ADDENDA_LOG_LOG_TO_OSTREAM_GUARD_HPP

#include <iosfwd>
#include <memory>

#include <spdlog/spdlog.h>

namespace cath {

	/// \brief RAII guard for temporarily switching to logging to an ostream
	class log_to_ostream_guard final {
	  private:
		/// A shared_ptr to the spdlog logger that's being switched out of spdlog's default logger for
		/// the lifetime of the log_to_ostream_guard
		::std::shared_ptr<::spdlog::logger> logger_shptr;

	  public:
		explicit log_to_ostream_guard( std::ostream & );
		~log_to_ostream_guard() noexcept;

		/// \brief Prevent any copying
		log_to_ostream_guard( const log_to_ostream_guard & ) = delete;
		/// \brief Prevent any moving
		log_to_ostream_guard( log_to_ostream_guard && ) noexcept = delete;
		/// \brief Prevent any copying
		log_to_ostream_guard &operator=( const log_to_ostream_guard & ) = delete;
		/// \brief Prevent any moving
		log_to_ostream_guard &operator=( log_to_ostream_guard && ) noexcept = delete;

		void reset_default_logger();
	};

} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_BOOST_ADDENDA_LOG_LOG_TO_OSTREAM_GUARD_HPP
