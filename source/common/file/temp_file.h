/// \file
/// \brief The temp_file class header

/// \copyright
/// Tony Lewis's Common C++ Library Code (here imported into the CATH Binaries project and then tweaked, eg namespaced in cath)
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

#ifndef TEMP_FILE_H_INCLUDED
#define TEMP_FILE_H_INCLUDED

//#include <boost/filesystem/path.hpp>
#include <boost/optional.hpp>

#include "common/path_type_aliases.h"

#include <string>

namespace cath {

	/// \brief A simple wrapper for creating a temporary file and then automatically cleaning it up on destruction
	///
	/// This can be constructed with an empty string, in which case has_filename() will return false,
	/// get_filename() will return an empty path and the destructor will do nothing.
	class temp_file final {
	private:
		opt_path filename;

		static boost::filesystem::path temp_filename_of_basename_pattern(const std::string &);

	public:
		explicit temp_file(const std::string &);
		~temp_file() noexcept;

		/// \brief Specify that the copy-ctor shouldn't be used
		///
		/// \todo Consider implementing a move-ctor that ensures the source filename gets wiped
		temp_file(const temp_file &) = delete;
		/// \brief Specify that the copy-assign shouldn't be used
		///
		/// \todo Consider implementing a move-assign that ensures the source filename gets wiped
		temp_file & operator=(const temp_file &) = delete;

		const opt_path & get_opt_filename() const;
	};

	bool has_filename(const temp_file &);
	boost::filesystem::path get_filename(const temp_file &);
}

#endif