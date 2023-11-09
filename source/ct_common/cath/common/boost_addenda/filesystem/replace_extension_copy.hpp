/// \file
/// \brief The replace_extension_copy header

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

#ifndef CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_BOOST_ADDENDA_FILESYSTEM_REPLACE_EXTENSION_COPY_HPP
#define CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_BOOST_ADDENDA_FILESYSTEM_REPLACE_EXTENSION_COPY_HPP

#include <filesystem>

namespace cath::common {

	/// \brief Return a copy of the specified path in which the extension has been
	///        replaced with the specified replacement
	inline ::std::filesystem::path replace_extension_copy(::std::filesystem::path        prm_file,
	                                                      const ::std::filesystem::path &prm_replacement = ::std::filesystem::path()
	                                                      ) {
		prm_file.replace_extension( prm_replacement );
		return prm_file;
	}

} // namespace cath::common

#endif // CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_BOOST_ADDENDA_FILESYSTEM_REPLACE_EXTENSION_COPY_HPP
