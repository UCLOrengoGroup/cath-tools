/// \file
/// \brief The not_implemented_exception class header

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

#ifndef _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_EXCEPTION_NOT_IMPLEMENTED_EXCEPTION_HPP
#define _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_EXCEPTION_NOT_IMPLEMENTED_EXCEPTION_HPP

#include <boost/exception/all.hpp>

namespace cath::common {

	/// \brief TODOCUMENT
	class not_implemented_exception : public boost::exception,
	                                  public std::logic_error {
	public:
		explicit not_implemented_exception(const std::string &);
	};

} // namespace cath::common

#endif // _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_EXCEPTION_NOT_IMPLEMENTED_EXCEPTION_HPP
