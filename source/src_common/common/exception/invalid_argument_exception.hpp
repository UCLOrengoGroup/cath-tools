/// \file
/// \brief The invalid_argument_exception class header

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

#ifndef _CATH_TOOLS_SOURCE_SRC_COMMON_COMMON_EXCEPTION_INVALID_ARGUMENT_EXCEPTION_H
#define _CATH_TOOLS_SOURCE_SRC_COMMON_COMMON_EXCEPTION_INVALID_ARGUMENT_EXCEPTION_H

#include <boost/exception/all.hpp>

namespace cath {
	namespace common {

		/// \brief TODOCUMENT
		class invalid_argument_exception : public boost::exception,
		                                   public std::invalid_argument {
		public:
			explicit invalid_argument_exception(const std::string &);
		};

	} // namespace common
} // namespace cath

#endif
