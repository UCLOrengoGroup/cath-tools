/// \file
/// \brief The out_of_range_exception class header

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

#ifndef OUT_OF_RANGE_EXCEPTION_H_INCLUDED
#define OUT_OF_RANGE_EXCEPTION_H_INCLUDED

#include <boost/exception/all.hpp>

namespace cath {
	namespace common {

		/// \brief TODOCUMENT
		class out_of_range_exception : public boost::exception,
		                               public std::out_of_range {
		public:
			explicit out_of_range_exception(const std::string &);
			virtual ~out_of_range_exception() noexcept = default;
		};

	}
}

#endif