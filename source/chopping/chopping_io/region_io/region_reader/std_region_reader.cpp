/// \file
/// \brief The std_region_reader class definitions

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

#include "std_region_reader.h"

#include "chopping/region/region.h"
#include "exception/not_implemented_exception.h" // ***** TEMPORARY *****

#include <iostream> // ***** TEMPORARY *****

using namespace cath::chop;
using namespace cath::common;
using namespace std;

/// \brief TODOCUMENT
region std_region_reader::do_read_region(const string &arg_region_string ///< TODOCUMENT
                                         ) const {
	cerr << "arg_region_string is " << arg_region_string << endl;

	BOOST_THROW_EXCEPTION(not_implemented_exception("std_region_reader::do_read_region()"));

	return region( 0, 0 );
}

/// \brief Ctor for std_region_reader
std_region_reader::std_region_reader(const std_region_io_spec &arg_std_region_io_spec ///< TODOCUMENT
                                     ) : region_io_spec( arg_std_region_io_spec ) {

}