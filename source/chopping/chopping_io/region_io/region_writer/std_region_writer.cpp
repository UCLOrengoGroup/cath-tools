/// \file
/// \brief The std_region_writer class definitions

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

#include "std_region_writer.hpp"

#include "exception/not_implemented_exception.hpp" // ***** TEMPORARY *****

#include <iostream> // ***** TEMPORARY *****

using namespace cath::chop;
using namespace cath::common;

using std::cerr; // ***** TEMPORARY *****
using std::string;

/// \brief Ctor for std_region_writer
std_region_writer::std_region_writer(const std_region_io_spec &arg_region_io_spec ///< TODOCUMENT
                                     ) : region_io_spec( arg_region_io_spec ) {
}

/// \brief TODOCUMENT
string std_region_writer::do_write_region(const region &/*arg_region*/ ///< TODOCUMENT
                                          ) const {
	cerr << "Writing region, but not yet implemented"<< "\n";

	BOOST_THROW_EXCEPTION(not_implemented_exception("std_region_writer::do_write_region()"));

	return "";
}
