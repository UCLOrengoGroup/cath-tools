/// \file
/// \brief The data_dirs_options_block class definitions

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

#include "data_file.hpp"

#include <boost/lexical_cast.hpp>
#include <boost/range/adaptor/transformed.hpp>
#include <boost/range/algorithm/max_element.hpp>

#include "common/exception/invalid_argument_exception.hpp"

using namespace cath::common;
using namespace cath::file;
using namespace std;

using boost::adaptors::transformed;
using boost::range::max_element;

/// \brief TODOCUMENT
///
/// \relates data_file
string cath::file::to_string(const data_file &prm_data_file ///< The data_file to output
                             ) {
	switch ( prm_data_file ) {
		case ( data_file::PDB  ) : { return "data_file::PDB"  ; }
		case ( data_file::DSSP ) : { return "data_file::DSSP" ; }
		case ( data_file::WOLF ) : { return "data_file::WOLF" ; }
		case ( data_file::SEC  ) : { return "data_file::SEC"  ; }
	}
	BOOST_THROW_EXCEPTION(invalid_argument_exception("Value of data_file not recognised whilst inserting into an ostream"));
}

/// \brief Simple insertion operator for data_file
///
/// \relates data_file
ostream & cath::file::operator<<(ostream         &prm_os,       ///< The ostream to which to output the data_file
                                 const data_file &prm_data_file ///< The data_file to output
                                 ) {
	prm_os << to_string( prm_data_file );
	return prm_os;
}

/// \brief The length of the string representing the specified data_file
///
/// \relates data_file
size_t cath::file::str_length_of_data_file(const data_file &prm_data_file ///< The data_file whose string length should be returned
                                           ) {
	return to_string( prm_data_file ).length();
}

/// \brief The maximum length of the string representation of all data_file values
///
/// \relates data_file
size_t cath::file::max_data_file_str_length() {
	return *max_element( detail::all_data_file_types | transformed( &str_length_of_data_file ) );
}
