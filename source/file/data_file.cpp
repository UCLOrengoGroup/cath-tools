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

#include "data_file.h"

//#include <boost/algorithm/string/case_conv.hpp>
//#include <boost/algorithm/string/classification.hpp>
//#include <boost/algorithm/string/join.hpp>
//#include <boost/algorithm/string/replace.hpp>
//#include <boost/algorithm/string/split.hpp>
#include <boost/lexical_cast.hpp>
//#include <boost/optional.hpp>
#include <boost/range/adaptor/transformed.hpp>
#include <boost/range/algorithm/max_element.hpp>

//#include "common/algorithm/contains.h"
//#include "common/boost_addenda/string_algorithm/split_build.h"
//#include "common/c++14/make_unique.h"
//#include "common/clone/make_uptr_clone.h"
#include "exception/invalid_argument_exception.h"
//#include "exception/runtime_error_exception.h"

//using namespace boost::algorithm;
//using namespace boost::filesystem;
//using namespace boost::program_options;
//using namespace cath;
using namespace cath::common;
using namespace cath::file;
//using namespace cath::opts::detail;
//using namespace cath::opts;
using namespace std;

using boost::adaptors::transformed;
//using boost::algorithm::is_any_of;
//using boost::algorithm::join;
//using boost::algorithm::replace_all_copy;
//using boost::algorithm::token_compress_on;
using boost::lexical_cast;
//using boost::none;
using boost::range::max_element;

/// \brief TODOCUMENT
///
/// \relates data_file
string cath::file::to_string(const data_file &arg_data_file ///< The data_file to output
                             ) {
	switch (arg_data_file) {
		case ( data_file::PDB  ) : { return "data_file::PDB"  ; }
		case ( data_file::DSSP ) : { return "data_file::DSSP" ; }
		case ( data_file::WOLF ) : { return "data_file::WOLF" ; }
		case ( data_file::SEC  ) : { return "data_file::SEC"  ; }
		default : {
			BOOST_THROW_EXCEPTION(invalid_argument_exception("Value of data_file not recognised whilst inserting into an ostream"));
			return ""; // Superfluous, post-throw return statement to appease Eclipse's syntax highlighter
		}
	}
}

/// \brief Simple insertion operator for data_file
///
/// \relates data_file
ostream & cath::file::operator<<(ostream         &arg_os,       ///< The ostream to which to output the data_file
                                 const data_file &arg_data_file ///< The data_file to output
                                 ) {
	switch (arg_data_file) {
		case ( data_file::PDB  ) : { arg_os << "data_file::PDB"  ; break; }
		case ( data_file::DSSP ) : { arg_os << "data_file::DSSP" ; break; }
		case ( data_file::WOLF ) : { arg_os << "data_file::WOLF" ; break; }
		case ( data_file::SEC  ) : { arg_os << "data_file::SEC"  ; break; }
		default : {
			BOOST_THROW_EXCEPTION(invalid_argument_exception("Value of data_file not recognised whilst inserting into an ostream"));
			break; // Superfluous, post-throw break statement to appease Eclipse's syntax highlighter
		}
	}
	return arg_os;
}

/// \brief TODOCUMENT
///
/// \relates data_file
size_t cath::file::str_length_of_data_file(const data_file &arg_data_file ///< TODOCUMENT
                                           ) {
	return to_string( arg_data_file ).length();
}

/// \brief TODOCUMENT
///
/// \relates data_file
size_t cath::file::max_data_file_str_length() {
	return *max_element( detail::all_data_file_types | transformed( &str_length_of_data_file ) );
}
