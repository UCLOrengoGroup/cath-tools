/// \file
/// \brief The chain_label class definitions

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

#include "chain_label.hpp"

#include <boost/algorithm/string/classification.hpp>

#include "exception/invalid_argument_exception.hpp"

using namespace boost::algorithm;
using namespace cath;
using namespace cath::common;
using namespace std;

/// \brief TODOCUMENT
void chain_label::sanity_check() const {
	if ( ! is_print()( chain_char ) ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Invalid, non-printing chain code in cath::write_pdb_file_entry()"));
	}
}

/// \brief Ctor for chain_label
chain_label::chain_label(const char &arg_chain_char ///< TODOCUMENT
                         ) : chain_char( arg_chain_char ) {
	sanity_check();
}

/// \brief TODOCUMENT
string chain_label::to_string() const {
	return { chain_char };
}

/// \brief TODOCUMENT
bool cath::operator==(const chain_label &arg_chain_label_a, ///< TODOCUMENT
                      const chain_label &arg_chain_label_b  ///< TODOCUMENT
                      ) {
	return ( arg_chain_label_a.get_char() == arg_chain_label_b.get_char() );
}

/// \brief TODOCUMENT
///
/// \relates chain_label
ostream & cath::operator<<(ostream           &arg_os,         ///< TODOCUMENT
                           const chain_label &arg_chain_label ///< TODOCUMENT
				           ) {
	arg_os << arg_chain_label.to_string();
	return arg_os;
}
