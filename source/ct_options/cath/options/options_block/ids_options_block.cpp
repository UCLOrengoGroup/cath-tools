/// \file
/// \brief The ids_options_block class definitions

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

#include "ids_options_block.hpp"

#include <boost/lexical_cast.hpp>
#include <boost/optional.hpp>
#include <boost/program_options.hpp>

#include "cath/common/clone/make_uptr_clone.hpp"
#include "cath/common/exception/invalid_argument_exception.hpp"
#include "cath/common/exception/out_of_range_exception.hpp"

#include <iostream>

using namespace ::cath;
using namespace ::cath::common;
using namespace ::cath::opts;
using namespace ::std;

using ::boost::none;
using ::boost::program_options::options_description;
using ::boost::program_options::value;
using ::boost::program_options::variables_map;

const string ids_options_block::PO_ID( "id" ); ///< The option name for the id option

/// \brief A standard do_clone() method to act as a virtual copy-ctor
///
/// This is a concrete definition of a virtual method that's pure in options_block
unique_ptr<options_block> ids_options_block::do_clone() const {
	return { make_uptr_clone( *this ) };
}

/// \brief Provide a name for the options block, as used in the options description text
///
/// This is a concrete definition of a virtual method that's pure in options_block
string ids_options_block::do_get_block_name() const {
	return "ID options";
}

/// \brief Add the block's non-hidden options to the provided options_description
///
/// This is a concrete definition of a virtual method that's pure in options_block
void ids_options_block::do_add_visible_options_to_description(options_description &prm_desc,           ///< The options_description to which the options are added
                                                              const size_t        &/*prm_line_length*/ ///< The line length to be used when outputting the description (not very clearly documented in Boost)
                                                              ) {
	prm_desc.add_options()
		( PO_ID.c_str(), value<str_vec>( &ids ), "Structure ids" );
}

///// \brief Add any hidden options to the provided options_description
//void ids_options_block::do_add_hidden_options_to_description(options_description &/*prm_desc*/ ///< The options_description to which the options are added
//                                                                  ) {
////	prm_desc.add_options()
////		( PO_ID.c_str(), value<str_vec>( &ids ), "Structure ids" );
//}

/// \brief Identify any conflicts that make the currently stored options invalid
///
/// \returns A string describing the conflict in the options or an empty string if there's none
str_opt ids_options_block::do_invalid_string(const variables_map &/*prm_variables_map*/ ///< The variables map, which options_blocks can use to determine which options were specified, defaulted etc
                                             ) const {
	return none;
}

/// \brief Return all options names for this block
str_vec ids_options_block::do_get_all_options_names() const {
	return {
		ids_options_block::PO_ID,
	};
}

// /// \brief TODOCUMENT
// size_t ids_options_block::num_ids_specified() const {
// 	return ids.size();
// }

// /// \brief TODOCUMENT
// string ids_options_block::get_id_of_index(const size_t &prm_index ///< TODOCUMENT
//                                                ) const {
// 	return ids[ prm_index ];
// }

/// \brief TODOCUMENT
const str_vec & ids_options_block::get_ids() const {
	return ids;
}

// /// \brief TODOCUMENT
// str_vec cath::opts::get_all_ids(const ids_options_block &prm_block ///< TODOCUMENT
//                                 ) {
// 	const size_t num_ids = prm_block.num_ids_specified();

// 	str_vec all_ids;
// 	all_ids.reserve( num_ids );
// 	for (const size_t &id_ctr : indices( num_ids ) ) {
// 		all_ids.push_back( prm_block.get_id_of_index( id_ctr ) );
// 	}

// 	return all_ids;
// }

/// \brief TODOCUMENT
bool cath::opts::ids_specified(const ids_options_block &prm_block ///< TODOCUMENT
                               ) {
	return ( ! prm_block.get_ids().empty() );
}

/// \brief Getter for the name of the first protein structure to be compared
string cath::opts::get_id_a(const ids_options_block &prm_block ///< TODOCUMENT
                            ) {
	if ( prm_block.get_ids().size() != 2 ) {
		BOOST_THROW_EXCEPTION(out_of_range_exception("Cannot get protein name a because two names have not (yet) been specified"));
	}
	return prm_block.get_ids()[ 0 ];
}

/// \brief Getter for the name of the second protein structure to be compared
string cath::opts::get_id_b(const ids_options_block &prm_block ///< TODOCUMENT
                            ) {
	if ( prm_block.get_ids().size() != 2 ) {
		BOOST_THROW_EXCEPTION(out_of_range_exception("Cannot get protein name b because two names have not (yet) been specified"));
	}
	return prm_block.get_ids()[ 1 ];
}

