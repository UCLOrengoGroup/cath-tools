/// \file
/// \brief The residue_name_align_map class definitions

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

#include "residue_name_align_map.hpp"

#include <boost/lexical_cast.hpp>
#include <boost/numeric/conversion/cast.hpp>

#include "common/algorithm/contains.hpp"
#include "common/type_aliases.hpp"
#include "exception/invalid_argument_exception.hpp"
#include "structure/residue_name.hpp"

using namespace cath;
using namespace cath::align::detail;
using namespace cath::common;
using namespace std;

using boost::lexical_cast;
using boost::numeric_cast;

/// \brief Ctor for residue_name_align_map
residue_name_align_map::residue_name_align_map(const str_vec &arg_residue_name_strings ///< The list of residue names to be indexed
                                               ) {
	for (const string &residue_name_string : arg_residue_name_strings) {
		if ( contains( index_of_residue_name, lexical_cast<string>( residue_name_string ) ) ) {
			BOOST_THROW_EXCEPTION(invalid_argument_exception(
				"Whilst residue-name aligning, arg_residue_names contains duplicate entry : \""
				+ residue_name_string
				+ "\""
			));
		}
		const auto index = index_of_residue_name.size();
		index_of_residue_name.insert( make_pair( residue_name_string, index ) );
	}
}

/// \brief Grab whether or not the residue_name_align_map contains a residue specified by its name
bool residue_name_align_map::contains_residue_name_string(const string &arg_residue_name ///< The name of the residue to queried
                                                          ) const {
	return ( contains( index_of_residue_name, arg_residue_name ) );
}

/// \brief Grab the index of a residue specified by its name
size_t residue_name_align_map::get_index_of_residue_name_string(const string &arg_residue_name ///< The name of the residue to queried
                                                                ) const {
	if ( ! contains (index_of_residue_name, arg_residue_name ) ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Whilst residue-name aligning, residue name \"" + arg_residue_name + "\" cannot be found"));
	}
	return index_of_residue_name.find(arg_residue_name)->second;
}

/// \brief Get a list of all the residue names (ie the keys)
str_vec residue_name_align_map::get_residue_name_strings() const {
	str_vec residue_names;
	residue_names.reserve(index_of_residue_name.size());
	for (const str_size_pair &residue_name_index_pair : index_of_residue_name) {
		residue_names.push_back(residue_name_index_pair.first);
	}
	return residue_names;
}

/// \brief TODOCUMENT
///
/// \relates residue_name_align_map
residue_name_align_map cath::align::detail::make_residue_name_align_map(const residue_name_vec &arg_residue_names ///< TODOCUMENT
                                                                        ) {
	str_vec residue_name_strings;
	residue_name_strings.reserve( arg_residue_names.size() );
	for (const residue_name &the_residue_name : arg_residue_names) {
		residue_name_strings.push_back( lexical_cast<string>( the_residue_name ) );
	}
	return residue_name_align_map( residue_name_strings );
}

/// \brief TODOCUMENT
///
/// \relates residue_name_align_map
bool cath::align::detail::contains_residue_name(const residue_name_align_map &arg_residue_name_align_map, ///< TODOCUMENT
                                                const residue_name           &arg_residue_name            ///< TODOCUMENT
                                                ) {
	return arg_residue_name_align_map.contains_residue_name_string(
		lexical_cast<string>( arg_residue_name )
	);
}

/// \brief TODOCUMENT
///
/// \relates residue_name_align_map
size_t cath::align::detail::get_index_of_residue_name(const residue_name_align_map &arg_residue_name_align_map, ///< TODOCUMENT
                                                      const residue_name           &arg_residue_name            ///< TODOCUMENT
                                                      ) {
	return arg_residue_name_align_map.get_index_of_residue_name_string(
		lexical_cast<string>( arg_residue_name )
	);
}

/// \brief Insertion operator to summarise a arg_residue_name_align_map
///
/// \relates residue_name_align_map
ostream & cath::align::detail::operator<<(ostream                      &arg_os,                    ///< ostream to which to dump the summary
                                          const residue_name_align_map &arg_residue_name_align_map ///< The residue_name_align_map to summarise
                                          ) {
	arg_os << "residue_name_align_map[";
	const str_vec residue_names = arg_residue_name_align_map.get_residue_name_strings();
	for (const string &residue_name_string : residue_names) {
		arg_os << " ";
		arg_os << residue_name_string;
		arg_os << " -> ";
		arg_os << arg_residue_name_align_map.get_index_of_residue_name_string( residue_name_string );
	}
	arg_os << "]";
	return arg_os;
}

