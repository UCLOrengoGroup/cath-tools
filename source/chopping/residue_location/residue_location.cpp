/// \file
/// \brief The residue_location class definitions

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

#include "residue_location.hpp"

#include "common/exception/invalid_argument_exception.hpp"
#include "structure/structure_type_aliases.hpp"

using namespace cath;
using namespace cath::common;
using namespace cath::chop;

/// \brief Ctor for residue_location
residue_location::residue_location(const residue_name &arg_residue_name ///< TODOCUMENT
                                   ) : the_residue_name( arg_residue_name )  {

}

/// \brief Ctor for residue_location
residue_location::residue_location(const residue_name &arg_residue_name, ///< TODOCUMENT
                                   const size_t       &arg_residue_index ///< TODOCUMENT
                                   ) : the_residue_name( arg_residue_name  ),
                                       residue_index   ( arg_residue_index ) {

}

/// \brief Ctor for residue_location
residue_location::residue_location(const size_t &arg_residue_index ///< TODOCUMENT
                                   ) : residue_index( arg_residue_index ) {

}

/// \brief TODOCUMENT
const residue_name_opt & residue_location::get_opt_residue_name() const {
	return the_residue_name;
}

/// \brief TODOCUMENT
const size_opt & residue_location::get_opt_residue_index() const {
	return residue_index;
}

/// \brief TODOCUMENT
///
/// \relates residue_location
bool cath::chop::has_residue_name(const residue_location &arg_residue_location ///< TODOCUMENT
                                  ) {
	return static_cast<bool>( arg_residue_location.get_opt_residue_name() );
}

/// \brief TODOCUMENT
///
/// \relates residue_location
bool cath::chop::has_residue_index(const residue_location &arg_residue_location ///< TODOCUMENT
                                   ) {
	return static_cast<bool>( arg_residue_location.get_opt_residue_index() );
}

/// \brief TODOCUMENT
///
/// \relates residue_location
residue_name cath::chop::get_residue_name(const residue_location &arg_residue_location ///< TODOCUMENT
                                          ) {
	return *arg_residue_location.get_opt_residue_name();
}

/// \brief TODOCUMENT
///
/// \relates residue_location
size_t cath::chop::get_residue_index(const residue_location &arg_residue_location ///< TODOCUMENT
                                     ) {
	return *arg_residue_location.get_opt_residue_index();
}

/// \brief TODOCUMENT
///
/// \relates residue_location
residue_locating cath::chop::get_residue_locating(const residue_location &arg_residue_location ///< TODOCUMENT
                                                  ) {
	return make_residue_locating_of_has_name_and_has_index(
		has_residue_name ( arg_residue_location ),
		has_residue_index( arg_residue_location )
	);
}

/// \brief Non-member equality operator for residue_location
///
/// \relates residue_location
bool cath::chop::operator==(const residue_location &arg_lhs, ///< The first  residue_location to compare
                            const residue_location &arg_rhs  ///< The second residue_location to compare
                            ) {
	return (
		arg_lhs.get_opt_residue_name()  == arg_rhs.get_opt_residue_name()
		&&
		arg_lhs.get_opt_residue_index() == arg_rhs.get_opt_residue_index()
	);
}

/// \brief TODOCUMENT
///
/// \relates residue_location
///
/// \pre Both residues must have index information
///      (ie has_residue_name( arg_residue_location_a) and has_residue_name( arg_residue_location_b) )
bool cath::chop::operator<(const residue_location &arg_residue_location_a, ///< TODOCUMENT
                           const residue_location &arg_residue_location_b  ///< TODOCUMENT
                           ) {
	// Check pre-condition that the both residue_locations have indices
	const bool has_residue_index_a = has_residue_index( arg_residue_location_a );
	const bool has_residue_index_b = has_residue_index( arg_residue_location_b );
	if ( ! has_residue_index_a && ! has_residue_index_b ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Unable to compare two residue locations because neither has an associated index"));
	}
	else if ( ! has_residue_index_a ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Unable to compare two residue locations because the first has no associated index"));
	}
	else if ( ! has_residue_index_b ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Unable to compare two residue locations because the second has no associated index"));
	}

	// Both have indices so return the result of comparing them
	const size_t residue_index_a = get_residue_index( arg_residue_location_a );
	const size_t residue_index_b = get_residue_index( arg_residue_location_b );
	return ( residue_index_a < residue_index_b );
}
