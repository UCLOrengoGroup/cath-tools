/// \file
/// \brief The pymol_tools class definitions

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

#include "pymol_tools.h"

#include <boost/algorithm/string/join.hpp>
#include <boost/algorithm/string/replace.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/range/adaptor/transformed.hpp>

#include "exception/invalid_argument_exception.h"
#include "structure/residue_name.h"

using namespace cath;
using namespace cath::common;
using namespace std::literals::string_literals;

using boost::adaptors::transformed;
using boost::algorithm::join;
using boost::algorithm::replace_all_copy;
using boost::numeric_cast;
using std::string;

/// \brief Calculates a sensible size for some PyMOL size by calculating fitting some a formula to two other values and using that.
///
/// The formula used is y = k / (x + c)
///
/// It turns out that this means that c and k can be calculated as follows:
///
/// c =         (x_1.y_1 - x_2.y_2) / (y_2 - y_1)
///
/// k = y_1.y_2.(x_1     - x_2    ) / (y_2 - y_1)
double pymol_tools::pymol_size(const size_t &arg_x_1,
                               const double       &arg_y_1,
                               const size_t &arg_x_2,
                               const double       &arg_y_2,
                               const size_t       &arg_x
                               ) {
	// Sanity check the inputs
	using boost::math::isfinite;
	if (!isfinite(arg_y_1) || !isfinite(arg_y_2)) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Both y values must be finite numbers"));
	}
	if (arg_y_1 == arg_y_2) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("The two y values must be distinct"));
	}

	// Convert the unsigned ints to doubles
	const double x_1(numeric_cast<double>(arg_x_1));
	const double x_2(numeric_cast<double>(arg_x_2));
	const double   x(numeric_cast<double>(arg_x));

	// Calculate the constants and the new value of y
	const double   c(                     (x_1 * arg_y_1 - x_2 * arg_y_2) / (arg_y_2 - arg_y_1) );
	const double   k( arg_y_1 * arg_y_2 * (x_1           - x_2          ) / (arg_y_2 - arg_y_1) );
	const double   y(k / (x + c));

	// Return the result
	return y;
}

/// \brief Escape a residue name string for use in PyMOL
string pymol_tools::parse_residue_name_for_pymol(const residue_name &arg_residue_name ///< The residue name string to escape
                                                 ) {
	return replace_all_copy( to_string( arg_residue_name ), "-", "\\-" );
}

/// \brief Escape residue name strings for use in PyMOL
str_vec pymol_tools::parse_residue_names_for_pymol(const residue_name_vec &arg_residue_names ///< The residue name strings to escape
                                                   ) {
	str_vec new_residue_names;
	new_residue_names.reserve( arg_residue_names.size() );
	for (const residue_name &the_residue_name : arg_residue_names) {
		new_residue_names.push_back( parse_residue_name_for_pymol( the_residue_name ) );
	}
	return new_residue_names;
}


/// \brief Generate the PyMOL string to select the specified residue(s)/atom(s) in the specified object
string pymol_tools::pymol_res_seln_str(const string           &arg_name,      ///< The object in which to select residues/atoms
                                       const residue_name_vec &arg_res_names, ///< The names of the residue(s) to select
                                       const str_opt          &arg_atom       ///< (optional) The name of the atom type to select (eg "CA")
                                       ) {
	return "/"
		+ arg_name
		+ "///"
		+ join(
			arg_res_names
				| transformed( [] (const residue_name &x) {
					return parse_residue_name_for_pymol( x );
				} ),
			"+"
		)
		+ "/"
		+ arg_atom.value_or( ""s )
		+ "/";
}
