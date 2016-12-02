/// \file
/// \brief The sec_file_record class definitions

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

#include "sec_file_record.hpp"

#include <boost/math/special_functions/fpclassify.hpp>

#include "common/difference.hpp"
#include "exception/invalid_argument_exception.hpp"
#include "structure/protein/residue.hpp"
#include "structure/protein/sec_struc.hpp"
#include "structure/protein/sec_struc_planar_angles.hpp"

#include <iostream>

using namespace boost::math;
using namespace cath;
using namespace cath::common;
using namespace cath::file;
using namespace cath::geom;
using namespace std;

using boost::math::isfinite;

/// \brief Check that the unit direction vector has a length vaguely close to one
void sec_file_record::check_unit_dirn_length() {
	if ( difference( 1.0, length( get_unit_dirn() ) ) > 0.1 ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("sec_file_record's unit direction "));
	}
}

/// \brief Check that the secondary structure type is valid
void sec_file_record::check_sec_struc_type() const {
	if (type != sec_struc_type::ALPHA_HELIX && type != sec_struc_type::BETA_STRAND) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("sec_file_record's type cannot be sec_struc_type::COIL"));
	}
}

/// \brief Ctor that directly populates all the fields with arguments
sec_file_record::sec_file_record(const size_t         &arg_start_residue_num, ///< The start residue of the secondary structure
                                 const size_t         &arg_stop_residue_num,  ///< The stop residue of the secondary structure
                                 const sec_struc_type &arg_sec_struc_type,    ///< The type of secondary structure: sec_struc_type::ALPHA_HELIX or sec_struc_type::BETA_STRAND
                                 const coord          &arg_midpoint,          ///< The coordinates of the secondary structure's midpoint
                                 const coord          &arg_unit_dirn          ///< A unit vector along the secondary structure
                                 ) : start_residue_num( arg_start_residue_num ),
                                     stop_residue_num ( arg_stop_residue_num  ),
                                     type             ( arg_sec_struc_type    ),
                                     midpoint         ( arg_midpoint          ),
                                     unit_dirn        ( arg_unit_dirn         ) {
	check_sec_struc_type();
	unit_dirn = normalise_copy( unit_dirn );
}

/// \brief Getter for the start residue number
size_t sec_file_record::get_start_residue_num() const {
	return start_residue_num;
}

/// \brief Getter for the stop residue number
size_t sec_file_record::get_stop_residue_num() const {
	return stop_residue_num;
}

/// \brief Getter for the secondary structure type
sec_struc_type sec_file_record::get_type() const {
	return type;
}

/// \brief Getter for the secondary structure's midpoint
coord sec_file_record::get_midpoint() const {
	return midpoint;
}

/// \brief Getter for the unit vector along the secondary structure
coord sec_file_record::get_unit_dirn() const {
	return unit_dirn;
}

/// \brief Convert a sec_file_record to a sec_struc
///
/// \relates sec_file_record
sec_struc cath::file::make_sec_struc(const sec_file_record &arg_sec_file_record ///< The sec_file_record to be converted
                                     ) {
	return sec_struc(
		arg_sec_file_record.get_start_residue_num(),
		arg_sec_file_record.get_stop_residue_num(),
		arg_sec_file_record.get_type(),
		arg_sec_file_record.get_midpoint(),
		arg_sec_file_record.get_unit_dirn()
	);
}
