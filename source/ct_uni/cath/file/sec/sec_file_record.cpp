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

#include "cath/common/difference.hpp"
#include "cath/common/exception/invalid_argument_exception.hpp"
#include "cath/structure/protein/residue.hpp"
#include "cath/structure/protein/sec_struc.hpp"
#include "cath/structure/protein/sec_struc_planar_angles.hpp"

#include <iostream>

using namespace cath;
using namespace cath::common;
using namespace cath::file;
using namespace cath::geom;
using namespace std;

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
sec_file_record::sec_file_record(const size_t         &prm_start_residue_num, ///< The start residue of the secondary structure
                                 const size_t         &prm_stop_residue_num,  ///< The stop residue of the secondary structure
                                 const sec_struc_type &prm_sec_struc_type,    ///< The type of secondary structure: sec_struc_type::ALPHA_HELIX or sec_struc_type::BETA_STRAND
                                 coord                 prm_midpoint,          ///< The coordinates of the secondary structure's midpoint
                                 coord                 prm_unit_dirn          ///< A unit vector along the secondary structure
                                 ) : start_residue_num { prm_start_residue_num      },
                                     stop_residue_num  { prm_stop_residue_num       },
                                     type              { prm_sec_struc_type         },
                                     midpoint          { std::move( prm_midpoint  ) },
                                     unit_dirn         { std::move( prm_unit_dirn ) } {
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

/// \brief Return whether the two specified sec_file_records are identical
///
/// \relates sec_file_record
bool cath::file::operator==(const sec_file_record &prm_lhs, ///< The first  sec_file_record to compare
                            const sec_file_record &prm_rhs  ///< The second sec_file_record to compare
                            ) {
	return (
		prm_lhs.get_start_residue_num() == prm_rhs.get_start_residue_num()
		&&
		prm_lhs.get_stop_residue_num()  == prm_rhs.get_stop_residue_num()
		&&
		prm_lhs.get_type()              == prm_rhs.get_type()
		&&
		prm_lhs.get_midpoint()          == prm_rhs.get_midpoint()
		&&
		prm_lhs.get_unit_dirn()         == prm_rhs.get_unit_dirn()
	);
}

/// \brief Generate a string describing the specified sec_file_record
///
/// \relates sec_file_record
string cath::file::to_string(const sec_file_record &prm_sec_file_record ///< The sec_file_record to describe
                             ) {
	return
		  "sec_file_record["
		+ std::to_string( prm_sec_file_record.get_start_residue_num() )
		+ ", "
		+ std::to_string( prm_sec_file_record.get_stop_residue_num()  )
		+ ", "
		+ to_string     ( prm_sec_file_record.get_type()              )
		+ ", "
		+ to_string     ( prm_sec_file_record.get_midpoint()          )
		+ ", "
		+ to_string     ( prm_sec_file_record.get_unit_dirn()         )
		+ "]";
}

/// \brief Insert a description of the specified sec_file_record into the specified ostream
///
/// \relates sec_file_record
ostream & cath::file::operator<<(ostream               &prm_os,             ///< The ostream into which the description should be inserted
                                 const sec_file_record &prm_sec_file_record ///< The sec_file_record to describe
                                 ) {
	prm_os << to_string( prm_sec_file_record );
	return prm_os;
}

/// \brief Convert a sec_file_record to a sec_struc
///
/// \relates sec_file_record
sec_struc cath::file::make_sec_struc(const sec_file_record &prm_sec_file_record ///< The sec_file_record to be converted
                                     ) {
	return sec_struc(
		prm_sec_file_record.get_start_residue_num(),
		prm_sec_file_record.get_stop_residue_num(),
		prm_sec_file_record.get_type(),
		prm_sec_file_record.get_midpoint(),
		prm_sec_file_record.get_unit_dirn()
	);
}
