/// \file
/// \brief The sec_file class definitions

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

#include "sec_file.hpp"

#include <boost/range/irange.hpp>

#include "common/algorithm/transform_build.hpp"
#include "common/boost_addenda/range/indices.hpp"
#include "common/cpp14/cbegin_cend.hpp"
#include "common/exception/invalid_argument_exception.hpp"
#include "common/size_t_literal.hpp"
#include "file/sec/sec_file_record.hpp"
#include "structure/protein/sec_struc.hpp"
#include "structure/protein/sec_struc_planar_angles.hpp"

using namespace cath;
using namespace cath::common;
using namespace cath::file;

using boost::irange;
using std::max;
using std::min;

/// \brief Ctor to populate the sec_file_records and planar_angles_lists
sec_file::sec_file(sec_file_record_vec             arg_sec_file_records,   ///< The list of sec_file_records
                   sec_struc_planar_angles_vec_vec arg_inter_planar_angles ///< The list of planar_angle_lists
                   ) : records             { std::move( arg_sec_file_records    ) },
                       inter_planar_angles { std::move( arg_inter_planar_angles ) } {
	if ( records.empty() ) {
		if ( ! inter_planar_angles.empty() ) {
			BOOST_THROW_EXCEPTION(invalid_argument_exception("There are inter_planar_angles despite there being no sec file records"));
		}
	}
	else if ( inter_planar_angles.size() + 1 != records.size() ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception(
			"Number of inter_planar_angles ("
			+ std::to_string( inter_planar_angles.size() )
			+ ") does not match one less than the number of sec file records ("
			+ std::to_string( records.size() )
			+ ")"
		));
	}
	for (const size_t &planar_ctr : indices( inter_planar_angles.size() ) ) {
		if ( inter_planar_angles.size() - planar_ctr != inter_planar_angles[planar_ctr].size() ) {
			BOOST_THROW_EXCEPTION(invalid_argument_exception("The number of entries in one of the inter_planar_angles does not match expected"));
		}
	}
}

/// \brief The number of sec_file_records
sec_file::size_type sec_file::size() const {
	return records.size();
}

/// \brief TODOCUMENT
const sec_struc_planar_angles & sec_file::get_planar_angles_of_indices(const size_t &arg_index_first, ///< TODOCUMENT
                                                                       const size_t &arg_index_second ///< TODOCUMENT
                                                                       ) const {
	// Sanity check the inputs
	if (arg_index_first >= arg_index_second) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("The first index must be strictly less than the second"));
	}
	if (arg_index_second >= records.size()) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("The indices must both be less than the number of sec records"));
	}

	// Return the desired inter_planar_angles
	return inter_planar_angles[arg_index_first][arg_index_second - arg_index_first - 1];

}

/// \brief Standard begin() operator to allow iteration over the sec_file_records
sec_file::const_iterator sec_file::begin() const {
	return common::cbegin( records );
}

/// \brief Standard end() operator to allow iteration over the sec_file_records
sec_file::const_iterator sec_file::end() const {
	return common::cend( records );
}

/// \brief Non-member, non-friend function for converting a sec_file to a sec_struc_vec
///
/// \relates sec_file
sec_struc_vec cath::file::make_sec_struc_list(const sec_file &arg_sec_file ///< The sec_file to convert
                                              ) {
	// Build a vector of sec_strucs from the vector of sec_file_records
	sec_struc_vec new_sec_strucs;
	new_sec_strucs.reserve(arg_sec_file.size());
	for (const sec_file_record &the_sec_file_record : arg_sec_file) {
		new_sec_strucs.push_back(make_sec_struc(the_sec_file_record));
	}

	// Populate the sec_struc_planar_angles of the sec_strucs
	//
	// This code duplicates the behaviour of the original SSAP code when populating the sec_struc_planar_angles for
	// each sec_struc:
	//  - as expected, for the sec_struc_planar_angles to later sec_strucs, use the entries read from the sec file
	//  - for the sec_struc_planar_angles to the sec_struc in question, use a dummy sec_struc_planar_angles
	//    (sec_struc_planar_angles::NULL_SEC_STRUC_PLANAR_ANGLES)
	//  - for the sec_struc_planar_angles to preceding sec_strucs, just copy (without alteration) the sec_struc_planar_angles from
	//    that sec_struc to this
	for (const size_t &sec_struc_ctr_a : indices( new_sec_strucs.size() ) ) {
		sec_struc_planar_angles_vec new_planar_angles;
		new_planar_angles.reserve(new_sec_strucs.size());
		for (const size_t &sec_struc_ctr_b : indices( new_sec_strucs.size() ) ) {
			if (sec_struc_ctr_a == sec_struc_ctr_b) {
				new_planar_angles.push_back(sec_struc_planar_angles::NULL_SEC_STRUC_PLANAR_ANGLES);
			}
			else {
				new_planar_angles.push_back(arg_sec_file.get_planar_angles_of_indices(
					min(sec_struc_ctr_a, sec_struc_ctr_b),
					max(sec_struc_ctr_a, sec_struc_ctr_b)
				));
			}
		}
		new_sec_strucs[sec_struc_ctr_a].set_planar_angles(new_planar_angles);
	}

	return new_sec_strucs;
}


/// \brief Calculate the planar angles corresponding to the specified vector of sec_file_records
///
/// The resulting data structure is as in sec_file, ie for each sec_file_record, there is a vector
/// of the angles between that sec_file_record and each of the following sec_file_records.
///
/// \relates sec_file
sec_struc_planar_angles_vec_vec cath::file::calc_planar_angles(const sec_file_record_vec &arg_sec_file_records ///< The sec_file_records for which to calculate the planar_angles
                                                               ) {
	if ( arg_sec_file_records.empty() ) {
		return {};
	}
	return transform_build<sec_struc_planar_angles_vec_vec>(
		indices( arg_sec_file_records.size() - 1 ),
		[&] (const size_t &x) {
			return transform_build<sec_struc_planar_angles_vec>(
				irange( x + 1_z, arg_sec_file_records.size() ),
				[&] (const size_t &y) {
					return make_planar_angles(
						arg_sec_file_records[ x ],
						arg_sec_file_records[ y ]
					);
				}
			);
		}
	);
}

/// \brief Make a sec_file from the specified vector of sec_file_records by calculating the planar_angles
///
/// \relates sec_file
sec_file cath::file::make_sec_file_with_calced_planar_angles(const sec_file_record_vec &arg_sec_file_records ///< The sec_file_records from which to build the sec_file
                                                             ) {
	return {
		arg_sec_file_records,
		calc_planar_angles( arg_sec_file_records )
	};
}
