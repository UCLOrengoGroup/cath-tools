/// \file
/// \brief The protein_list class definitions

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

#include "protein_list.h"

#include <boost/algorithm/minmax_element.hpp>

#include "common/algorithm/transform_build.h"
#include "common/boost_addenda/minmax_element.h"
#include "common/c++14/cbegin_cend.h"
//#include "file/pdb/pdb.h"
//#include "file/pdb/pdb_atom.h"
//#include "file/pdb/pdb_residue.h"
#include "structure/protein/protein.h"
#include "structure/protein/residue.h"
#include "structure/protein/sec_struc.h"
#include "structure/protein/sec_struc_planar_angles.h"

using namespace cath;
using namespace cath::common;
using namespace std;

/// \brief TODOCUMENT
void protein_list::push_back(const protein &arg_pdb ///< TODOCUMENT
                             ) {
	proteins.push_back(arg_pdb);
}

/// \brief TODOCUMENT
void protein_list::reserve(const size_t &arg_size ///< TODOCUMENT
                           ) {
	proteins.reserve(arg_size);
}

/// \brief TODOCUMENT
size_t protein_list::size() const noexcept {
	return proteins.size();
}

/// \brief TODOCUMENT
size_t protein_list::max_size() const noexcept {
	return proteins.max_size();
}

/// \brief TODOCUMENT
bool protein_list::empty() const {
	return proteins.empty();
}

/// \brief TODOCUMENT
protein & protein_list::operator[](const size_t &arg_index ///< TODOCUMENT
                                   ) {
	return proteins[arg_index];
}

/// \brief TODOCUMENT
const protein & protein_list::operator[](const size_t &arg_index ///< TODOCUMENT
                                         ) const {
	return proteins[arg_index];
}

/// \brief TODOCUMENT
protein_list::const_iterator protein_list::begin() const {
	return common::cbegin( proteins );
}
/// \brief TODOCUMENT
protein_list::const_iterator protein_list::end() const {
	return common::cend( proteins );
}

/// \brief TODOCUMENT
protein_list cath::make_protein_list(const protein_vec &arg_proteins ///< TODOCUMENT
                                     ) {
	protein_list new_protein_list;
	new_protein_list.reserve( arg_proteins.size() );
	for ( const protein &the_protein : arg_proteins ) {
		new_protein_list.push_back( the_protein );
	}
	return new_protein_list;
}

/// \brief TODOCUMENT
///
/// \todo Consider retiring this and changing Boost Range adaptor style approach.
///       (At present, the only calls to this are commented out.)
protein_list cath::make_subset_protein_list(const protein_list &arg_proteins, ///< TODOCUMENT
                                            const size_vec     &arg_indices   ///< TODOCUMENT
                                            ) {
	const size_t num_proteins = arg_proteins.size();

	protein_vec new_proteins;
	new_proteins.reserve( arg_indices.size() );

	for (const size_t &index : arg_indices) {
		if ( index >= num_proteins ) {
			BOOST_THROW_EXCEPTION(invalid_argument_exception(""));
		}
		new_proteins.push_back( arg_proteins[ index ] );
	}
	return make_protein_list( new_proteins );
}

/// \brief TODOCUMENT
///
/// \relates protein_list
amino_acid_vec_vec cath::get_amino_acid_lists(const protein_list &arg_proteins ///< TODOCUMENT
                                              ) {
	return transform_build<amino_acid_vec_vec>(
		arg_proteins,
		[] (const protein &x) { return get_amino_acid_list( x ); }
	);
}

/// \brief TODOCUMENT
///
/// \relates protein_list
size_vec cath::get_protein_lengths(const protein_list &arg_proteins ///< TODOCUMENT
                                   ) {
	return transform_build<size_vec>(
		arg_proteins,
		[] (const protein &x) { return x.get_length(); }
	);
}

/// \brief TODOCUMENT
///
/// \relates protein_list
size_size_pair cath::min_max_protein_length(const protein_list &arg_protein_list ///< TODOCUMENT
                                            ) {
	if ( arg_protein_list.empty() ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot calculate min_max_protein_length() for empty protein_list"));
	}
	const auto minmax_iters = common::minmax_element(
		arg_protein_list,
		[] (const protein &x, const protein &y) { return x.get_length() < y.get_length(); }
	);
	return make_pair(
		minmax_iters.first->get_length(),
		minmax_iters.second->get_length()
	);
}
