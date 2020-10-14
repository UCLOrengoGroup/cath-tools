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

#include "protein_list.hpp"

#include <boost/algorithm/minmax_element.hpp>

#include "common/algorithm/transform_build.hpp"
#include "common/boost_addenda/minmax_element.hpp"
#include "common/cpp14/cbegin_cend.hpp"
//#include "file/pdb/pdb.hpp"
//#include "file/pdb/pdb_atom.hpp"
//#include "file/pdb/pdb_residue.hpp"
#include "structure/protein/protein.hpp"
#include "structure/protein/residue.hpp"
#include "structure/protein/sec_struc.hpp"
#include "structure/protein/sec_struc_planar_angles.hpp"

using namespace cath;
using namespace cath::common;
using namespace std;

/// \brief TODOCUMENT
void protein_list::push_back(const protein &prm_pdb ///< TODOCUMENT
                             ) {
	proteins.push_back(prm_pdb);
}

/// \brief TODOCUMENT
void protein_list::reserve(const size_t &prm_size ///< TODOCUMENT
                           ) {
	proteins.reserve(prm_size);
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
protein & protein_list::operator[](const size_t &prm_index ///< TODOCUMENT
                                   ) {
	return proteins[prm_index];
}

/// \brief TODOCUMENT
const protein & protein_list::operator[](const size_t &prm_index ///< TODOCUMENT
                                         ) const {
	return proteins[prm_index];
}

/// \brief Standard non-const begin() operator to provide range access
auto protein_list::begin() -> iterator{
	return std::begin( proteins );
}

/// \brief Standard non-const end() operator to provide range access
auto protein_list::end() -> iterator{
	return std::end( proteins );
}

/// \brief Standard const begin() operator to provide range access
auto protein_list::begin() const -> const_iterator {
	return common::cbegin( proteins );
}

/// \brief Standard const end() operator to provide range access
auto protein_list::end() const -> const_iterator {
	return common::cend( proteins );
}

/// \brief TODOCUMENT
protein_list cath::make_protein_list(const protein_vec &prm_proteins ///< TODOCUMENT
                                     ) {
	protein_list new_protein_list;
	new_protein_list.reserve( prm_proteins.size() );
	for ( const protein &the_protein : prm_proteins ) {
		new_protein_list.push_back( the_protein );
	}
	return new_protein_list;
}

/// \brief TODOCUMENT
///
/// \todo Consider retiring this and changing Boost Range adaptor style approach.
///       (At present, the only calls to this are commented out.)
protein_list cath::make_subset_protein_list(const protein_list &prm_proteins, ///< TODOCUMENT
                                            const size_vec     &prm_indices   ///< TODOCUMENT
                                            ) {
	const size_t num_proteins = prm_proteins.size();

	protein_vec new_proteins;
	new_proteins.reserve( prm_indices.size() );

	for (const size_t &index : prm_indices) {
		if ( index >= num_proteins ) {
			BOOST_THROW_EXCEPTION(invalid_argument_exception(""));
		}
		new_proteins.push_back( prm_proteins[ index ] );
	}
	return make_protein_list( new_proteins );
}

/// \brief TODOCUMENT
///
/// \relates protein_list
amino_acid_vec_vec cath::get_amino_acid_lists(const protein_list &prm_proteins ///< TODOCUMENT
                                              ) {
	return transform_build<amino_acid_vec_vec>(
		prm_proteins,
		[] (const protein &x) { return get_amino_acid_list( x ); }
	);
}

/// \brief TODOCUMENT
///
/// \relates protein_list
size_vec cath::get_protein_lengths(const protein_list &prm_proteins ///< TODOCUMENT
                                   ) {
	return transform_build<size_vec>(
		prm_proteins,
		[] (const protein &x) { return x.get_length(); }
	);
}

/// \brief TODOCUMENT
///
/// \relates protein_list
size_size_pair cath::min_max_protein_length(const protein_list &prm_protein_list ///< TODOCUMENT
                                            ) {
	if ( prm_protein_list.empty() ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot calculate min_max_protein_length() for empty protein_list"));
	}
	const auto minmax_iters = common::minmax_element(
		prm_protein_list,
		[] (const protein &x, const protein &y) { return x.get_length() < y.get_length(); }
	);
	return make_pair(
		minmax_iters.first->get_length(),
		minmax_iters.second->get_length()
	);
}
