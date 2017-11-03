/// \file
/// \brief The view_cache class definitions

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

#include "view_cache.hpp"


#include "common/boost_addenda/range/indices.hpp"
#include "ssap/context_res.hpp"
#include "ssap/ssap.hpp"
#include "structure/protein/protein.hpp"

using namespace cath::common;
using namespace cath::geom;
using namespace cath::index;


/// \brief Private static method that implements the process of building the views from proteins
coord_vec_vec view_cache::build_views(const protein &arg_protein ///< The protein which the view_cache should be built to represent
                                      ) {
	// Grab the number of residues and prepare the views accordingly
	const size_t num_residues = arg_protein.get_length();
	coord_vec_vec new_views( num_residues );
	for (coord_vec &view_of : new_views) {
		view_of.reserve( num_residues );
	}

	// Loop over the all from-versus-to residue pairs and add the resulting views
	for (const size_t &from_res_ctr : indices( num_residues ) ) {
		for (const size_t &to_res_ctr : indices( num_residues ) ) {
			coord_vec &view_of_from = new_views[ from_res_ctr ];
			view_of_from.push_back( view_vector_of_residue_pair(
				arg_protein.get_residue_ref_of_index( from_res_ctr ),
				arg_protein.get_residue_ref_of_index( to_res_ctr   )
			) );
		}
	}
	return new_views;
}

/// \brief Ctor for view_cache from a protein that the view_cache should represent
view_cache::view_cache(const protein &arg_protein ///< The protein which the view_cache should be built to represent
                       ) : views( build_views( arg_protein ) ) {
}

