/// \file
/// \brief The view_cache_list class definitions

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

#include "view_cache_list.hpp"

#include "ssap/context_res.hpp"
#include "ssap/ssap.hpp"
#include "structure/protein/protein.hpp"
#include "structure/protein/protein_list.hpp"
#include "structure/protein/sec_struc.hpp"
#include "structure/protein/sec_struc_planar_angles.hpp"
#include "structure/view_cache/view_cache.hpp"

using namespace cath::index;

/// \brief Build a vector of view_cache objects from a protein_list object
view_cache_vec view_cache_list::build_caches(const protein_list &prm_protein_list ///< The list of proteins from which the view_cache_vec should be constructed
                                             ) {
	const size_t num_proteins = prm_protein_list.size();

	view_cache_vec new_caches;
	new_caches.reserve( num_proteins );

	for (const protein &the_protein : prm_protein_list) {
		new_caches.push_back( view_cache( the_protein ) );
	}

	return new_caches;
}

/// \brief Ctor for view_cache_list from a protein_list
view_cache_list::view_cache_list(const protein_list &prm_protein_list ///< The list of proteins from which the view_cache_list should be constructed
                                 ) : view_caches( build_caches( prm_protein_list ) ) {
}
