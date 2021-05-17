/// \file
/// \brief The spatial_index class definitions

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

#include "spatial_index.hpp"

#include <boost/range/adaptor/transformed.hpp>

#include "cath/common/boost_addenda/range/indices.hpp"
#include "cath/scan/detail/scan_index_store/scan_index_lattice_store.hpp"
#include "cath/scan/detail/scan_index_store/scan_index_store_helper.hpp"
#include "cath/scan/res_pair_keyer/res_pair_keyer.hpp"
#include "cath/scan/res_pair_keyer/res_pair_keyer_part/detail/axis_keyer_part.hpp"
#include "cath/scan/res_pair_keyer/res_pair_keyer_part/res_pair_view_x_keyer_part.hpp"
#include "cath/scan/res_pair_keyer/res_pair_keyer_part/res_pair_view_y_keyer_part.hpp"
#include "cath/scan/res_pair_keyer/res_pair_keyer_part/res_pair_view_z_keyer_part.hpp"

using namespace ::cath::common;
using namespace ::cath::file;
using namespace ::cath::scan;
using namespace ::cath::scan::detail;

using ::boost::adaptors::transformed;
using ::std::vector;

namespace  {

	/// \brief TODOCUMENT
	template <sod Sod, typename Cell>
	struct simple_spatial_lattice_store_maker final {

		/// \brief TODOCUMENT
		template <typename Rng>
		auto operator()(const Rng   &prm_rng,       ///< TODOCUMENT
		                const float &prm_cell_size, ///< TODOCUMENT
		                const float &prm_max_dist   ///< TODOCUMENT
		                ) {
			const auto keyer = make_res_pair_keyer(
				simple_locn_x_keyer_part{ prm_cell_size },
				simple_locn_y_keyer_part{ prm_cell_size },
				simple_locn_z_keyer_part{ prm_cell_size },
				res_pair_from_to_index_keyer_part{}
			);

			using store_type = scan_index_lattice_store<decltype( keyer )::key_index_tuple_type, Cell>;
			return store_maker<Sod, store_type>{}(
				prm_rng,
				keyer,
				simple_locn_crit{ prm_max_dist * prm_max_dist }
			);
		}
	};

} // namespace


static_assert( std::is_copy_assignable_v<cath::scan::simple_locn_index>, "" );

/// \brief TODOCUMENT
locn_index_store cath::scan::make_sparse_lattice(const protein &prm_protein,   ///< TODOCUMENT
                                                 const float   &prm_cell_size, ///< TODOCUMENT
                                                 const float   &prm_max_dist   ///< TODOCUMENT
                                                 ) {
	return simple_spatial_lattice_store_maker< sod::SPARSE, vector< simple_locn_index > >{}(
		indices( prm_protein.get_length() )
			| transformed( [&] (const size_t &x) {
				return make_simple_locn_index_of_ca( prm_protein.get_residue_ref_of_index( x ), debug_numeric_cast<unsigned int>( x ) );
			} ),
		prm_cell_size,
		prm_max_dist
	);
}

/// \brief TODOCUMENT
locn_index_store cath::scan::make_dense_lattice(const protein &prm_protein,   ///< TODOCUMENT
                                                const float   &prm_cell_size, ///< TODOCUMENT
                                                const float   &prm_max_dist   ///< TODOCUMENT
                                                ) {
	return simple_spatial_lattice_store_maker< sod::DENSE, vector< simple_locn_index > >{}(
		indices( prm_protein.get_length() )
			| transformed( [&] (const size_t &x) {
				return make_simple_locn_index_of_ca( prm_protein.get_residue_ref_of_index( x ), debug_numeric_cast<unsigned int>( x ) );
			} ),
		prm_cell_size,
		prm_max_dist
	);
}

/// \brief TODOCUMENT
locn_index_store cath::scan::make_sparse_lattice(const pdb   &prm_pdb,       ///< TODOCUMENT
                                                 const float &prm_cell_size, ///< TODOCUMENT
                                                 const float &prm_max_dist   ///< TODOCUMENT
                                                 ) {
	return simple_spatial_lattice_store_maker< sod::SPARSE, vector< simple_locn_index > >{}(
		indices( prm_pdb.get_num_residues() )
			| transformed( [&] (const size_t &x) {
				return make_simple_locn_index_of_ca( prm_pdb.get_residue_of_index__backbone_unchecked( x ), debug_numeric_cast<unsigned int>( x ) );
			} ),
		prm_cell_size,
		prm_max_dist
	);
}

/// \brief TODOCUMENT
locn_index_store cath::scan::make_dense_lattice(const pdb   &prm_pdb,       ///< TODOCUMENT
                                                const float &prm_cell_size, ///< TODOCUMENT
                                                const float &prm_max_dist   ///< TODOCUMENT
                                                ) {
	return simple_spatial_lattice_store_maker< sod::DENSE, vector< simple_locn_index > >{}(
		indices( prm_pdb.get_num_residues() )
			| transformed( [&] (const size_t &x) {
				return make_simple_locn_index_of_ca( prm_pdb.get_residue_of_index__backbone_unchecked( x ), debug_numeric_cast<unsigned int>( x ) );
			} ),
		prm_cell_size,
		prm_max_dist
	);
}
