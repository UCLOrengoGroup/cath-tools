/// \file
/// \brief The dssp_hbond_calc class definitions

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

#include "dssp_hbond_calc.h"

#include <boost/range/irange.hpp>

#include "common/debug_numeric_cast.h"
#include "common/size_t_literal.h"
#include "file/pdb/pdb_atom.h"
#include "structure/sec_struc_calc/bifur_hbond_list.h"

using namespace cath::common;
using namespace cath::sec;
using namespace cath::file;

using boost::irange;
using std::make_pair;

/// \brief Calculate the bifur_hbond_list list of (possibly bifurcating) hbonds between
///        the residues in the specified PDB
///
/// This calls backbone_complete_subset_of_pdb() on the PDB. If the results are being generated
/// in a context where that is required for other tasks, it's better to call backbone_complete_subset_of_pdb()
/// outside and then use calc_bifur_hbonds_of_backbone_complete_pdb() instead.
bifur_hbond_list dssp_hbond_calc::calc_bifur_hbonds_of_pdb__recalc_backbone_residues(const pdb             &arg_pdb,            ///< The PDB to query
                                                                                     const ostream_ref_opt &arg_ostream_ref_opt ///< An optional reference to an ostream to which any logging should be sent
                                                                                     ) {
	// constexpr bool SKIP_RESIDUES_DSSP_WOULD_SKIP = false;
	constexpr bool SKIP_RESIDUES_DSSP_WOULD_SKIP = true;
	return calc_bifur_hbonds_of_backbone_complete_pdb(
		backbone_complete_subset_of_pdb(
			arg_pdb,
			arg_ostream_ref_opt,
			SKIP_RESIDUES_DSSP_WOULD_SKIP
		)
	);
}

/// \brief Calculate the bifur_hbond_list list of (possibly bifurcating) hbonds between
///        the residues in the specified PDB
///
/// \pre The specified PDB must already contain only backbone-complete residues
///
/// For simplicity, calc_bifur_hbonds_of_pdb__recalc_backbone_residues() can be used
/// with a non backbone-complete PDB
bifur_hbond_list dssp_hbond_calc::calc_bifur_hbonds_of_backbone_complete_pdb(const pdb &arg_pdb ///< The PDB to query
                                                                             ) {
	const size_t num_pdb_residues = arg_pdb.get_num_residues();

	bifur_hbond_list results{ num_pdb_residues };

	for (const size_t &i : irange( 0_z, num_pdb_residues ) ) {

		for (const size_t &j : irange( i + 1, num_pdb_residues ) ) {

			for (const auto &indices : { make_pair( i, j ), make_pair( j, i ) } ) {
				const auto &index_1 = indices.first;
				const auto &index_2 = indices.second;

				if ( dssp_hbond_calc::has_hbond_energy_asymm( arg_pdb, index_1, index_2 ) ) {
					const auto energy = dssp_hbond_calc::get_hbond_energy_asymm(
						arg_pdb,
						index_1,
						index_2
					);
					if ( energy < 0.0 ) {
						results.update_with_nh_idx_co_idx_energy(
							debug_numeric_cast<hbond_partner_t>( index_1  ),
							debug_numeric_cast<hbond_partner_t>( index_2 ),
							energy
						);
					}
				}
			}
		}
	}

	return results;
}
