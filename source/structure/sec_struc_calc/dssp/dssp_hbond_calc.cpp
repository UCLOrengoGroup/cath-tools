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

#include "dssp_hbond_calc.hpp"

#include <boost/range/irange.hpp>

#include "common/debug_numeric_cast.hpp"
#include "common/size_t_literal.hpp"
#include "file/pdb/pdb_atom.hpp"
#include "scan/spatial_index/spatial_index.hpp"
#include "structure/sec_struc_calc/dssp/bifur_hbond_list.hpp"

using namespace cath::common;
using namespace cath::file;
using namespace cath::scan;
using namespace cath::sec;

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
	return calc_bifur_hbonds_of_backbone_complete_pdb(
		backbone_complete_subset_of_pdb(
			arg_pdb,
			arg_ostream_ref_opt,
			dssp_skip_res_skipping::SKIP
		).first
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
	// Note that this is set a little bit higher because sometimes float rounding
	// errors take the answer over the cutoff. This can be fixed using doubles
	// in the distance calculation (along with cast_pdb_coord_float_to_double() )
	// but it's simpler to just increase this number slightly, which wont break
	// results because they'll still all be checked against the correct value of
	// MIN_NO_HBOND_CA_DIST in has_hbond_energy()
	constexpr float MAX_DIST  =  9.03125; // 9 + 1/32
	constexpr float CELL_SIZE = 18.0;
	const size_t num_pdb_residues = arg_pdb.get_num_residues();

	bifur_hbond_list results{ num_pdb_residues };

	const auto lattice = make_sparse_lattice( arg_pdb, CELL_SIZE );

	scan_sparse_lattice(
		lattice,
		arg_pdb,
		CELL_SIZE,
		MAX_DIST,
		[&] (const simple_locn_index &x, const simple_locn_index &y) {
			if ( x.index != y.index ) {
				if ( dssp_hbond_calc::has_hbond_energy_asymm( arg_pdb, x.index, y.index ) ) {
					const auto energy = dssp_hbond_calc::get_hbond_energy_asymm( arg_pdb, x.index, y.index );
					if ( energy < 0.0 ) {
						results.update_with_nh_idx_co_idx_energy(
							x.index,
							y.index,
							energy
						);
					}
				}
			}
		}
	);

	return results;
}
