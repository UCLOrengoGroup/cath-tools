/// \file
/// \brief The dssp_hbond_calc class header

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

#ifndef _CATH_TOOLS_SOURCE_UNI_STRUCTURE_SEC_STRUC_CALC_DSSP_DSSP_HBOND_CALC_H
#define _CATH_TOOLS_SOURCE_UNI_STRUCTURE_SEC_STRUC_CALC_DSSP_DSSP_HBOND_CALC_H

#include <boost/algorithm/clamp.hpp>

#include "file/pdb/pdb.hpp"
#include "file/pdb/pdb_residue.hpp"
#include "structure/geometry/coord.hpp"
#include "structure/protein/amino_acid.hpp"

namespace cath { namespace sec { class bifur_hbond_list; } }

// Notes
// -----
//
// Refs:
//
//  * http://swift.cmbi.ru.nl/gv/dssp/
//  * https://en.wikipedia.org/wiki/DSSP_(hydrogen_bond_estimation_algorithm)
//  * https://pdfs.semanticscholar.org/2fcd/a810f3f43cd1aaeea04782ae4234cb61a3f8.pdf
//
// "Continuum secondary structure captures protein flexibility" :
//
//  * http://www.cell.com/structure/fulltext/S0969-2126(02)00700-1
//
// Place a bond between C==O of i and N-H of j if
// H bond if E < -0.5kcal/mol:
// Where E = 0.0084 {1/r_on + 1/r_ch - 1/r_oh - 1/r_cn} . 332 kcal/mol

namespace cath {
	namespace sec {

		/// \brief The type to use for hbond energy calculations
		///
		/// This should be kept as double so that the behaviour directly replicates DSSP's
		using hbond_energy_t = double;

		/// \brief Provide static methods for performing DSSP hbond calculations
		class dssp_hbond_calc final {
		public:
			dssp_hbond_calc() = delete;
			~dssp_hbond_calc() = delete;

			static hbond_energy_t get_hbond_energy(const geom::coord &,
			                                       const geom::coord &,
			                                       const geom::coord &,
			                                       const geom::coord &);

			static hbond_energy_t get_hbond_energy(const boost::optional<file::pdb_residue> &,
			                                       const file::pdb_residue &,
			                                       const file::pdb_residue &);

			static hbond_energy_t get_hbond_energy_asymm(const file::pdb &,
			                                             const size_t &,
			                                             const size_t &);

			static bool has_hbond_energy(const boost::optional<file::pdb_residue> &,
			                             const file::pdb_residue &,
			                             const file::pdb_residue &);

			static bool has_hbond_energy_asymm(const file::pdb &,
			                                   const size_t &,
			                                   const size_t &);

			static bifur_hbond_list calc_bifur_hbonds_of_pdb__recalc_backbone_residues(const file::pdb &,
			                                                                           const ostream_ref_opt & = boost::none);

			static bifur_hbond_list calc_bifur_hbonds_of_backbone_complete_pdb(const file::pdb &);
		};

		/// \brief Calculate the DSSP hbond energy between the specified N & H coords of one
		///        residue and the C & O coords of another
		inline hbond_energy_t dssp_hbond_calc::get_hbond_energy(const geom::coord &arg_n, ///< The N coord of one residue
		                                                        const geom::coord &arg_h, ///< The H coord of one residue
		                                                        const geom::coord &arg_c, ///< The C coord of another residue
		                                                        const geom::coord &arg_o  ///< The O coord of another residue
		                                                        ) {
			constexpr double         MIN_DISTANCE      = 0.5;
			constexpr double         ENERGY_MULTIPLIER = 0.42 * 0.2 * 332;
			constexpr hbond_energy_t MIN_ENERGY        = static_cast<hbond_energy_t>(   -9.9 );
			constexpr hbond_energy_t MAX_ENERGY        = static_cast<hbond_energy_t>(    0.0 );
			constexpr hbond_energy_t ROUNDING_FACTOR   = static_cast<hbond_energy_t>( 1000.0 );

			const double dist_no = distance_between_points( arg_n, arg_o );
			const double dist_hc = distance_between_points( arg_h, arg_c );
			const double dist_ho = distance_between_points( arg_h, arg_o );
			const double dist_nc = distance_between_points( arg_n, arg_c );

			if ( dist_ho < MIN_DISTANCE || dist_hc < MIN_DISTANCE || dist_nc < MIN_DISTANCE || dist_no < MIN_DISTANCE ) {
				return MIN_ENERGY;
			}

			return std::round(
				ROUNDING_FACTOR * boost::algorithm::clamp(
					ENERGY_MULTIPLIER * (
						  ( 1.0 / dist_no )
						+ ( 1.0 / dist_hc )
						- ( 1.0 / dist_ho )
						- ( 1.0 / dist_nc )
					),
					MIN_ENERGY,
					MAX_ENERGY
				)
			) / ROUNDING_FACTOR;
		}

		/// \brief Calculate the DSSP hbond energy between the specified two residues (with supporting info
		///        from the residue that precedes the first one)
		inline hbond_energy_t dssp_hbond_calc::get_hbond_energy(const boost::optional<file::pdb_residue> &arg_residue_i_prev, ///< The residue that precedes the first residue
		                                                        const file::pdb_residue                  &arg_residue_i,      ///< The first residue
		                                                        const file::pdb_residue                  &arg_residue_j       ///< The second residue
		                                                        ) {
			const auto pseudo_h = [&] {
				if ( arg_residue_i_prev ) {
					const auto prev_c_to_o = get_oxygen_coord( *arg_residue_i_prev )
					                         -
					                         get_carbon_coord( *arg_residue_i_prev );
					return get_nitrogen_coord( arg_residue_i ) - ( prev_c_to_o / length( prev_c_to_o ) );
				}
				return get_nitrogen_coord( arg_residue_i );
			} ();
			return get_hbond_energy(
				get_nitrogen_coord( arg_residue_i ),
				pseudo_h,
				get_carbon_coord  ( arg_residue_j ),
				get_oxygen_coord  ( arg_residue_j )
			);
		}

		/// \brief Calculate the DSSP hbond energy between the two residues at the specified indices of the specified PDB
		///
		/// As indicated by the name, this result is asymmetrical; swapping arg_i and arg_j will likely
		/// change the result
		inline hbond_energy_t dssp_hbond_calc::get_hbond_energy_asymm(const file::pdb &arg_pdb, ///< The PDB containing the residues in question
		                                                              const size_t    &arg_i,   ///< The index of the first  residue in the PDB
		                                                              const size_t    &arg_j    ///< The index of the second residue in the PDB
		                                                              ) {
			const auto prev_index = index_of_preceding_residue_in_same_chain( arg_pdb, arg_i );
			return get_hbond_energy(
				prev_index
					? boost::make_optional( arg_pdb.get_residue_of_index__backbone_unchecked( *prev_index ) )
					: boost::none,
				arg_pdb.get_residue_of_index__backbone_unchecked( arg_i     ),
				arg_pdb.get_residue_of_index__backbone_unchecked( arg_j     )
			);
		}

		/// \brief Whether DSSP might assign a valid hbond energy between the specified two residues
		///        (with supporting info from the residue that precedes the first one)
		///
		/// If calling with a PDB and indices, prefer to use the wrapper function below because that can
		/// make additional checks on the indices
		inline bool dssp_hbond_calc::has_hbond_energy(const boost::optional<file::pdb_residue> &arg_residue_i_prev, ///< The residue that precedes the first residue
		                                              const file::pdb_residue                  &arg_residue_i,      ///< The first residue
		                                              const file::pdb_residue                  &arg_residue_j       ///< The second residue
		                                              ) {
			constexpr double MIN_NO_HBOND_CA_DIST = 9.0;

			return (
				arg_residue_i.has_carbon_alpha()
				&&
				arg_residue_j.has_carbon_alpha()
				&&
				(
					distance_between_points(
						get_carbon_alpha_coord( arg_residue_i ),
						get_carbon_alpha_coord( arg_residue_j )
					)
					< MIN_NO_HBOND_CA_DIST
				)
				&&
				(
					! arg_residue_i_prev
					||
					(
						arg_residue_i_prev->has_carbon()
						&&
						arg_residue_i_prev->has_oxygen()
					)
				)
				&&
				arg_residue_i.has_nitrogen()
				&&
				arg_residue_j.has_carbon()
				&&
				arg_residue_j.has_oxygen()
				&&
				(
					! is_proper_amino_acid( arg_residue_i.get_amino_acid() )
					||
					arg_residue_i.get_amino_acid() != amino_acid{ 'P' } // Proline has side-chain on N
				)
			);
		}

		/// \brief Whether DSSP might assign a valid hbond energy between the specified two residues
		///        (with supporting info from the residue that precedes the first one)
		inline bool dssp_hbond_calc::has_hbond_energy_asymm(const file::pdb &arg_pdb, ///< The PDB containing the residues in question
		                                                    const size_t    &arg_i,   ///< The index of the first  residue in the PDB
		                                                    const size_t    &arg_j    ///< The index of the second residue in the PDB
		                                                    ) {
			const auto prev_index = index_of_preceding_residue_in_same_chain( arg_pdb, arg_i );
			return (
				(
					arg_i < arg_j
					||
					arg_i > arg_j + 1
				)
				&&
				has_hbond_energy(
					prev_index
						? boost::make_optional( arg_pdb.get_residue_of_index__backbone_unchecked( *prev_index ) )
						: boost::none,
					arg_pdb.get_residue_of_index__backbone_unchecked( arg_i     ),
					arg_pdb.get_residue_of_index__backbone_unchecked( arg_j     )
				)
			);
		}


	} // namespace sec
} // namespace cath

#endif
