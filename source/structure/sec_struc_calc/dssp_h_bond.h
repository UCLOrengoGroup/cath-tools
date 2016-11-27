/// \file
/// \brief The dssp_h_bond class header

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

#ifndef _CATH_TOOLS_SOURCE_STRUCTURE_SEC_STRUC_CALC_DSSP_H_BOND_H
#define _CATH_TOOLS_SOURCE_STRUCTURE_SEC_STRUC_CALC_DSSP_H_BOND_H

#include "file/pdb/pdb.h"
#include "file/pdb/pdb_residue.h"
#include "structure/geometry/coord.h"
#include "structure/protein/amino_acid.h"

namespace cath {
	namespace sec {

		/// \brief TODOCUMENT
		//
		//  * http://swift.cmbi.ru.nl/gv/dssp/
		//  * https://en.wikipedia.org/wiki/DSSP_(hydrogen_bond_estimation_algorithm)
		//  * https://pdfs.semanticscholar.org/2fcd/a810f3f43cd1aaeea04782ae4234cb61a3f8.pdf
		//
		// "Continuum secondary structure captures protein flexibility" :
		//
		//  * http://www.cell.com/structure/fulltext/S0969-2126(02)00700-1
		//
		//
		//
		// Place a bond between C==O of i and N-H of j if
		// H bond if E < -0.5kcal/mol:
		// Where E = 0.0084 {1/r_on + 1/r_ch - 1/r_oh - 1/r_cn} . 332 kcal/mol
		//
		// 27888 / 10000 {1/r_on + 1/r_ch - 1/r_oh - 1/r_cn} < -0.5
		//  1743 /   625 {1/r_on + 1/r_ch - 1/r_oh - 1/r_cn} < -   1 /    2
		//               {1/r_on + 1/r_ch - 1/r_oh - 1/r_cn} < - 625 / 3486
		//    15 1013 A G  S  > S-     0   0   12    105,-0.5     4,-2.7    24,-0.0     5,-0.3  -0.182  92.8 -74.3 -78.3-173.8   27.0  135.0   33.1
		//    16 1014 A V  H  > S+     0   0    4      1,-0.2     4,-2.1     2,-0.2     5,-0.2   0.849 134.0  50.4 -58.2 -33.5   24.8  132.6   31.2
		//    17 1015 A I  H  > S+     0   0   17      2,-0.2     4,-2.6     1,-0.2     5,-0.2   0.955 112.8  43.1 -71.6 -46.9   24.4  130.4   34.2
		//    18 1016 A G  H  > S+     0   0    0    102,-0.3     4,-2.7     1,-0.2    -2,-0.2   0.919 116.9  48.7 -62.2 -43.1   23.4  133.1   36.7
		//    19 1017 A L  H  X S+     0   0    0     -4,-2.7     4,-2.3     2,-0.2    -1,-0.2   0.929 113.3  45.6 -62.6 -47.2   21.0  134.6   34.1
		//    20 1018 A S  H  X S+     0   0    0     -4,-2.1     4,-2.2    -5,-0.3     5,-0.2   0.929 114.7  48.1 -64.7 -42.1   19.4  131.3   33.1
		//    21 1019 A S  H  X S+     0   0    0     -4,-2.6     4,-2.2     1,-0.2    -2,-0.2   0.945 112.2  49.5 -62.0 -47.2   19.0  130.3   36.8
		//    22 1020 A A  H  X S+     0   0    0     -4,-2.7     4,-2.8    -5,-0.2    -1,-0.2   0.888 109.3  53.1 -56.7 -42.5   17.5  133.7   37.6
		//    23 1021 A L  H  X S+     0   0    3     -4,-2.3     4,-2.4     2,-0.2    -1,-0.2   0.930 108.9  45.8 -67.3 -45.9   15.1  133.5   34.7
		//    24 1022 A I  H  X S+     0   0   16     -4,-2.2     4,-1.6     1,-0.2    -1,-0.2   0.936 114.8  50.0 -62.5 -41.6   13.6  130.1   35.7
		//    25 1023 A L  H  <>S+     0   0    0     -4,-2.2     5,-2.4    -5,-0.2     4,-0.3   0.930 110.0  49.7 -60.8 -44.2   13.3  131.3   39.3
		//    26 1024 A A  H ><5S+     0   0   11     -4,-2.8     3,-1.5     1,-0.2    -1,-0.2   0.917 109.1  52.2 -61.8 -40.4   11.6  134.5   38.2
		//    27 1025 A R  H 3<5S+     0   0  149     -4,-2.4    -1,-0.2     1,-0.3    -2,-0.2   0.833 105.3  56.3 -62.7 -32.0    9.2  132.5   36.1
		//    28 1026 A K  T 3<5S-     0   0   85     -4,-1.6    -1,-0.3    -5,-0.2    -2,-0.2   0.434 127.3-100.6 -80.4  -0.9    8.4  130.4   39.2
		//
		//
		//    18 1016 A G  H  > S+     0   0    0    102,-0.3     4,-2.7     1,-0.2    -2,-0.2   0.919 116.9  48.7 -62.2 -43.1   23.4  133.1   36.7
		//    22 1020 A A  H  X S+     0   0    0     -4,-2.7     4,-2.8    -5,-0.2    -1,-0.2   0.888 109.3  53.1 -56.7 -42.5   17.5  133.7   37.6
		//
		// {1/r_on + 1/r_ch - 1/r_oh - 1/r_cn} = - 1125 / 1162 = -0.968158347676
		//
		// ATOM    131  N   GLY A1016      24.098 132.249  35.739  1.00  9.10
		// ATOM    132  CA  GLY A1016      23.354 133.075  36.666  1.00  8.87
		// ATOM    133  C   GLY A1016      22.093 133.611  36.006  1.00  8.83
		// ATOM    134  O   GLY A1016      21.020 133.627  36.607  1.00  9.33

		// ATOM    155  N   ALA A1020      18.386 132.630  37.199  1.00  9.12
		// ATOM    156  CA  ALA A1020      17.516 133.744  37.580  1.00  9.37
		// ATOM    157  C   ALA A1020      16.225 133.709  36.827  1.00  9.61
		// ATOM    158  O   ALA A1020      15.148 133.956  37.406  1.00 10.41
		// ATOM    159  CB  ALA A1020      18.226 135.061  37.384  1.00  8.97

		// // assign the Hydrogen
		// mH = GetN();
		// 
		// if (mType != kProline and mPrev != nullptr)
		// {
		// 	const MAtom& pc = mPrev->GetC();
		// 	const MAtom& po = mPrev->GetO();
		// 	
		// 	double CODistance = Distance(pc, po);
		// 	
		// 	mH.mLoc.mX += (pc.mLoc.mX - po.mLoc.mX) / CODistance; 
		// 	mH.mLoc.mY += (pc.mLoc.mY - po.mLoc.mY) / CODistance; 
		// 	mH.mLoc.mZ += (pc.mLoc.mZ - po.mLoc.mZ) / CODistance; 
		// }


		// kCouplingConstant  = -27.888, // = -332 * 0.42 * 0.2
		// kMinHBondEnergy    =  -9.9
		// kMinimalCADistance =   9.0
		// kMinimalDistance   =   0.5


		// // Calculate the HBond energies
		// for (uint32 i = 0; i + 1 < inResidues.size(); ++i)
		// {
		// 	MResidue* ri = inResidues[i];
		// 	
		// 	for (uint32 j = i + 1; j < inResidues.size(); ++j)
		// 	{
		// 		MResidue* rj = inResidues[j];
		// 		
		// 		if (Distance(ri->GetCAlpha(), rj->GetCAlpha()) < kMinimalCADistance)
		// 		{
		// 			MResidue::CalculateHBondEnergy(*ri, *rj);
		// 			if (j != i + 1)
		// 				MResidue::CalculateHBondEnergy(*rj, *ri);
		// 		}
		// 	}
		// }
		// first < second || first > second + 1
		// second is bigger than first or first


		//
		// TODO: use the angle to improve bond energy calculation.
		// double MResidue::CalculateHBondEnergy(MResidue& inDonor, MResidue& inAcceptor)
		// {
		// 	double result = 0;
		// 	
		// 	if (inDonor.mType != kProline)
		// 	{
		// 		double distanceHO = Distance(inDonor.GetH(), inAcceptor.GetO());
		// 		double distanceHC = Distance(inDonor.GetH(), inAcceptor.GetC());
		// 		double distanceNC = Distance(inDonor.GetN(), inAcceptor.GetC());
		// 		double distanceNO = Distance(inDonor.GetN(), inAcceptor.GetO());
		// 		
		// 		if (distanceHO < kMinimalDistance or distanceHC < kMinimalDistance or distanceNC < kMinimalDistance or distanceNO < kMinimalDistance)
		// 			result = kMinHBondEnergy;
		// 		else
		// 			result = kCouplingConstant / distanceHO - kCouplingConstant / distanceHC + kCouplingConstant / distanceNC - kCouplingConstant / distanceNO;
		// 
		// 		// DSSP compatibility mode:
		// 		result = bm::round(result * 1000) / 1000;
		// 
		// 		if (result < kMinHBondEnergy)
		// 			result = kMinHBondEnergy;
		// 	}
		// 
		// 	// update donor
		// 	if (result < inDonor.mHBondAcceptor[0].energy)
		// 	{
		// 		inDonor.mHBondAcceptor[1] = inDonor.mHBondAcceptor[0];
		// 		inDonor.mHBondAcceptor[0].residue = &inAcceptor;
		// 		inDonor.mHBondAcceptor[0].energy = result;
		// 	}
		// 	else if (result < inDonor.mHBondAcceptor[1].energy)
		// 	{
		// 		inDonor.mHBondAcceptor[1].residue = &inAcceptor;
		// 		inDonor.mHBondAcceptor[1].energy = result;
		// 	}		
		// 
		// 	// and acceptor
		// 	if (result < inAcceptor.mHBondDonor[0].energy)
		// 	{
		// 		inAcceptor.mHBondDonor[1] = inAcceptor.mHBondDonor[0];
		// 		inAcceptor.mHBondDonor[0].residue = &inDonor;
		// 		inAcceptor.mHBondDonor[0].energy = result;
		// 	}
		// 	else if (result < inAcceptor.mHBondDonor[1].energy)
		// 	{
		// 		inAcceptor.mHBondDonor[1].residue = &inDonor;
		// 		inAcceptor.mHBondDonor[1].energy = result;
		// 	}		
		// 	
		// 	return result;
		// }
		class dssp_h_bond final {
		public:
			dssp_h_bond() = delete;
			~dssp_h_bond() = delete;

			static double get_h_bond_energy(const geom::coord &,
			                                const geom::coord &,
			                                const geom::coord &,
			                                const geom::coord &);

			static double get_h_bond_energy(const file::pdb_residue &,
			                                const file::pdb_residue &,
			                                const file::pdb_residue &);

			static double get_h_bond_energy_asymm(const file::pdb &,
			                                      const size_t &,
			                                      const size_t &);

			static bool has_h_bond_energy(const file::pdb_residue &,
			                              const file::pdb_residue &,
			                              const file::pdb_residue &);

			static bool has_h_bond_energy_asymm(const file::pdb &,
			                                    const size_t &,
			                                    const size_t &);
		};

		/// \brief TODOCUMENT
		inline double dssp_h_bond::get_h_bond_energy(const geom::coord &arg_c, ///< TODOCUMENT
		                                             const geom::coord &arg_o, ///< TODOCUMENT
		                                             const geom::coord &arg_n, ///< TODOCUMENT
		                                             const geom::coord &arg_h  ///< TODOCUMENT
		                                             ) {
			constexpr double ENERGY_MULTIPLIER = 0.42 * 0.2 * 332;
			constexpr double MIN_ENERGY        = -9.9;
			const double energy = ENERGY_MULTIPLIER * (
				  ( 1.0 / distance_between_points( arg_o, arg_n ) ) 
				+ ( 1.0 / distance_between_points( arg_c, arg_h ) )
				- ( 1.0 / distance_between_points( arg_o, arg_h ) )
				- ( 1.0 / distance_between_points( arg_c, arg_n ) )
			);
			return ( energy > MIN_ENERGY ) ? std::round( energy * 1000.0 ) / 1000
			                               : MIN_ENERGY;
		}

		/// \brief TODOCUMENT
		inline double dssp_h_bond::get_h_bond_energy(const file::pdb_residue &arg_residue_i,      ///< TODOCUMENT
		                                             const file::pdb_residue &arg_residue_j_prev, ///< TODOCUMENT
		                                             const file::pdb_residue &arg_residue_j       ///< TODOCUMENT
		                                             ) {
			const geom::coord prev_c_to_o = get_oxygen_coord( arg_residue_j_prev )
			                                -
			                                get_carbon_coord( arg_residue_j_prev );
			return get_h_bond_energy(
				get_carbon_coord  ( arg_residue_i ),
				get_oxygen_coord  ( arg_residue_i ),
				get_nitrogen_coord( arg_residue_j ),
				get_nitrogen_coord( arg_residue_j ) - ( prev_c_to_o / length( prev_c_to_o ) )
			);
		}

		/// \brief TODOCUMENT
		inline double dssp_h_bond::get_h_bond_energy_asymm(const file::pdb &arg_pdb, ///< TODOCUMENT
		                                                   const size_t    &arg_i,   ///< TODOCUMENT
		                                                   const size_t    &arg_j    ///< TODOCUMENT
		                                                   ) {
			return get_h_bond_energy(
				arg_pdb.get_residue_cref_of_index__backbone_unchecked( arg_i     ),
				arg_pdb.get_residue_cref_of_index__backbone_unchecked( arg_j - 1 ),
				arg_pdb.get_residue_cref_of_index__backbone_unchecked( arg_j     )
			);
		}

		/// \brief TODOCUMENT
		inline bool dssp_h_bond::has_h_bond_energy(const file::pdb_residue &arg_residue_i,      ///< TODOCUMENT
		                                           const file::pdb_residue &arg_residue_j_prev, ///< TODOCUMENT
		                                           const file::pdb_residue &arg_residue_j       ///< TODOCUMENT
		                                           ) {
			constexpr double MIN_NO_H_BOND_CA_DIST = 9.0;
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
					< MIN_NO_H_BOND_CA_DIST
				)
				&&
				arg_residue_i.has_carbon()
				&&
				arg_residue_i.has_oxygen()
				&&
				arg_residue_j_prev.has_carbon()
				&&
				arg_residue_j_prev.has_oxygen()
				&&
				arg_residue_j.has_nitrogen()
				&&
				get_amino_acid( arg_residue_j ).get_letter() != 'P' // Proline has side-chain on N
			);
		}

		inline bool dssp_h_bond::has_h_bond_energy_asymm(const file::pdb &arg_pdb, ///< TODOCUMENT
		                                                 const size_t    &arg_i,   ///< TODOCUMENT
		                                                 const size_t    &arg_j    ///< TODOCUMENT
		                                                 ) {
			return (
				arg_j >= 1
				&&
				(
					arg_i < arg_j
					||
					arg_i > arg_j + 1
				)
				&&
				has_h_bond_energy(
					arg_pdb.get_residue_cref_of_index__backbone_unchecked( arg_i     ),
					arg_pdb.get_residue_cref_of_index__backbone_unchecked( arg_j - 1 ),
					arg_pdb.get_residue_cref_of_index__backbone_unchecked( arg_j     )
				)
			);
		}


	}
} // namespace cath

#endif
