/// \file
/// \brief The dssp_dupl_fixture header

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

#ifndef _CATH_TOOLS_SOURCE_STRUCTURE_SEC_STRUC_CALC_DSSP_DUPL_FIXTURE_H
#define _CATH_TOOLS_SOURCE_STRUCTURE_SEC_STRUC_CALC_DSSP_DUPL_FIXTURE_H

#include <boost/filesystem/path.hpp>
#include <boost/optional.hpp>

#include "structure/residue_name.h"

#include <utility>
#include <vector>

namespace cath {
	namespace sec {

		using dsspfile_h_bond     = std::pair<int, double>;

		using dsspfile_h_bond_opt = boost::optional<dsspfile_h_bond>;

		struct dssp_dupl_res final {
		// private:
			std::pair<dsspfile_h_bond_opt, dsspfile_h_bond_opt> hbonds_this_nh_1st_2nd;
			std::pair<dsspfile_h_bond_opt, dsspfile_h_bond_opt> hbonds_this_co_1st_2nd;

			size_t residue_index;

			residue_name pdb_residue_name;

		public:
		};

		using dssp_dupl_res_vec = std::vector<dssp_dupl_res>;

		/// \brief TODOCUMENT
		class dssp_dupl_fixture {
		protected:
			dssp_dupl_res_vec parse_dssp_for_calc_testing(const boost::filesystem::path &);
			dssp_dupl_res_vec parse_dssp_for_calc_testing(std::istream &);
			std::pair<size_t, dssp_dupl_res> parse_dssp_residue_line(const std::string &);

			dsspfile_h_bond_opt parse_dsspfile_bond(const std::string &);
		};

	} // namespace sec
} // namespace cath
#endif
