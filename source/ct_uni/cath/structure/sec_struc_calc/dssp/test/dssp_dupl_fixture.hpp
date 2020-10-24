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

#ifndef _CATH_TOOLS_SOURCE_CT_UNI_CATH_STRUCTURE_SEC_STRUC_CALC_DSSP_TEST_DSSP_DUPL_FIXTURE_HPP
#define _CATH_TOOLS_SOURCE_CT_UNI_CATH_STRUCTURE_SEC_STRUC_CALC_DSSP_TEST_DSSP_DUPL_FIXTURE_HPP

#include <boost/filesystem/path.hpp>
#include <boost/optional.hpp>

#include "cath/biocore/residue_name.hpp"
#include "cath/common/type_aliases.hpp"

#include <utility>
#include <vector>

namespace cath { namespace sec { class bifur_hbond; } }
namespace cath { namespace sec { class bifur_hbond_list; } }
namespace cath { namespace sec { struct hbond_half; } }
namespace cath { namespace sec { using hbond_half_opt = boost::optional<hbond_half>; } }

namespace cath {
	namespace sec {

		/// \brief Represent an h-bond as parsed from a DSSP file
		using dsspfile_hbond     = std::pair<int, double>;

		/// \brief Type alias for an optional dsspfile_hbond
		using dsspfile_hbond_opt = boost::optional<dsspfile_hbond>;

		/// \brief Represent the hbond information stored on one DSSP line for one residue
		struct dssp_dupl_res final {

			/// \brief The first and second best hbonds from the NH of this residue
			std::pair<dsspfile_hbond_opt, dsspfile_hbond_opt> hbonds_this_nh_1st_2nd;

			/// \brief The first and second best hbonds from the CO of this residue
			std::pair<dsspfile_hbond_opt, dsspfile_hbond_opt> hbonds_this_co_1st_2nd;

			/// \brief The index of this residue in the DSSP file
			size_t residue_index;

			/// \brief The residue_name of this residue
			residue_name pdb_residue_name;
		};

		/// \brief Make a null dssp_dupl_res
		inline dssp_dupl_res make_null_dssp_dupl_res() {
			return {
				std::make_pair( boost::none, boost::none ),
				std::make_pair( boost::none, boost::none ),
				0,
				residue_name{},
			};
		}

		/// \brief Type alias for a vector of dssp_dupl_res
		using dssp_dupl_res_vec = std::vector<dssp_dupl_res>;

		hbond_half_opt mapped_dsspfile_hbond(const dsspfile_hbond_opt &,
		                                     const size_t &,
		                                     const size_vec &);

		str_opt difference_string(const std::string &,
		                          const hbond_half_opt &,
		                          const hbond_half_opt &);
		str_opt difference_string(const dssp_dupl_res &,
		                          const bifur_hbond &,
		                          const size_vec &);
		str_opt difference_string(const dssp_dupl_res_vec &,
		                          const bifur_hbond_list &);

		/// \brief Fixture to assist in testing duplication of DSSP functionality
		class dssp_dupl_fixture {
		protected:
			~dssp_dupl_fixture() noexcept = default;

			dssp_dupl_res_vec parse_dssp_for_calc_testing(const boost::filesystem::path &);
			dssp_dupl_res_vec parse_dssp_for_calc_testing(std::istream &);
			std::pair<size_t, dssp_dupl_res> parse_dssp_residue_line(const std::string &);

			dsspfile_hbond_opt parse_dsspfile_bond(const std::string &);

			static const boost::filesystem::path & DSSP_ROOT_TEST_DATA_DIR();
			static const boost::filesystem::path & DSSP_HBOND_TEST_DATA_DIR();
			static const boost::filesystem::path & DSSP_SS_TEST_DATA_DIR();
		};

	} // namespace sec
} // namespace cath
#endif // _CATH_TOOLS_SOURCE_CT_UNI_CATH_STRUCTURE_SEC_STRUC_CALC_DSSP_TEST_DSSP_DUPL_FIXTURE_HPP
