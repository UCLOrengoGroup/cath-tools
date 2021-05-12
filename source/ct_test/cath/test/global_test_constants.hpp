/// \file
/// \brief The global_test_constants class header

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

#ifndef _CATH_TOOLS_SOURCE_CT_TEST_CATH_TEST_GLOBAL_TEST_CONSTANTS_HPP
#define _CATH_TOOLS_SOURCE_CT_TEST_CATH_TEST_GLOBAL_TEST_CONSTANTS_HPP

#include <filesystem>

#include "cath/common/path_type_aliases.hpp"
#include "cath/common/type_aliases.hpp"

namespace cath {

	/// \brief TODOCUMENT
	class global_test_constants {
	private:
		static void check_required_files_exist();
		static void cleanup_temporary_files();

	protected:
		global_test_constants();
		~global_test_constants() noexcept;

		/// \brief Specify that the copy-ctor shouldn't be used
		global_test_constants(const global_test_constants &) = delete;
		/// \brief Specify that the copy-assign shouldn't be used
		global_test_constants & operator=(const global_test_constants &) = delete;

	public:
		static const ::std::filesystem::path              & TEST_BASIC_FILE_TEST_DATA_DIR();
		static const ::std::filesystem::path              & TEST_MULTI_SSAP_SUPERPOSE_DIR();
		static const ::std::filesystem::path              & TEST_RESIDUE_IDS_DATA_DIR();
		static const ::std::filesystem::path              & TEST_SOURCE_DATA_DIR();
		static const ::std::filesystem::path              & TEST_EXAMPLE_PDBS_DATA_DIR();
		static const ::std::filesystem::path              & TEST_SSAP_REGRESSION_DATA_DIR();
		static const ::std::filesystem::path              & TEST_SSAP_ALIGNMENT_GLUING_DATA_DIR();
		static const ::std::filesystem::path              & TEST_SUP_JSON_DIR();
		static const ::std::filesystem::path              & TEST_SVM_DIR();
		static const ::std::filesystem::path              & TEST_OUTPUT_DIRECTORY();
		static const ::std::filesystem::path              & CRH_TEST_DATA_DIR();
		static const ::std::filesystem::path              & CRH_CATH_DC_HANDLING_DATA_DIR();
		static const ::std::filesystem::path              & CRH_HMM_COVERAGE_DATA_DIR();
		static const ::std::filesystem::path              & CRH_HMMSCAN_DATA_DIR();

		static const ::std::filesystem::path              & NONEXISTENT_FILE();

		static const std::string                          & EXAMPLE_A_PDB_STEMNAME();
		static const std::string                          & EXAMPLE_B_PDB_STEMNAME();

		static const ::std::filesystem::path              & EXAMPLE_A_PDB_FILENAME();
		static const ::std::filesystem::path              & EXAMPLE_B_PDB_FILENAME();

		static const ::std::filesystem::path              & EXAMPLE_A_DSSP_FILENAME();
		static const ::std::filesystem::path              & EXAMPLE_B_DSSP_FILENAME();

		static const ::std::filesystem::path              & EXAMPLE_A_WOLF_FILENAME();
		static const ::std::filesystem::path              & EXAMPLE_B_WOLF_FILENAME();

		static const ::std::filesystem::path              & EXAMPLE_A_SEC_FILENAME();
		static const ::std::filesystem::path              & EXAMPLE_B_SEC_FILENAME();

		static const ::std::filesystem::path              & ALIGNMENT_FILE();
		static const ::std::filesystem::path              & TEST_PAIR_SUPPOSN_XML();

		static const ::std::filesystem::path              & MODIFIED_PDB_FILENAME();
		static const ::std::filesystem::path              & TEST_OUTPUT_FILENAME();

		static const ::std::filesystem::path              & EXAMPLE_DOUBLES_FILENAME();
		static const ::std::filesystem::path              & EXAMPLE_TUPLES_FILENAME();

		static const ::std::filesystem::path              & CRH_EG_DOMTBL_IN_FILENAME();
		static const ::std::filesystem::path              & CRH_EG_DOMTBL_OUT_FILENAME();
		static const ::std::filesystem::path              & CRH_EG_DOMTBL_LIMIT_2_OUT_FILENAME();
		static const ::std::filesystem::path              & CRH_EG_DOMTBL_JSON_OUT_FILENAME();
		static const ::std::filesystem::path              & CRH_EG_HMMSEARCH_IN_FILENAME();
		static const ::std::filesystem::path              & CRH_EG_HMMSEARCH_OUT_FILENAME();
		static const ::std::filesystem::path              & CRH_EG_HMMSEARCH_LIMIT_2_OUT_FILENAME();
		static const ::std::filesystem::path              & CRH_EG_HMMSEARCH_BIG_GAP_OUT_FILENAME();
		static const ::std::filesystem::path              & CRH_EG_HMMSEARCH_SMALL_GAP_OUT_FILENAME();
		static const ::std::filesystem::path              & CRH_EG_HMMSEARCH_TRIMMED_OUT_FILENAME();
		static const ::std::filesystem::path              & CRH_EG_HMMSEARCH_BIG_TRIM_OUT_FILENAME();
		static const ::std::filesystem::path              & CRH_EG_HMMSEARCH_HMMSEARCH_ALN_OUT_FILENAME();
		static const ::std::filesystem::path              & CRH_EG_HMMSEARCH_SUMMARISE_OUT_FILENAME();
		static const ::std::filesystem::path              & CRH_EG_HMMSEARCH_HTML_OUT_FILENAME();
		static const ::std::filesystem::path              & CRH_EG_HMMSEARCH_JSON_OUT_FILENAME();
		static const ::std::filesystem::path              & CRH_EG_RAW_EVALUE_IN_FILENAME();
		static const ::std::filesystem::path              & CRH_EG_RAW_EVALUE_OUT_FILENAME();
		static const ::std::filesystem::path              & CRH_EG_RAW_EVALUE_LIMIT_2_OUT_FILENAME();
		static const ::std::filesystem::path              & CRH_EG_RAW_SCORE_IN_FILENAME();
		static const ::std::filesystem::path              & CRH_EG_RAW_SCORE_OUT_FILENAME();
		static const ::std::filesystem::path              & CRH_EG_RAW_SCORE_LIMIT_2_OUT_FILENAME();

		static const double                               & LOOSER_ACCURACY_PERCENTAGE();
		static const double                               & ACCURACY_PERCENTAGE();

		template <typename T>
		static const T                                    & LOOSER_ACCURACY_PERCENTAGE_TMPL();
		template <typename T>
		static const T                                    & ACCURACY_PERCENTAGE_TMPL();

		static constexpr double                             DOUBLE_INFINITY      = std::numeric_limits<double>::infinity();
		static constexpr double                             DOUBLE_QUIET_NAN     = std::numeric_limits<double>::quiet_NaN();
		static constexpr double                             DOUBLE_SIGNALING_NAN = std::numeric_limits<double>::signaling_NaN();
		static constexpr double                             DOUBLE_DENORM_MIN    = std::numeric_limits<double>::denorm_min();
		static constexpr double                             DOUBLE_MIN           = std::numeric_limits<double>::min();
		static constexpr double                             DOUBLE_MAX           = std::numeric_limits<double>::max();

		static const doub_vec                             & INVALID_DOUBLES();

		static const path_vec & REQUIRED_PREEXISTING_FILES();
		static const path_vec & TEMPORARY_FILES();
	};

	template <> const float       & global_test_constants::LOOSER_ACCURACY_PERCENTAGE_TMPL<float>();
	template <> const double      & global_test_constants::LOOSER_ACCURACY_PERCENTAGE_TMPL<double>();
	template <> const long double & global_test_constants::LOOSER_ACCURACY_PERCENTAGE_TMPL<long double>();
	template <> const float       & global_test_constants::ACCURACY_PERCENTAGE_TMPL<float>();
	template <> const double      & global_test_constants::ACCURACY_PERCENTAGE_TMPL<double>();
	template <> const long double & global_test_constants::ACCURACY_PERCENTAGE_TMPL<long double>();
} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_TEST_CATH_TEST_GLOBAL_TEST_CONSTANTS_HPP
