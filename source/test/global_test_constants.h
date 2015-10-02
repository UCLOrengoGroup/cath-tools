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

#ifndef GLOBAL_TEST_CONSTANTS_H_INCLUDED
#define GLOBAL_TEST_CONSTANTS_H_INCLUDED

#include <boost/filesystem.hpp>

#include "structure/entry_querier/residue_querier.h"
#include "structure/entry_querier/sec_struc_querier.h"

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
		static const cath::residue_querier                & EXAMPLE_RESIDUE_QUERIER();
		static const cath::sec_struc_querier              & EXAMPLE_SEC_STRUC_QUERIER();

		static const boost::filesystem::path              & TEST_BASIC_FILE_TEST_DATA_DIR();
		static const boost::filesystem::path              & TEST_MULTI_SSAP_SUPERPOSE_DIR();
		static const boost::filesystem::path              & TEST_RESIDUE_NAMES_DATA_DIR();
		static const boost::filesystem::path              & TEST_SOURCE_DATA_DIR();
		static const boost::filesystem::path              & TEST_EXAMPLE_PDBS_DATA_DIR();
		static const boost::filesystem::path              & TEST_SSAP_REGRESSION_DATA_DIR();
		static const boost::filesystem::path              & TEST_SSAP_ALIGNMENT_GLUING_DATA_DIR();
		static const boost::filesystem::path              & TEST_SUP_JSON_DIR();
		static const boost::filesystem::path              & TEST_OUTPUT_DIRECTORY();

		static const boost::filesystem::path              & NONEXISTENT_FILE();

		static const std::string                          & EXAMPLE_A_PDB_STEMNAME();
		static const std::string                          & EXAMPLE_B_PDB_STEMNAME();

		static const boost::filesystem::path              & EXAMPLE_A_PDB_FILENAME();
		static const boost::filesystem::path              & EXAMPLE_B_PDB_FILENAME();

		static const boost::filesystem::path              & EXAMPLE_A_DSSP_FILENAME();
		static const boost::filesystem::path              & EXAMPLE_B_DSSP_FILENAME();

		static const boost::filesystem::path              & EXAMPLE_A_WOLF_FILENAME();
		static const boost::filesystem::path              & EXAMPLE_B_WOLF_FILENAME();

		static const boost::filesystem::path              & EXAMPLE_A_SEC_FILENAME();
		static const boost::filesystem::path              & EXAMPLE_B_SEC_FILENAME();

		static const boost::filesystem::path              & ALIGNMENT_FILE();
		static const boost::filesystem::path              & TEST_PAIR_SUPPOSN_XML();

		static const boost::filesystem::path              & MODIFIED_PDB_FILENAME();
		static const boost::filesystem::path              & TEST_OUTPUT_FILENAME();

		static const boost::filesystem::path              & EXAMPLE_DOUBLES_FILENAME();
		static const boost::filesystem::path              & EXAMPLE_TUPLES_FILENAME();

		static const double                               & LOOSER_ACCURACY_PERCENTAGE();
		static const double                               & ACCURACY_PERCENTAGE();

		template <typename T>
		static const T                                    & LOOSER_ACCURACY_PERCENTAGE_TMPL();
		template <typename T>
		static const T                                    & ACCURACY_PERCENTAGE_TMPL();

		static const double                               & DOUBLE_INFINITY();
		static const double                               & DOUBLE_QUIET_NAN();
		static const double                               & DOUBLE_SIGNALING_NAN();
		static const double                               & DOUBLE_DENORM_MIN();
		static const double                               & DOUBLE_MIN();
		static const double                               & DOUBLE_MAX();

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
}

#endif
