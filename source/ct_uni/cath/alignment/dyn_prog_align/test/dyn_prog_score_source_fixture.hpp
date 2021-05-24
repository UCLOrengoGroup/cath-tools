/// \file
/// \brief The prototypes for functions that generate examples of each of the dyn_prog_score_source

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

#ifndef _CATH_TOOLS_SOURCE_CT_UNI_CATH_ALIGNMENT_DYN_PROG_ALIGN_TEST_DYN_PROG_SCORE_SOURCE_FIXTURE_HPP
#define _CATH_TOOLS_SOURCE_CT_UNI_CATH_ALIGNMENT_DYN_PROG_ALIGN_TEST_DYN_PROG_SCORE_SOURCE_FIXTURE_HPP

#include <string_view>

#include "cath/alignment/dyn_prog_align/dyn_prog_score_source/mask_dyn_prog_score_source.hpp"
#include "cath/alignment/dyn_prog_align/dyn_prog_score_source/new_matrix_dyn_prog_score_source.hpp"
#include "cath/alignment/dyn_prog_align/dyn_prog_score_source/old_matrix_dyn_prog_score_source.hpp"
#include "cath/alignment/dyn_prog_align/dyn_prog_score_source/sequence_string_dyn_prog_score_source.hpp"
#include "cath/common/container/vector_of_vector.hpp"
#include "cath/test/global_test_constants.hpp"

namespace cath::align {

	/// \brief A test fixture for dyn_prog_score_source tests that adds some extras to global_test_constants
	class dyn_prog_score_source_fixture : protected global_test_constants {
	private:
//		const align::entry_querier_dyn_prog_score_source &   make_example_entry_querier_dyn_prog_score_source();
		const align::mask_dyn_prog_score_source &            make_example_mask_dyn_prog_score_source();
		const align::old_matrix_dyn_prog_score_source &      make_example_old_matrix_dyn_prog_score_source();
		const align::new_matrix_dyn_prog_score_source &      make_example_new_matrix_dyn_prog_score_source();
		const align::sequence_string_dyn_prog_score_source & make_example_sequence_string_dyn_prog_score_source();

		/// \brief A first  sequence string for making an example sequence_string_dyn_prog_score_source
		static constexpr ::std::string_view SEQUENCE_STRING_A = "CC";

		/// \brief A second sequence string for making an example sequence_string_dyn_prog_score_source
		static constexpr ::std::string_view SEQUENCE_STRING_B = "ACCD";

	public:
		template <typename DPSS>
		DPSS make_example_dyn_prog_score_source();

	};

	/// \brief A template function for constructing an example of each type of dyn_prog_score_source
	///
	/// This default just returns a default constructed instance of the dyn_prog_score_source
	/// but this must be specialised for any dyn_prog_score_source without a default ctor
	template <typename DPSS>
	DPSS dyn_prog_score_source_fixture::make_example_dyn_prog_score_source() {
		return DPSS();
	}

//	template <>
//	align::entry_querier_dyn_prog_score_source dyn_prog_score_source_fixture::make_example_dyn_prog_score_source<align::entry_querier_dyn_prog_score_source>();

	template <>
	align::mask_dyn_prog_score_source dyn_prog_score_source_fixture::make_example_dyn_prog_score_source<align::mask_dyn_prog_score_source>();

	template <>
	align::old_matrix_dyn_prog_score_source dyn_prog_score_source_fixture::make_example_dyn_prog_score_source<align::old_matrix_dyn_prog_score_source>();

	template <>
	align::new_matrix_dyn_prog_score_source dyn_prog_score_source_fixture::make_example_dyn_prog_score_source<align::new_matrix_dyn_prog_score_source>();

	template <>
	align::sequence_string_dyn_prog_score_source dyn_prog_score_source_fixture::make_example_dyn_prog_score_source<align::sequence_string_dyn_prog_score_source>();

} // namespace cath::align

#endif // _CATH_TOOLS_SOURCE_CT_UNI_CATH_ALIGNMENT_DYN_PROG_ALIGN_TEST_DYN_PROG_SCORE_SOURCE_FIXTURE_HPP
