/// \file
/// \brief The prototypes for functions that generate examples of each of the alignment

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

#ifndef _CATH_TOOLS_SOURCE_CT_UNI_CATH_ALIGNMENT_TEST_ALIGNMENT_FIXTURE_HPP
#define _CATH_TOOLS_SOURCE_CT_UNI_CATH_ALIGNMENT_TEST_ALIGNMENT_FIXTURE_HPP

#include <optional>

#include "cath/alignment/alignment.hpp"
#include "cath/common/cpp20/make_array.hpp"
#include "cath/test/global_test_constants.hpp"

// Add a template factory function, make_example_alignment(), that is specialised
// for constructing each of these types.

namespace cath {
	namespace align {

		/// \brief A test fixture for alignment tests that adds some extras to global_test_constants
		class alignment_fixture {
		  protected:
			// A normal list
			static constexpr auto aln_list_a = common::make_array<aln_posn_opt>( 0, 1, 2, 3 );

			// A normal list
			static constexpr auto aln_list_b = common::make_array<aln_posn_opt>( 0, 1, 2, ::std::nullopt );

			// A list that is longer than the normal lists
			static constexpr auto aln_list_long = common::make_array<aln_posn_opt>( 0, 1, 2, 3, 4, 5 );

			static constexpr auto example_scores = common::make_array<score_opt>( 3.6, 6.8, 2.1, 999.999 );

			static alignment aln_a_a();
			static alignment aln_a_b();
			static alignment aln_b_a();
			static alignment aln_long_long();
		};

	} // namespace align
} // namespace cath
#endif // _CATH_TOOLS_SOURCE_CT_UNI_CATH_ALIGNMENT_TEST_ALIGNMENT_FIXTURE_HPP
