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
///        and a template factory function, make_example_alignment(), that is specialised
///        for constructing each of these types.

#ifndef ALIGNMENT_FIXTURE_H_INCLUDED
#define ALIGNMENT_FIXTURE_H_INCLUDED

#include "alignment/alignment.h"
#include "test/global_test_constants.h"

namespace cath {
	namespace align {

		/// \brief A test fixture for alignment tests that adds some extras to global_test_constants
		class alignment_fixture {
		protected:
			// A normal list
			static const opt_aln_posn_vec aln_list_a;
			// A normal list
			static const opt_aln_posn_vec aln_list_b;
			// A list that is longer than the normal lists
			static const opt_aln_posn_vec aln_list_long;

			static const alignment aln_a_a;
			static const alignment aln_a_b;
			static const alignment aln_b_a;
			static const alignment aln_long_long;

			static const cath::opt_score_vec example_scores;
		};

	}
}
#endif
