/// \file
/// \brief The cath_score_align_program_exception_wrapper definitions

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

#include "cath_score_align/cath_align_scorer.h"
#include "cath_score_align/options/cath_score_align_options.h"
#include "common/program_exception_wrapper.h"

using namespace cath::common;
using namespace cath::opts;
using namespace std;

namespace cath {

	/// \brief A concrete program_exception_wrapper that implements do_run_program() to parse the options and then pass them to cath_align_scorer::score()
	///
	/// Using program_exception_wrapper allows the program to be wrapped in standard last-chance exception handling.
	class cath_score_align_program_exception_wrapper final : public program_exception_wrapper {
		virtual string do_get_program_name() const override final {
			return "cath-score-align";
		}

		/// \brief Parse the options and then pass them to cath_score_alignr::superpose()
		virtual void do_run_program(int argc, char * argv[]) override final {
			const auto the_cath_score_align_options = make_and_parse_options<cath_score_align_options>( argc, argv );
			cath_align_scorer::score(the_cath_score_align_options);
		}
	};
} // namespace cath

/// \brief A main function for cath_score_align that just calls run_program() on a cath_score_align_program_exception_wrapper
int main(int argc, char * argv[] ) {
	return cath::cath_score_align_program_exception_wrapper().run_program( argc, argv );
}
