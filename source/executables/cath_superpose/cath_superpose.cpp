/// \file
/// \brief The cath_superpose_program_exception_wrapper definitions

/// \copyright
/// CATH Binaries - Protein structure comparison tools such as SSAP and SNAP
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

#include "cath_superpose/cath_superposer.h"
#include "exception/program_exception_wrapper.h"
#include "options/executable/cath_superpose_options/cath_superpose_options.h"

using namespace cath::common;
using namespace cath::opts;
using namespace std;

namespace cath {

	/// \brief A concrete program_exception_wrapper that implements do_run_program() to parse the options and then pass them to cath_superposer::superpose()
	///
	/// Using program_exception_wrapper allows the program to be wrapped in standard last-chance exception handling.
	class cath_superpose_program_exception_wrapper final : public program_exception_wrapper {
		virtual string do_get_program_name() const override final {
			return "cath-superpose";
		}

		/// \brief Parse the options and then pass them to cath_superposer::superpose()
		virtual void do_run_program(int argc, char * argv[]) override final {
			cath_superpose_options the_cath_superpose_options;
			the_cath_superpose_options.parse_options(argc, argv);
			cath_superposer::superpose(the_cath_superpose_options);
		}
	};
}

/// \brief A main function for cath_superpose that just calls run_program() on a cath_superpose_program_exception_wrapper
int main(int argc, char * argv[] ) {
	return cath::cath_superpose_program_exception_wrapper().run_program( argc, argv );
}