/// \file
/// \brief The cath_superpose_program_exception_wrapper definitions

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

#include <boost/log/trivial.hpp>

#include "cath/cath_superpose/cath_superposer.hpp"
#include "cath/cath_superpose/options/cath_superpose_options.hpp"
#include "cath/chopping/domain/domain.hpp"
#include "cath/common/program_exception_wrapper.hpp"

using namespace ::cath;
using namespace ::cath::common;
using namespace ::cath::opts;
using namespace ::std;

using ::boost::log::trivial::info;
using ::boost::log::trivial::severity;

namespace {

	/// \brief A concrete program_exception_wrapper that implements do_run_program() to parse the options and then pass them to cath_superposer::superpose()
	///
	/// Using program_exception_wrapper allows the program to be wrapped in standard last-chance exception handling.
	class cath_superpose_program_exception_wrapper final : public program_exception_wrapper {
		string do_get_program_name() const final {
			return "cath-superpose";
		}

		/// \brief Parse the options and then pass them to cath_superposer::superpose()
		void do_run_program(int argc, char * argv[]) final {
			const auto the_cath_superpose_options = make_and_parse_options<cath_superpose_options>( argc, argv );
			get_sink_ptr()->set_filter( severity >= info );
			cath_superposer::superpose(the_cath_superpose_options);
		}
	};

} // namespace

/// \brief A main function for cath_superpose that just calls run_program() on a cath_superpose_program_exception_wrapper
int main(int argc, char * argv[] ) {
	return cath_superpose_program_exception_wrapper().run_program( argc, argv );
}
