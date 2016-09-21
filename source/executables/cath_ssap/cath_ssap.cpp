/// \file
/// \brief The ssap_program_exception_wrapper definitions

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

#include <boost/log/expressions.hpp>
#include <boost/log/trivial.hpp>

#include "exception/program_exception_wrapper.h"
#include "options/executable/cath_ssap_options/cath_ssap_options.h"
#include "ssap/ssap.h"

using namespace boost::log::trivial;
using namespace cath::common;
using namespace cath::opts;
using namespace std;

namespace cath {

	/// \brief A concrete program_exception_wrapper that implements do_run_program() to parse the options and then pass them to call_ssap()
	///
	/// Using program_exception_wrapper allows the program to be wrapped in standard last-chance exception handling.
	class ssap_program_exception_wrapper final : public program_exception_wrapper {
		virtual string do_get_program_name() const override final {
			return "cath-ssap";
		}

		/// \brief Parse the options and then pass them to run_ssap()
		virtual void do_run_program(int argc, char * argv[]) override final {
			const auto the_cath_ssap_options = make_and_parse_options<cath_ssap_options>( argc, argv );

			// If using Boost Log, and if debug requested, set the default filter to : permit all messages of severity info or more
			if ( the_cath_ssap_options.get_old_ssap_options().get_debug() ) {
//				get_sink_ptr()->set_filter( severity >= debug );
				get_sink_ptr()->set_filter( severity >= trace );
			}

			run_ssap( the_cath_ssap_options );
		}
	};
}

/// \brief A main function for SSAP that just calls run_program() on a ssap_program_exception_wrapper
int main(int argc, char * argv[]) {
	return cath::ssap_program_exception_wrapper().run_program( argc, argv );
}
