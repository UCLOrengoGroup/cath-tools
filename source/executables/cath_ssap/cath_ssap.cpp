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

#include <spdlog/spdlog.h>

#include "cath/chopping/domain/domain.hpp"
#include "cath/common/program_exception_wrapper.hpp"
#include "cath/ssap/options/cath_ssap_options.hpp"
#include "cath/ssap/ssap.hpp"

using namespace ::cath;
using namespace ::cath::common;
using namespace ::cath::opts;
using namespace ::std;

namespace {

	/// \brief A concrete program_exception_wrapper that implements do_run_program() to parse the options and then pass them to call_ssap()
	///
	/// Using program_exception_wrapper allows the program to be wrapped in standard last-chance exception handling.
	class ssap_program_exception_wrapper final : public program_exception_wrapper {
		[[nodiscard]] string_view do_get_program_name() const final {
			return "cath-ssap";
		}

		/// \brief Parse the options and then pass them to run_ssap()
		void do_run_program(int argc, char * argv[]) final {
			const auto the_cath_ssap_options = make_and_parse_options<cath_ssap_options>( argc, argv );

			// Set log level to: trace if debug requested, or info otherwise
			::spdlog::set_level( the_cath_ssap_options.get_old_ssap_options().get_debug() ? ::spdlog::level::trace
			                                                                              : ::spdlog::level::info );

			run_ssap( the_cath_ssap_options );
		}
	};

} // namespace

/// \brief A main function for SSAP that just calls run_program() on a ssap_program_exception_wrapper
int main(int argc, char * argv[]) {
	return ssap_program_exception_wrapper().run_program( argc, argv );
}
