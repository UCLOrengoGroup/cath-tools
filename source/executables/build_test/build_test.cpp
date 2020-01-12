/// \file
/// \brief The main test file that defines the test module "Cath Tools Master Test Suite"

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

#define BOOST_TEST_MODULE Cath Tools Master Test Suite
#define BOOST_AUTO_TEST_MAIN

#include <boost/exception/diagnostic_information.hpp>
#include <boost/log/core.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/trivial.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_monitor.hpp>

#include "common/test_or_exe_run_mode.hpp"
#include "test/predicate/bootstrap_mode.hpp"

#include <iostream>

using namespace cath::common;
using namespace cath::test;

using boost::log::trivial::info;
using boost::log::trivial::severity;
using boost::unit_test::unit_test_monitor;
using std::cerr;
using std::endl;

namespace {

	/// \brief TODOCUMENT
	void boost_exception_translator(const boost::exception &prm_excptn ///< TODOCUMENT
									) {
		cerr << "The execution_monitor caught a boost::exception and passed it to the boost_exception_translator :\n"
				<< diagnostic_information( prm_excptn ) << "\n"
				<< endl;
		throw;
	}

	/// \brief TODOCUMENT
	class prepare_for_test_global_fixture final {
	public:
		prepare_for_test_global_fixture();
		~prepare_for_test_global_fixture();

		prepare_for_test_global_fixture( const prepare_for_test_global_fixture & ) = delete;
		prepare_for_test_global_fixture( prepare_for_test_global_fixture && ) noexcept = delete;
		prepare_for_test_global_fixture & operator=( const prepare_for_test_global_fixture & ) = delete;
		prepare_for_test_global_fixture & operator=( prepare_for_test_global_fixture && ) noexcept = delete;

		static void warn_if_bootstrapping();
	};

	/// \brief TODOCUMENT
	prepare_for_test_global_fixture::prepare_for_test_global_fixture() {
		boost::log::core::get()->set_filter(
			severity >= info
//				severity >= trace
		);
		unit_test_monitor.register_exception_translator<boost::exception>( &boost_exception_translator );

		// Try to warn if in bootstrapping mode
		warn_if_bootstrapping();

		// Record that this is running in test mode
		run_mode_flag::value = run_mode::TEST;
	}

	/// \brief In the dtor, try to warn if in bootstrapping mode
	prepare_for_test_global_fixture::~prepare_for_test_global_fixture() {
		try {
			warn_if_bootstrapping();
		}
		catch (...) {
		}
	}

	/// \brief Print a conspicuous warning if bootstrapping is turned on
	void prepare_for_test_global_fixture::warn_if_bootstrapping() {
		if ( bootstrap_env_var_is_set() ) {
			BOOST_LOG_TRIVIAL( warning )
				<< "\033[1m Environment variable "
				<< get_bootstrap_env_var()
				<<" is set, so test files are being bootstrapped\033[0m";
		}
	}

} // namespace

/// \brief TODOCUMENT
BOOST_GLOBAL_FIXTURE( prepare_for_test_global_fixture );
