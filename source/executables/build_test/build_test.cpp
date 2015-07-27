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
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_monitor.hpp>

#include <iostream>

using namespace boost::log;
using namespace boost::log::trivial;
using namespace boost::unit_test;
using namespace std;

namespace cath {
	namespace test {

		/// \brief TODOCUMENT
		void boost_exception_translator(const boost::exception &arg_excptn ///< TODOCUMENT
		                                ) {
			cerr << "The execution_monitor caught a boost::exception and passed it to the boost_exception_translator :\n"
			     << diagnostic_information( arg_excptn ) << "\n"
			     << endl;
			throw;
		}

		/// \brief TODOCUMENT
		class prepare_for_test_global_fixture final {
		public:
			prepare_for_test_global_fixture();
		};

		/// \brief TODOCUMENT
		prepare_for_test_global_fixture::prepare_for_test_global_fixture() {
			boost::log::core::get()->set_filter(
				severity >= info
//				severity >= trace
			);
			unit_test_monitor.register_exception_translator<boost::exception>( &boost_exception_translator );
		}

		/// \brief TODOCUMENT
		BOOST_GLOBAL_FIXTURE( prepare_for_test_global_fixture )

	}
}
