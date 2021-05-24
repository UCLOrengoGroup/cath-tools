/// \file
/// \brief The boost_check_no_throw_diag header

/// \copyright
/// Tony Lewis's Common C++ Library Code (here imported into the CATH Tools project and then tweaked, eg namespaced in cath)
/// Copyright (C) 2007, Tony Lewis
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
///
/// This provides an extra testing tool that's like BOOST_CHECK_NO_THROW()
/// but that re-throws any exception thrown by the tested code after it has
/// registered the failure. This allows the exception to be caught higher up and
/// handled by any exception translators.
///
/// This is because I think BOOST_CHECK_NO_THROW() wastes time when there
/// is a failure because it tells the user nothing about the exception that
/// was thrown.

#ifndef _CATH_TOOLS_SOURCE_CT_TEST_CATH_TEST_BOOST_ADDENDA_BOOST_CHECK_NO_THROW_DIAG_HPP
#define _CATH_TOOLS_SOURCE_CT_TEST_CATH_TEST_BOOST_ADDENDA_BOOST_CHECK_NO_THROW_DIAG_HPP

#include <boost/test/tools/detail/fwd.hpp>
#include <boost/test/tools/interface.hpp>

namespace cath::test {

	#define BOOST_CHECK_NO_THROW_DIAG_IMPL( S, TL )                                                      \
	    try {                                                                                            \
	        BOOST_TEST_PASSPOINT();                                                                      \
	        S;                                                                                           \
	        BOOST_TEST_TOOL_DIRECT_IMPL( true,  TL, "no exceptions thrown by " BOOST_STRINGIZE( S ) );   \
	    }                                                                                                \
	    catch( ... ) {                                                                                   \
	        BOOST_TEST_TOOL_DIRECT_IMPL( false, TL, "exception thrown by "     BOOST_STRINGIZE( S ) );   \
	        throw;                                                                                       \
	    }                                                                                                \
	    /**/

	#define BOOST_WARN_NO_THROW_DIAG( S )            BOOST_CHECK_NO_THROW_DIAG_IMPL( S, WARN    )
	#define BOOST_CHECK_NO_THROW_DIAG( S )           BOOST_CHECK_NO_THROW_DIAG_IMPL( S, CHECK   )
	#define BOOST_REQUIRE_NO_THROW_DIAG( S )         BOOST_CHECK_NO_THROW_DIAG_IMPL( S, REQUIRE )

} // namespace cath::test

#endif // _CATH_TOOLS_SOURCE_CT_TEST_CATH_TEST_BOOST_ADDENDA_BOOST_CHECK_NO_THROW_DIAG_HPP
