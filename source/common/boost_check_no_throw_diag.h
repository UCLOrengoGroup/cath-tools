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

#ifndef BOOST_CHECK_NO_THROW_DIAG_H_INCLUDED
#define BOOST_CHECK_NO_THROW_DIAG_H_INCLUDED

#define BOOST_CHECK_NO_THROW_DIAG_IMPL( S, TL )                                                          \
    try {                                                                                           \
        S;                                                                                          \
        BOOST_CHECK_IMPL( true, "no exceptions thrown by " BOOST_STRINGIZE( S ), TL, CHECK_MSG ); } \
    catch( ... ) {                                                                                  \
        BOOST_CHECK_IMPL( false, "exception thrown by " BOOST_STRINGIZE( S ), TL, CHECK_MSG );      \
        throw;                                                                                      \
    }                                                                                               \

#define BOOST_WARN_NO_THROW_DIAG( S )            BOOST_CHECK_NO_THROW_DIAG_IMPL( S, WARN )
#define BOOST_CHECK_NO_THROW_DIAG( S )           BOOST_CHECK_NO_THROW_DIAG_IMPL( S, CHECK )
#define BOOST_REQUIRE_NO_THROW_DIAG( S )         BOOST_CHECK_NO_THROW_DIAG_IMPL( S, REQUIRE )

#endif
