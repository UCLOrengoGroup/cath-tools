/// \file
/// \brief The argc_argv_faker class header

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

#ifndef _CATH_TOOLS_SOURCE_SRC_COMMON_COMMON_ARGC_ARGV_FAKER_HPP
#define _CATH_TOOLS_SOURCE_SRC_COMMON_COMMON_ARGC_ARGV_FAKER_HPP

#include <boost/shared_array.hpp>

#include "common/type_aliases.hpp"

#include <iosfwd>
#include <vector>

namespace cath {

	/// \brief A simple class to allow the faking of argc and argv values from a vector of strings
	///        for those times you need to access an old-school interface directly.
	///
	/// This can also be useful for taking a local copy of const argv style arguments so that they
	/// that can be passed to a parser that doesn't respect that const
	///
	/// Of course, the argc_argv_faker must persist at least as long as the argv result is used
	/// because the data pointed to by argv will be deleted along with the argc_argv_faker.
	///
	/// Example usage:
	///
	/// \code
	/// const str_vec my_args = { "my", "argument", "list" };
	/// argc_argv_faker my_argc_argv_faker( my_args );
	/// GetOpt(my_argc_argv_faker.get_argc(), my_argc_argv_faker.get_argv(), ....
	/// \endcode
	class argc_argv_faker final {
	private:
		/// \brief TODOCUMENT
		std::vector<boost::shared_array<char> > arguments;

		/// \brief TODOCUMENT
		std::vector<char *                    > argument_ptrs;

		/// \brief TODOCUMENT
		int argc;

		void init(const str_vec &);

	public:
		explicit argc_argv_faker(const str_vec &);
		argc_argv_faker(const int &,
		                const char * const []);

		int & get_argc();
		const int & get_argc() const;
		char * * get_argv();
		char * const * get_argv() const;
	};

	std::ostream & operator<<(std::ostream &,
	                          const argc_argv_faker &);

} // namespace cath

#endif
