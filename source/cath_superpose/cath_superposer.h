/// \file
/// \brief The cath_superposer class header

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

#ifndef CATH_SUPERPOSER_H_INCLUDED
#define CATH_SUPERPOSER_H_INCLUDED

#include <iostream>

namespace cath { namespace opts { class cath_superpose_options; } }
namespace cath { namespace sup { class superposition_context; } }

namespace cath {

	/// \brief A class containing a public static method to do the real work of a call to the cath-superpose command
	///
	/// The details of the job to be done are passed in a cath_superpose_options object.
	class cath_superposer {
	private:
		cath_superposer() = delete;

	public:
		static void superpose(const opts::cath_superpose_options &,
		                      std::istream &arg_istream = std::cin,
		                      std::ostream &arg_stdout = std::cout,
		                      std::ostream &arg_stderr = std::cerr);

		static sup::superposition_context get_superposition_context(const opts::cath_superpose_options &,
		                                                            std::istream &,
		                                                            std::ostream &);
	};

}

#endif
