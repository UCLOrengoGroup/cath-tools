/// \file
/// \brief The html_align_outputter class header

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

#ifndef HTML_ALIGN_OUTPUTTER_H_INCLUDED
#define HTML_ALIGN_OUTPUTTER_H_INCLUDED

#include "common/type_aliases.h"

#include <iosfwd>

namespace cath { namespace align { class alignment; } }
namespace cath { class display_colourer; }
namespace cath { namespace file { class pdb_list; } }

namespace cath {
	namespace align {

		/// \brief Simple wrapper class for outputting an alignment to an ostream in horizontal format
		///
		/// Use like this:
		///
		///    cerr << html_align_outputter( my_alignment, pdbs, names, colourer ) << endl;
		///
		/// This provides a convenient way for the user to choose a different format
		/// when outputting an alignment to an ostream via the insertion operator
		class html_align_outputter final {
		private:
			/// \brief A const-reference to the alignment to be output
			const alignment        &the_alignment;

			/// \brief TODOCUMENT
			const file::pdb_list   &pdbs;

			/// \brief TODOCUMENT
			const str_vec          &names;

			/// \brief TODOCUMENT
			const display_colourer &colourer;

		public:
			html_align_outputter(const alignment &,
			                     const file::pdb_list &,
			                     const str_vec &,
			                     const display_colourer &);

			const alignment & get_alignment() const;
			const file::pdb_list & get_pdbs() const;
			const str_vec & get_names() const;
			const display_colourer & get_display_colourer() const;
		};

		std::ostream & operator<<(std::ostream &,
		                          const html_align_outputter &);

	}
}

#endif
