/// \file
/// \brief The horiz_align_outputter class header

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

#ifndef CATH_TOOLS_SOURCE_CT_UNI_CATH_ALIGNMENT_IO_OUTPUTTER_HORIZ_ALIGN_OUTPUTTER_HPP
#define CATH_TOOLS_SOURCE_CT_UNI_CATH_ALIGNMENT_IO_OUTPUTTER_HORIZ_ALIGN_OUTPUTTER_HPP

#include <iosfwd>

// clang-format off
namespace cath::align { class alignment; }
// clang-format on

namespace cath::align {

	/// \brief Simple wrapper class for outputting an alignment to an ostream in horizontal format
	///
	/// Use like this:
	///
	///    cerr << horiz_align_outputter( my_alignment ) << endl;
	///
	/// This provides a convenient way for the user to choose a different format
	/// when outputting an alignment to an ostream via the insertion operator
	///
	/// \todo Consider writing a human_readable_align_outputter that could display:
	///        * pair scores with extended ASCII codes ' ', '░', '▒', '▓' and '█'
	///        * secondary structures with '\' for helix, '=' for strand and '—' for coil, eg: '———————=====————\\\\\\\\—————'
	class horiz_align_outputter final {
	private:
		/// \brief A const-reference to the alignment to be output
		const alignment &the_alignment;

	public:
		explicit horiz_align_outputter(const alignment &);

		[[nodiscard]] const alignment &get_alignment() const;
	};

	std::ostream & operator<<(std::ostream &,
	                          const horiz_align_outputter &);

} // namespace cath::align

#endif // CATH_TOOLS_SOURCE_CT_UNI_CATH_ALIGNMENT_IO_OUTPUTTER_HORIZ_ALIGN_OUTPUTTER_HPP
