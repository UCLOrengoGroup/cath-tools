/// \file
/// \brief The dyn_prog_aligner class header

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

#ifndef _CATH_TOOLS_SOURCE_ALIGNMENT_DYN_PROG_ALIGN_DYN_PROG_ALIGNER_H
#define _CATH_TOOLS_SOURCE_ALIGNMENT_DYN_PROG_ALIGN_DYN_PROG_ALIGNER_H

#include "alignment/align_type_aliases.hpp"
#include "common/type_aliases.hpp"
#include "ssap/windowed_matrix.hpp"

#include <memory>

namespace cath { namespace align { class dyn_prog_score_source; } }
namespace cath { namespace align { namespace gap { class gap_penalty; } } }

namespace cath {
	namespace align {

		/// \brief ABC defining interface for classes that align generic dyn_prog_score_sources with dynamic-programming
		class dyn_prog_aligner {
		public:
			/// \brief TODOCUMENT
			using size_type = windowed_matrix::size_type;

		private:
			/// \brief A standard do_clone() method to act as a virtual copy-ctor
			///
			/// This is a pure virtual function (so must be overridden by any concrete, derived classes).
			virtual std::unique_ptr<dyn_prog_aligner> do_clone() const = 0;

			/// \brief TODOCUMENT
			virtual score_alignment_pair do_align(const dyn_prog_score_source &,
			                                      const gap::gap_penalty &,
			                                      const size_type &) const = 0;

		public:
			dyn_prog_aligner() = default;
			std::unique_ptr<dyn_prog_aligner> clone() const;
			virtual ~dyn_prog_aligner() noexcept = default;

			dyn_prog_aligner(const dyn_prog_aligner &) = default;
			dyn_prog_aligner(dyn_prog_aligner &&) noexcept = default;
			dyn_prog_aligner & operator=(const dyn_prog_aligner &) = default;
			dyn_prog_aligner & operator=(dyn_prog_aligner &&) noexcept = default;

			score_alignment_pair align(const dyn_prog_score_source &,
			                           const gap::gap_penalty &,
			                           const size_type &) const;
		};

	} // namespace align
} // namespace cath

#endif
