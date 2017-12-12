/// \file
/// \brief The gen_dyn_prog_string_aligner class header

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

#ifndef _CATH_TOOLS_SOURCE_UNI_ALIGNMENT_DYN_PROG_ALIGN_DETAIL_STRING_ALIGNER_GEN_DYN_PROG_STRING_ALIGNER_H
#define _CATH_TOOLS_SOURCE_UNI_ALIGNMENT_DYN_PROG_ALIGN_DETAIL_STRING_ALIGNER_GEN_DYN_PROG_STRING_ALIGNER_H

#include "alignment/dyn_prog_align/detail/string_aligner/string_aligner.hpp"

#include <memory>

namespace cath {
	namespace align {
		class dyn_prog_aligner;

		namespace detail {

			/// \brief A string_aligner that uses a (more general) dyn_prog_aligner to perform the string alignment
			///
			/// The gen_dyn_prog_string_aligner takes a dyn_prog_aligner argument when it's being constructed
			/// and stores a clone() of it.
			class gen_dyn_prog_string_aligner final : public string_aligner {
			private:
				/// \brief A smart pointer to hold the clone() of the dyn_prog_aligner passed to the ctor
				const std::unique_ptr<dyn_prog_aligner> dyn_prog_aligner_ptr;

				void check_dyn_prog_aligner_ptr() const;
				dyn_prog_aligner & get_dyn_prog_aligner();
				const dyn_prog_aligner & get_dyn_prog_aligner() const;

				str_str_pair do_align(const std::string &,
				                      const std::string &,
				                      const gap::gap_penalty &) const final;

			public:
				explicit gen_dyn_prog_string_aligner(const dyn_prog_aligner &);
			};

		} // namespace detail
	} // namespace align
} // namespace cath

#endif
