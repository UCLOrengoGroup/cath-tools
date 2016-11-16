/// \file
/// \brief The ssap_code_dyn_prog_aligner class header

/// \copyright
/// CATH Tools - Protein structure comparison tools such as SSAP and SNAP
/// Copyright (C) 1989, Orengo Group, University College London
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

#ifndef _CATH_TOOLS_SOURCE_ALIGNMENT_DYN_PROG_ALIGN_SSAP_CODE_DYN_PROG_ALIGNER_H
#define _CATH_TOOLS_SOURCE_ALIGNMENT_DYN_PROG_ALIGN_SSAP_CODE_DYN_PROG_ALIGNER_H

#include <boost/tuple/tuple.hpp>

#include "alignment/dyn_prog_align/dyn_prog_aligner.h"
#include "common/type_aliases.h"

namespace cath {
	namespace align {
		class alignment;
		class dyn_prog_score_source;

		/// \brief TODOCUMENT
		///
		/// WARNING: This uses static data so using separate instances is not enough
		///          to make this thread-safe.
		class ssap_code_dyn_prog_aligner final : public dyn_prog_aligner {
		private:
			virtual std::unique_ptr<dyn_prog_aligner> do_clone() const override final;

			using size_size_int_int_score_tuple = std::tuple<size_t,
			                                                 size_t,
			                                                 int,
			                                                 int,
			                                                 score_type>;

			static size_size_int_int_score_tuple score_matrix(const dyn_prog_score_source &,
			                                                  const score_type &,
			                                                  const size_type &,
			                                                  int_vec_vec &);

			static alignment traceback(const size_t &,
			                           const size_t &,
			                           const int &,
			                           const int &,
			                           const int &,
			                           const int &,
			                           const int_vec_vec &);

			static void traceback_recursive(alignment &,
			                                const int &,
			                                const int &,
			                                const int &,
			                                const int &,
			                                const int &,
			                                const int &,
			                                const int_vec_vec &);

			virtual score_alignment_pair do_align(const dyn_prog_score_source &,
			                                      const gap::gap_penalty &,
			                                      const size_type &) const override final;
		};
	} // namespace align
} // namespace cath

#endif
