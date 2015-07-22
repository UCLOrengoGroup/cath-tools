/// \file
/// \brief The entry_querier_dyn_prog_score_source class header

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

#ifndef ENTRY_QUERIER_DYN_PROG_SCORE_SOURCE_H_INCLUDED
#define ENTRY_QUERIER_DYN_PROG_SCORE_SOURCE_H_INCLUDED

#include "alignment/dyn_prog_align/dyn_prog_score_source/dyn_prog_score_source.h"

namespace cath {
	class entry_querier;
	class protein;

	namespace align {
		/// \brief TODOCUMENT
		class entry_querier_dyn_prog_score_source final : public dyn_prog_score_source {
		private:
			/// \brief TODOCUMENT
			const entry_querier &the_entry_querier;

			/// \brief TODOCUMENT
			const protein       &protein_a;

			/// \brief TODOCUMENT
			const protein       &protein_b;

			/// \brief TODOCUMENT
			const size_t         view_from_index_a;

			/// \brief TODOCUMENT
			const size_t         view_from_index_b;

			virtual size_t do_get_length_a() const override final;
			virtual size_t do_get_length_b() const override final;
			virtual score_type do_get_score(const size_t &,
			                                const size_t &) const override final;

		public:
			entry_querier_dyn_prog_score_source(const entry_querier &,
			                                    const protein &,
			                                    const protein &,
			                                    const size_t &,
			                                    const size_t &);
			virtual ~entry_querier_dyn_prog_score_source() noexcept = default;
		};
	}
}
#endif