/// \file
/// \brief The residue_querier class header

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

#ifndef RESIDUE_QUERIER_H_INCLUDED
#define RESIDUE_QUERIER_H_INCLUDED

#include "structure/entry_querier/entry_querier.h"

namespace cath {

	/// \brief TODOCUMENT
	class residue_querier final : public entry_querier {
	private:
		virtual size_t           do_get_length(const cath::protein &) const override final;
		virtual double           do_get_gap_penalty_ratio() const override final;
		virtual size_t           do_num_excluded_on_either_size() const override final;
		virtual float_score_type do_optimum_single_score() const override final;

		virtual std::string  do_get_entry_name() const override final;
		virtual score_type   do_distance_score__offset_1(const cath::protein &,
		                                                 const cath::protein &,
		                                                 const size_t &,
		                                                 const size_t &,
		                                                 const size_t &,
		                                                 const size_t &) const override final;

		virtual bool         do_are_comparable__offset_1(const cath::protein &,
		                                                 const cath::protein &,
		                                                 const size_t &,
		                                                 const size_t &,
		                                                 const size_t &,
		                                                 const size_t &) const override final;

		virtual bool         do_are_similar__offset_1(const cath::protein &,
		                                              const cath::protein &,
		                                              const size_t &,
		                                              const size_t &) const override final;

		virtual bool         do_temp_hacky_is_residue() const override final;

	public:
		virtual ~residue_querier() noexcept = default;

		static const float_score_type RESIDUE_A_VALUE;
		static const float_score_type RESIDUE_B_VALUE;
		static const float_score_type RESIDUE_MIN_SCORE_CUTOFF;
		static const float_score_type RESIDUE_MAX_DIST_SQ_CUTOFF;
	};

}

#endif
