/// \file
/// \brief The sec_struc_querier class header

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

#ifndef _CATH_TOOLS_SOURCE_STRUCTURE_ENTRY_QUERIER_SEC_STRUC_QUERIER_H
#define _CATH_TOOLS_SOURCE_STRUCTURE_ENTRY_QUERIER_SEC_STRUC_QUERIER_H

#include "structure/entry_querier/entry_querier.h"

namespace cath {

	/// \brief TODOCUMENT.
	///
	///
	class sec_struc_querier final : public entry_querier {
	private:
		virtual size_t       do_get_length(const cath::protein &) const override final;
		virtual double       do_get_gap_penalty_ratio() const override final;
		virtual size_t       do_num_excluded_on_either_size() const override final;
		virtual double       do_optimum_single_score() const override final;

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
		/// \brief The value a used in the SSAP paper (for secondary structures)
		///
		/// Note that this scaled by INTEGER_SCALING^2 ( = 10 * 10 = 100), which matches
		/// the result of the distance being scaled by INTEGER_SCALING
		static constexpr size_t SEC_STRUC_A_VALUE = 800 * entry_querier::INTEGER_SCALING * entry_querier::INTEGER_SCALING;

		/// \brief The value b used in the SSAP paper (for secondary structures)
		///
		/// Note that this scaled by INTEGER_SCALING^2 ( = 10 * 10 = 100), which matches
		/// the result of the distance being scaled by INTEGER_SCALING
		static constexpr size_t SEC_STRUC_B_VALUE=   6 * entry_querier::INTEGER_SCALING * entry_querier::INTEGER_SCALING;

		/// \brief The minimum score that the algorithm will consider for secondary structures
		///        (often referred to as c)
		///
		/// Note that this does not need scaling
		static constexpr size_t SEC_STRUC_MIN_SCORE_CUTOFF = 10;

		/// \brief The maximum squared-distance value that the algorithm will consider (for secondary structures)
		///
		/// This is calculated as \f$ \frac{a}{c} - b\f$,
		/// which comes from solving \f$ c = \frac{a}{ d^2 + b} \f$ for \f$ s^2 \f$ )
		///
		/// NOTE: This has an implied scaling ( * INTEGER_SCALING * INTEGER_SCALING) built into it.
		///       This line will not need changing as long as A, B and the distances have their scaling
		///       removed simultaneously.
		///
		/// With current (unscaled) values of a=800, b=6 and c=10, this gives a maximum (unscaled) value
		/// of 74 for d^2, which implies a maximum value for d of sqrt(74) ~= 8.6A
		static constexpr size_t SEC_STRUC_MAX_DIST_SQ_CUTOFF = SEC_STRUC_A_VALUE / SEC_STRUC_MIN_SCORE_CUTOFF - SEC_STRUC_B_VALUE;
	};

} // namespace cath

#endif
