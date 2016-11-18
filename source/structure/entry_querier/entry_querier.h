/// \file
/// \brief The entry_querier class header

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

#ifndef _CATH_TOOLS_SOURCE_STRUCTURE_ENTRY_QUERIER_ENTRY_QUERIER_H
#define _CATH_TOOLS_SOURCE_STRUCTURE_ENTRY_QUERIER_ENTRY_QUERIER_H

#include "common/type_aliases.h"

#include <string>

namespace cath {
	class protein;

	/// \brief TODOCUMENT
	class entry_querier {
	private:
		virtual size_t       do_get_length(const cath::protein &) const = 0;
		virtual double       do_get_gap_penalty_ratio() const = 0;
		virtual size_t       do_num_excluded_on_either_size() const = 0;
		virtual double       do_optimum_single_score() const = 0;
		virtual std::string  do_get_entry_name() const = 0;
		virtual score_type   do_distance_score__offset_1(const cath::protein &,
		                                                 const cath::protein &,
		                                                 const size_t &,
		                                                 const size_t &,
		                                                 const size_t &,
		                                                 const size_t &) const = 0;
		virtual bool         do_are_comparable__offset_1(const cath::protein &,
		                                                 const cath::protein &,
		                                                 const size_t &,
		                                                 const size_t &,
		                                                 const size_t &,
		                                                 const size_t &) const = 0;
		virtual bool         do_are_similar__offset_1(const cath::protein &,
		                                              const cath::protein &,
		                                              const size_t &,
		                                              const size_t &) const = 0;
		virtual bool         do_temp_hacky_is_residue() const = 0;

	public:
		virtual ~entry_querier() noexcept = default;

		size_t       get_length(const cath::protein &) const;
		double       get_gap_penalty_ratio() const;
		size_t       num_excluded_on_either_size() const;
		double       optimum_single_score() const;
		std::string  get_entry_name() const;
		score_type   distance_score__offset_1(const cath::protein &,
		                                      const cath::protein &,
		                                      const size_t &,
		                                      const size_t &,
		                                      const size_t &,
		                                      const size_t &) const;
		bool         are_comparable__offset_1(const cath::protein &,
		                                      const cath::protein &,
		                                      const size_t &,
		                                      const size_t &,
		                                      const size_t &,
		                                      const size_t &) const;
		bool         are_similar__offset_1(const cath::protein &,
		                                   const cath::protein &,
		                                   const size_t &,
		                                   const size_t &) const;
		bool         temp_hacky_is_residue() const;

		/// \brief This appears to be used to multiply up int values to achieve one decimal place
		///        for score calculations
		static constexpr size_t INTEGER_SCALING = 10;
	};

	std::string get_plural_name(const entry_querier &);
	double get_gap_penalty(const entry_querier &);
	bool pair_is_not_excluded(const entry_querier &,
	                          const size_t &,
	                          const size_t &);
	bool pair_is_not_excluded(const size_t &,
	                          const size_t &,
	                          const size_t &);
	size_t num_comparable(const entry_querier &,
	                      const size_t &);
	size_t _num_comparable_impl(const size_t &,
	                            const size_t &);


	/// \brief TODOCUMENT
	inline size_t entry_querier::get_length(const protein &arg_protein ///< TODOCUMENT
	                                        ) const {
		return do_get_length(arg_protein);
	}

	/// \brief Return the ratio of the gap penalty to the optimum single score
	///
	/// If you just want the gap penalty, just use the non-member function get_gap_penalty().
	///
	/// Example: this value is 2, the optimum single score is 50 => gap penalty is 100.
	inline double entry_querier::get_gap_penalty_ratio() const {
		return do_get_gap_penalty_ratio();
	}

	/// \brief TODOCUMENT
	inline double entry_querier::optimum_single_score() const {
		return do_optimum_single_score();
	}

	/// \brief TODOCUMENT
	inline size_t entry_querier::num_excluded_on_either_size() const {
		return do_num_excluded_on_either_size();
	}

	/// \brief TODOCUMENT
	inline std::string entry_querier::get_entry_name() const {
		return do_get_entry_name();
	}

	/// \brief TODOCUMENT
	inline score_type entry_querier::distance_score__offset_1(const protein &arg_protein_a,         ///< TODOCUMENT
	                                                          const protein &arg_protein_b,         ///< TODOCUMENT
	                                                          const size_t  &arg_a_view_from_index, ///< TODOCUMENT
	                                                          const size_t  &arg_b_view_from_index, ///< TODOCUMENT
	                                                          const size_t  &arg_a_dest_to_index,   ///< TODOCUMENT
	                                                          const size_t  &arg_b_dest_to_index    ///< TODOCUMENT
	                                                          ) const {
		return do_distance_score__offset_1(
			arg_protein_a,
			arg_protein_b,
			arg_a_view_from_index,
			arg_b_view_from_index,
			arg_a_dest_to_index,
			arg_b_dest_to_index
		);
	}

	/// \brief TODOCUMENT
	inline bool entry_querier::are_comparable__offset_1(const protein &arg_protein_a,         ///< TODOCUMENT
	                                                    const protein &arg_protein_b,         ///< TODOCUMENT
	                                                    const size_t  &arg_a_view_from_index, ///< TODOCUMENT
	                                                    const size_t  &arg_b_view_from_index, ///< TODOCUMENT
	                                                    const size_t  &arg_a_dest_to_index,   ///< TODOCUMENT
	                                                    const size_t  &arg_b_dest_to_index    ///< TODOCUMENT
	                                                    ) const {
		const bool a_index_pair_not_excluded = pair_is_not_excluded(*this, arg_a_view_from_index, arg_a_dest_to_index);
		const bool b_index_pair_not_excluded = pair_is_not_excluded(*this, arg_b_view_from_index, arg_b_dest_to_index);
		if (!a_index_pair_not_excluded || !b_index_pair_not_excluded) {
			return false;
		}
		return do_are_comparable__offset_1(
			arg_protein_a,
			arg_protein_b,
			arg_a_view_from_index,
			arg_b_view_from_index,
			arg_a_dest_to_index,
			arg_b_dest_to_index
		);
	}

	/// \brief TODOCUMENT
	inline bool entry_querier::are_similar__offset_1(const protein &arg_protein_a,         ///< TODOCUMENT
	                                                 const protein &arg_protein_b,         ///< TODOCUMENT
	                                                 const size_t  &arg_index_a__offset_1, ///< TODOCUMENT
	                                                 const size_t  &arg_index_b__offset_1  ///< TODOCUMENT
	                                                 ) const {
		return do_are_similar__offset_1(arg_protein_a, arg_protein_b, arg_index_a__offset_1, arg_index_b__offset_1);
	}

	/// \brief TODOCUMENT
	///
	/// \relates entry_querier
	///
	/// \todo Remove (the need for) this temporary hacky workaround
	inline bool entry_querier::temp_hacky_is_residue() const {
		return do_temp_hacky_is_residue();
	}

} // namespace cath

#endif
