/// \file
/// \brief The overlap_score class header

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

#ifndef CATH_TOOLS_SOURCE_CT_UNI_CATH_SCORE_ALIGNED_PAIR_SCORE_OVERLAP_SCORE_HPP
#define CATH_TOOLS_SOURCE_CT_UNI_CATH_SCORE_ALIGNED_PAIR_SCORE_OVERLAP_SCORE_HPP

#include <boost/numeric/conversion/cast.hpp>
#include <boost/range/join.hpp>

#include "cath/common/algorithm/constexpr_find.hpp"
#include "cath/common/algorithm/copy_build.hpp"
#include "cath/common/clone/clone_ptr.hpp"
#include "cath/common/clone/make_uptr_clone.hpp"
#include "cath/score/aligned_pair_score/aligned_pair_score.hpp"
#include "cath/score/aligned_pair_score/overlap_type.hpp"
#include "cath/score/length_getter/num_aligned_length_getter.hpp"

namespace cath::score {

	/// \brief Calculate (and represent) the number of aligned residues in an alignment.
	///
	/// \todo Update the names and description based on the details of the score_common_coord_handler.
	template <overlap_type OL>
	class overlap_score : public aligned_pair_score {
	private:
		friend class boost::serialization::access;

		[[nodiscard]] std::unique_ptr<aligned_pair_score> do_clone() const final;

		[[nodiscard]] boost::logic::tribool do_higher_is_better() const final;
		[[nodiscard]] score_value do_calculate( const align::alignment &, const protein &, const protein & ) const final;
		[[nodiscard]] std::string       do_description() const final;
		[[nodiscard]] std::string       do_id_name() const final;
		[[nodiscard]] str_bool_pair_vec do_short_name_suffixes() const final;
		[[nodiscard]] std::string       do_long_name() const final;

		// std::unique_ptr<aligned_pair_score> do_build_from_short_name_spec(const std::string &) const final;

		[[nodiscard]] bool do_less_than_with_same_dynamic_type( const aligned_pair_score & ) const final;

	  public:
		/// \brief TODOCUMENT
		using numerator_type   = detail::length_getter_of_enum_t< std::get<1>( common::constexpr_find( all_overlap_types, OL ) ) >;

		/// \brief TODOCUMENT
		using denominator_type = detail::length_getter_of_enum_t< std::get<2>( common::constexpr_find( all_overlap_types, OL ) ) >;
	};


	/// \brief A standard do_clone method.
	template <overlap_type OL>
	std::unique_ptr<aligned_pair_score> overlap_score<OL>::do_clone() const {
		return { common::make_uptr_clone( *this ) };
	}

	/// \brief Concrete implementation that records that more aligned residues is generally better
	template <overlap_type OL>
	boost::tribool overlap_score<OL>::do_higher_is_better() const {
		return true;
	}

	/// \brief Concrete implementation for calculating the number of common residues defined by this alignment
	///
	/// This uses the score_common_coord_handler (and hence its policies)
	template <overlap_type OL>
	score_value overlap_score<OL>::do_calculate(const align::alignment &prm_alignment, ///< The pair alignment to be scored
	                                            const protein          &prm_protein_a, ///< The protein associated with the first  half of the alignment
	                                            const protein          &prm_protein_b  ///< The protein associated with the second half of the alignment
	                                            ) const {
		const size_t numerator_value   = numerator_type  ().get_length( prm_alignment, prm_protein_a, prm_protein_b );
		const size_t denominator_value = denominator_type().get_length( prm_alignment, prm_protein_a, prm_protein_b );
		return 100.0 * boost::numeric_cast<score_value>( numerator_value ) / boost::numeric_cast<score_value>( denominator_value );
	}

	/// \brief Concrete implementation that describes what this score means
	template <overlap_type OL>
	std::string overlap_score<OL>::do_description() const {
		return { "Percentage overlap calculated as 100.0 * "
			+ numerator_type  ().full_short_name()
			+ " divided by "
			+ denominator_type().full_short_name()
		};
	}

	/// \brief TODOCUMENT
	template <overlap_type OL>
	std::string overlap_score<OL>::do_id_name() const {
		return "overlap";
	}

	/// \brief TODOCUMENT
	template <overlap_type OL>
	str_bool_pair_vec overlap_score<OL>::do_short_name_suffixes() const {
		const str_bool_pair_vec numerator_id_suffix  { { numerator_type  ().id_name(), true } };
		const str_bool_pair_vec denominator_id_suffix{ { denominator_type().id_name(), true } };
		return common::copy_build<str_bool_pair_vec>(
			boost::range::join(
				boost::range::join(
					numerator_id_suffix,
					numerator_type  ().short_name_suffixes()
				),
				boost::range::join(
					denominator_id_suffix,
					denominator_type().short_name_suffixes()
				)
			)
		);
	}

	/// \brief Concrete implementation providing long name
	template <overlap_type OL>
	std::string overlap_score<OL>::do_long_name() const {
		return { "Overlap of "
			+ numerator_type  ().human_friendly_short_name()
			+ " over "
			+ denominator_type().human_friendly_short_name()
		};
	}

	/// \brief TODOCUMENT
	template <overlap_type OL>
	bool overlap_score<OL>::do_less_than_with_same_dynamic_type(const aligned_pair_score &prm_aligned_pair_score ///< TODOCUMENT
	                                                        ) const {
		const auto &casted_aligned_pair_score = dynamic_cast< decltype( *this ) >( prm_aligned_pair_score );
		return ( *this < casted_aligned_pair_score );
	}

	/// \brief TODOCUMENT
	///
	/// \relates overlap_score
	template <overlap_type OL>
	bool operator<(const overlap_score<OL> &/*prm_overlap_score_a*/, ///< TODOCUMENT
	               const overlap_score<OL> &/*prm_overlap_score_b*/  ///< TODOCUMENT
	               ) {
		return false;
	}

	/// \brief TODOCUMENT
	using naive_overlap  = overlap_score< overlap_type::SHORTER_OVER_LONGER >;

	/// \brief TODOCUMENT
	using local_overlap  = overlap_score< overlap_type::NUM_ALIGNED_OVER_SHORTER >;

	/// \brief TODOCUMENT
	using global_overlap = overlap_score< overlap_type::NUM_ALIGNED_OVER_LONGER >;

} // namespace cath::score

#endif // CATH_TOOLS_SOURCE_CT_UNI_CATH_SCORE_ALIGNED_PAIR_SCORE_OVERLAP_SCORE_HPP
