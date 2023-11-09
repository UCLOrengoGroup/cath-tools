/// \file
/// \brief The pseudo_string_score class header

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

#ifndef CATH_TOOLS_SOURCE_CT_UNI_CATH_SCORE_ALIGNED_PAIR_SCORE_PSEUDO_STRING_SCORE_HPP
#define CATH_TOOLS_SOURCE_CT_UNI_CATH_SCORE_ALIGNED_PAIR_SCORE_PSEUDO_STRING_SCORE_HPP

#include "cath/score/aligned_pair_score/aligned_pair_score.hpp"
#include "cath/score/aligned_pair_score/detail/score_common_coord_handler.hpp"

namespace cath::score {

	/// \brief .
	///
	/// \ingroup cath_score_aligned_pair_score_group
	class pseudo_string_score : public aligned_pair_score {
	private:
		friend class boost::serialization::access;

		template<class archive> void serialize(archive &ar,
		                                       const size_t /*version*/
		                                       ) {
			ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP( aligned_pair_score );
			ar & BOOST_SERIALIZATION_NVP( score_name             );
			ar & BOOST_SERIALIZATION_NVP( higher_is_better_value );
		}

		/// \brief TODOCUMENT
		std::string score_name;

		/// \brief TODOCUMENT
		boost::logic::tribool higher_is_better_value;

		[[nodiscard]] std::unique_ptr<aligned_pair_score> do_clone() const final;

		[[nodiscard]] boost::logic::tribool do_higher_is_better() const final;
		[[nodiscard]] score_value do_calculate( const align::alignment &, const protein &, const protein & ) const final;
		[[nodiscard]] std::string       do_description() const final;
		[[nodiscard]] std::string       do_id_name() const final;
		[[nodiscard]] str_bool_pair_vec do_short_name_suffixes() const final;
		[[nodiscard]] std::string       do_long_name() const final;
		[[nodiscard]] std::string       do_reference() const final;

		// std::unique_ptr<aligned_pair_score> do_build_from_short_name_spec(const std::string &) const final;

		[[nodiscard]] bool do_less_than_with_same_dynamic_type( const aligned_pair_score & ) const final;

	  public:
		explicit pseudo_string_score(std::string,
		                             boost::logic::tribool);
	};

	bool operator<(const pseudo_string_score &,
	               const pseudo_string_score &);

} // namespace cath::score

#endif // CATH_TOOLS_SOURCE_CT_UNI_CATH_SCORE_ALIGNED_PAIR_SCORE_PSEUDO_STRING_SCORE_HPP
