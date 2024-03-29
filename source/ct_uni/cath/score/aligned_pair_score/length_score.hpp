/// \file
/// \brief The length_score class header

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

#ifndef CATH_TOOLS_SOURCE_CT_UNI_CATH_SCORE_ALIGNED_PAIR_SCORE_LENGTH_SCORE_HPP
#define CATH_TOOLS_SOURCE_CT_UNI_CATH_SCORE_ALIGNED_PAIR_SCORE_LENGTH_SCORE_HPP

#include "cath/common/clone/clone_ptr.hpp"
#include "cath/score/aligned_pair_score/aligned_pair_score.hpp"
#include "cath/score/length_getter/length_getter_make_clone.hpp"
#include "cath/score/length_getter/num_aligned_length_getter.hpp"

namespace cath::score {

	/// \brief Calculate (and represent) the number of aligned residues in an alignment.
	///
	/// \todo Update the names and description based on the details of the score_common_coord_handler.
	class length_score : public aligned_pair_score {
	private:
		friend class boost::serialization::access;

		template<class archive> void serialize(archive &ar,
		                                       const size_t /*version*/
		                                       ) {
			ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP( aligned_pair_score );
			ar & BOOST_SERIALIZATION_NVP( length_getter_ptr );
		}

		/// \brief TODOCUMENT
		common::clone_ptr<length_getter> length_getter_ptr;

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
		explicit length_score(const length_getter &);

		[[nodiscard]] const length_getter &get_length_getter() const;

		[[nodiscard]] std::string description_brackets_string() const;
	};

	bool operator<(const length_score &,
	               const length_score &);

} // namespace cath::score

#endif // CATH_TOOLS_SOURCE_CT_UNI_CATH_SCORE_ALIGNED_PAIR_SCORE_LENGTH_SCORE_HPP
