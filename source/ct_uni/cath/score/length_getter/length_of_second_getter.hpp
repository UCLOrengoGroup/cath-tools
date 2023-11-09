/// \file
/// \brief The length_of_second_getter class header

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

#ifndef CATH_TOOLS_SOURCE_CT_UNI_CATH_SCORE_LENGTH_GETTER_LENGTH_OF_SECOND_GETTER_HPP
#define CATH_TOOLS_SOURCE_CT_UNI_CATH_SCORE_LENGTH_GETTER_LENGTH_OF_SECOND_GETTER_HPP

#include "cath/score/length_getter/sym_protein_only_length_getter.hpp"

namespace cath::score {

	/// \brief TODOCUMENT
	class length_of_second_getter final : public protein_only_length_getter {
	  private:
		[[nodiscard]] std::unique_ptr<protein_only_length_getter> do_protein_only_clone() const final;

		[[nodiscard]] size_t do_get_length( const protein &, const protein & ) const final;

		[[nodiscard]] length_getter_category do_get_length_getter_category() const final;

		[[nodiscard]] std::string do_get_choice_adjective() const final;

		[[nodiscard]] bool do_less_than_with_same_dynamic_type( const length_getter & ) const final;
	};

} // namespace cath::score

#endif // CATH_TOOLS_SOURCE_CT_UNI_CATH_SCORE_LENGTH_GETTER_LENGTH_OF_SECOND_GETTER_HPP
