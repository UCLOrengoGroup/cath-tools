/// \file
/// \brief The make_clone header

/// \copyright
/// Tony Lewis's Common C++ Library Code (here imported into the CATH Tools project and then tweaked, eg namespaced in cath)
/// Copyright (C) 2007, Tony Lewis
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

#ifndef CATH_TOOLS_SOURCE_CT_UNI_CATH_SCORE_LENGTH_GETTER_LENGTH_GETTER_MAKE_CLONE_HPP
#define CATH_TOOLS_SOURCE_CT_UNI_CATH_SCORE_LENGTH_GETTER_LENGTH_GETTER_MAKE_CLONE_HPP

#include "cath/common/clone/detail/make_clone.hpp"

// clang-format off
namespace cath::score { class protein_only_length_getter; }
namespace cath::score { class sym_protein_only_length_getter; }
// clang-format on

namespace cath::common::detail {

	template <>
	std::unique_ptr<score::protein_only_length_getter> make_clone( const score::protein_only_length_getter & );

	template <>
	std::unique_ptr<score::sym_protein_only_length_getter> make_clone( const score::sym_protein_only_length_getter & );

} // namespace cath::common::detail

#endif // CATH_TOOLS_SOURCE_CT_UNI_CATH_SCORE_LENGTH_GETTER_LENGTH_GETTER_MAKE_CLONE_HPP
