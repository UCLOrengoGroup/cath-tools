/// \file
/// \brief The supn_regions_context class header

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

#ifndef CATH_TOOLS_SOURCE_CT_UNI_CATH_SUPERPOSITION_SUPN_REGIONS_CONTEXT_HPP
#define CATH_TOOLS_SOURCE_CT_UNI_CATH_SUPERPOSITION_SUPN_REGIONS_CONTEXT_HPP

#include <boost/any.hpp>

#include "cath/chopping/chopping_type_aliases.hpp"
#include "cath/common/algorithm/constexpr_is_uniq.hpp"
#include "cath/common/type_aliases.hpp"

// clang-format off
namespace cath::sup { class superposition_content_spec; }
// clang-format on

namespace cath::sup {

	/// \brief The context to include in superpositions when showing the specified region(s) of structure
	enum class supn_regions_context : char {
		ALONE,    ///< Only include the specified region(s) of structure
		IN_CHAIN, ///< Include the full chain of the specified region(s) of structure
		IN_PDB    ///< Include the full PDB of the specified region(s) of structure
	};

	/// \brief A constexpr list of all supn_regions_contexts
	static constexpr std::array<supn_regions_context, 3> all_supn_regions_contexts { {
		supn_regions_context::ALONE,
		supn_regions_context::IN_CHAIN,
		supn_regions_context::IN_PDB,
	} };

	// Compile-time check that there aren't any duplicates in all_supn_regions_contexts
	static_assert( common::constexpr_is_uniq( all_supn_regions_contexts ), "all_supn_regions_contexts shouldn't contain repeated values" );

	/// \brief Store a constexpr record of the number of supn_regions_contexts
	inline constexpr size_t num_supn_regions_contexts = std::tuple_size_v< decltype( all_supn_regions_contexts ) >;

	namespace detail {

		/// \brief Class with static getter for a map from name to hits_input_format_tag
		struct supn_regions_context_by_name final {
			static std::map<std::string, supn_regions_context> get();
		};

	} // namespace detail

	/// \brief Class with static getter for a list of all the supn_regions_context names
	struct all_supn_regions_context_names final {
		static str_vec get();
	};

	chop::region_vec_opt get_regions_expanded_for_context(const chop::region_vec &,
	                                                      const supn_regions_context &);

	chop::region_vec_opt get_regions_expanded_for_context(const chop::region_vec &,
	                                                      const superposition_content_spec &);

	std::string to_string(const supn_regions_context &);

	std::ostream & operator<<(std::ostream &,
	                          const supn_regions_context &);

	std::istream & operator>>(std::istream &,
	                          supn_regions_context &);

	std::string description_of_supn_regions_context(const supn_regions_context &);

	void validate(boost::any &,
	              const str_vec &,
	              supn_regions_context *,
	              int);

} // namespace cath::sup

#endif // CATH_TOOLS_SOURCE_CT_UNI_CATH_SUPERPOSITION_SUPN_REGIONS_CONTEXT_HPP
