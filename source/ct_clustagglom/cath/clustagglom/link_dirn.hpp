/// \file
/// \brief The link_dirn header

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

#ifndef _CATH_TOOLS_SOURCE_CT_CLUSTAGGLOM_CATH_CLUSTAGGLOM_LINK_DIRN_HPP
#define _CATH_TOOLS_SOURCE_CT_CLUSTAGGLOM_CATH_CLUSTAGGLOM_LINK_DIRN_HPP

#include <boost/any.hpp>

#include "cath/common/algorithm/constexpr_is_uniq.hpp"
#include "cath/common/cpp20/make_array.hpp"
#include "cath/common/detail/maybe_unused_namespace_scope_constexpr.hpp"
#include "cath/common/type_aliases.hpp"

namespace cath {
	namespace clust {

		/// \brief Whether the link between items/clusters is a measure of dissimilarity or strength (ie similarity)
		enum class link_dirn : bool {
			DISSIMILARITY, ///< The link is a measure of dissimilarity
			STRENGTH,      ///< The link is a measure of strength (ie similarity)
		};

		/// \brief A constexpr list of all link_dirns
		static constexpr auto all_link_dirns = common::make_array(
			link_dirn::DISSIMILARITY,
			link_dirn::STRENGTH
		);

		// Compile-time check that there aren't any duplicates in all_link_dirns
		static_assert( common::constexpr_is_uniq( all_link_dirns ), "all_link_dirns shouldn't contain repeated values" );

		/// \brief Store a constexpr record of the number of link_dirns
		static constexpr size_t num_link_dirns = std::tuple_size_v< decltype( all_link_dirns ) >;
		MAYBE_UNUSED_NAMESPACE_SCOPE_CONSTEXPR( num_link_dirns )

		namespace detail {

			/// \brief Class with static getter for a map from name to link_dirn
			struct link_dirn_by_name final {
				static std::map<std::string, link_dirn> get();
			};

		} // namespace detail

		std::string to_string(const link_dirn &);

		std::istream & operator>>(std::istream &,
		                          link_dirn &);

		std::string description_of_link_dirn(const link_dirn &);

		void validate(boost::any &,
		              const str_vec &,
		              link_dirn *,
		              int);

	} // namespace clust
} // namespace cath


#endif // _CATH_TOOLS_SOURCE_CT_CLUSTAGGLOM_CATH_CLUSTAGGLOM_LINK_DIRN_HPP
