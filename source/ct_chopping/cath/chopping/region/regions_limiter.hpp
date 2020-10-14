/// \file
/// \brief The regions_limiter class header

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

#ifndef _CATH_TOOLS_SOURCE_CHOPPING_REGION_REGIONS_LIMITER_HPP
#define _CATH_TOOLS_SOURCE_CHOPPING_REGION_REGIONS_LIMITER_HPP

#include <boost/optional/optional.hpp>

#include "cath/chopping/chopping_type_aliases.hpp"
#include "cath/common/type_aliases.hpp"

namespace cath { class residue_id; }

namespace cath {
	namespace chop {

		namespace detail {

			/// \brief Represent whether a region has been seen during the use of a regions_limiter
			enum class region_seen : bool {
				YES, ///< Yes, the region has     been seen
				NO   ///< No,  the region has not been seen
			};

			/// \brief Type alias for a vector of region_seen values
			using region_seen_vec = std::vector<region_seen>;

			/// \brief Type alias for an optional region_seen_vec
			using region_seen_vec_opt = boost::optional<region_seen_vec>;

		} // namespace chop

		/// \brief Implement limiting a series of residue_ids to an (optional) list of regions
		///
		/// The regions are stored by reference.
		///
		/// See regions_limiter_test_suite for examples of use.
		///
		/// \invariant If regions are specified, they must all satisfy: `[] (const region &x) { return has_chain_label( x ); }`
		///            else an invalid_argument_exception will be thrown on exception
		class regions_limiter final {
		private:
			/// \brief An optional reference_wrapper to a vector of regions
			region_vec_cref_opt regions;

			/// \brief Store which of the specified regions have been seen (or none if no regions were specified)
			detail::region_seen_vec_opt regions_seen;

			/// \brief The currently active region, or none if none is active
			size_opt active_region_idx;

			void sanity_check() const;

		public:
			regions_limiter() = default;
			explicit regions_limiter(const region_vec &);
			explicit regions_limiter(const region_vec_opt &);
			explicit regions_limiter(const region_vec_opt &&) = delete;

			/// \brief Prevent construction from an rvalue because the region_vec is stored by reference
			regions_limiter(const region_vec &&) = delete;

			bool update_residue_is_included(const residue_id &);

			region_vec unseen_regions() const;
		};

		str_opt warn_str_if_specified_regions_remain_unseen(const regions_limiter &);

		void warn_if_specified_regions_remain_unseen(const regions_limiter &);

	} // namespace chop
} // namespace cath

#endif
