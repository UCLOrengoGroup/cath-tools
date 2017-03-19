/// \file
/// \brief The superposition_content_spec class header

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

#ifndef _CATH_TOOLS_SOURCE_SUPERPOSITION_SUPERPOSITION_CONTENT_SPEC_H
#define _CATH_TOOLS_SOURCE_SUPERPOSITION_SUPERPOSITION_CONTENT_SPEC_H

#include <boost/optional.hpp>

#include "superposition/supn_regions_context.hpp"

namespace cath {

	/// \todo Move this to common/type_aliases.hpp
	using doub_opt = boost::optional<double>;

	namespace sup {

		/// \brief Specify what should be included in superpositions other than the specified region(s) of structure
		class superposition_content_spec final {
		private:
			/// \brief The context to include when showing the specified region(s) of structure
			supn_regions_context regions_context     = DEFAULT_REGIONS_CONTEXT;

			/// \brief The maximum distance from the specified region(s) to DNA/RNA for the that to be included (or none to always exclude)
			doub_opt include_dna_within_distance     = DEFAULT_DNA_MAX_DIST;

			/// \brief The maximum distance from the specified region(s) to ligands for the that to be included (or none to always exclude)
			doub_opt include_organic_within_distance = DEFAULT_ORGANIC_MAX_DIST;

		public:
			/// \brief The default value for the context to include when showing the specified region(s) of structure
			static constexpr supn_regions_context DEFAULT_REGIONS_CONTEXT  = supn_regions_context::ALONE;

			/// \brief The default value for the maximum distance from the specified region(s) to DNA/RNA for the that to be included (or none to always exclude)
			static constexpr double               DEFAULT_DNA_MAX_DIST     =  4.0;

			/// \brief The default value for the maximum distance from the specified region(s) to ligands for the that to be included (or none to always exclude)
			static constexpr double               DEFAULT_ORGANIC_MAX_DIST = 10.0;

			superposition_content_spec() = default;

			explicit superposition_content_spec(const supn_regions_context &,
			                                    const doub_opt & = DEFAULT_DNA_MAX_DIST,
			                                    const doub_opt & = DEFAULT_ORGANIC_MAX_DIST);

			const supn_regions_context & get_regions_context() const;
			const doub_opt & get_include_dna_within_distance() const;
			const doub_opt & get_include_organic_within_distance() const;

			superposition_content_spec & set_regions_context(const supn_regions_context &);
			superposition_content_spec & set_include_dna_within_distance(const doub_opt &);
			superposition_content_spec & set_include_organic_within_distance(const doub_opt &);
		};

	} // namespace sup
} // namespace cath

#endif
