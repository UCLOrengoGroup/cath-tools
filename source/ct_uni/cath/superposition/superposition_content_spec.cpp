/// \file
/// \brief The superposition_content_spec class definitions

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

#include "superposition_content_spec.hpp"

#include "cath/superposition/options/superposition_content_options_block.hpp"

#include <utility>

using namespace ::cath;
using namespace ::cath::sup;

using ::boost::none;

constexpr supn_regions_context superposition_content_spec::DEFAULT_REGIONS_CONTEXT;
constexpr double               superposition_content_spec::DEFAULT_DNA_MAX_DIST;
constexpr double               superposition_content_spec::DEFAULT_ORGANIC_MAX_DIST;

/// \brief Ctor
superposition_content_spec::superposition_content_spec(const supn_regions_context &prm_regions_context,                ///< The context to include when showing the specified region(s) of structure
                                                       doub_opt                    prm_include_dna_within_distance,    ///< The maximum distance from the specified region(s) to DNA/RNA for the that to be included (or none to always exclude)
                                                       doub_opt                    prm_include_organic_within_distance ///< The maximum distance from the specified region(s) to ligands for the that to be included (or none to always exclude)
                                                       ) : regions_context                { prm_regions_context                              },
                                                           include_dna_within_distance    { std::move( prm_include_dna_within_distance     ) },
                                                           include_organic_within_distance{ std::move( prm_include_organic_within_distance ) } {
}

/// \brief Getter for the context to include when showing the specified region(s) of structure
const supn_regions_context & superposition_content_spec::get_regions_context() const {
	return regions_context;
}

/// \brief Getter for the maximum distance from the specified region(s) to DNA/RNA for the that to be included (or none to always exclude)
const doub_opt & superposition_content_spec::get_include_dna_within_distance() const {
	return include_dna_within_distance;
}

/// \brief Getter for the maximum distance from the specified region(s) to ligands for the that to be included (or none to always exclude)
const doub_opt & superposition_content_spec::get_include_organic_within_distance() const {
	return include_organic_within_distance;
}

/// \brief Setter for the context to include when showing the specified region(s) of structure
superposition_content_spec & superposition_content_spec::set_regions_context(const supn_regions_context &prm_regions_context ///< The context to include when showing the specified region(s) of structure
                                                                             ) {
	regions_context = prm_regions_context;
	return *this;
}

/// \brief Setter for the maximum distance from the specified region(s) to DNA/RNA for the that to be included (or none to always exclude)
superposition_content_spec & superposition_content_spec::set_include_dna_within_distance(const doub_opt &prm_include_dna_within_distance ///< The maximum distance from the specified region(s) to DNA/RNA for the that to be included (or none to always exclude)
                                                                                         ) {
	include_dna_within_distance = prm_include_dna_within_distance;
	return *this;
}

/// \brief Setter for the maximum distance from the specified region(s) to ligands for the that to be included (or none to always exclude)
superposition_content_spec & superposition_content_spec::set_include_organic_within_distance(const doub_opt &prm_include_organic_within_distance ///< The maximum distance from the specified region(s) to ligands for the that to be included (or none to always exclude)
                                                                                             ) {
	include_organic_within_distance = prm_include_organic_within_distance;
	return *this;
}

/// \brief Return a string explaining why the specified superposition_content_spec is invalid or none if it isn't
///
/// \relates superposition_content_spec
str_opt cath::sup::get_invalid_description(const superposition_content_spec &prm_spec ///< The superposition_content_spec to query
                                           ) {
	if ( prm_spec.get_include_dna_within_distance() && prm_spec.get_include_dna_within_distance() < 0.0 ) {
		return
			"The "
			+ superposition_content_options_block::PO_INCLUDE_DNA_WITHIN_DISTANCE
			+ " distance cannot be negative";
	}
	if ( prm_spec.get_include_organic_within_distance() && prm_spec.get_include_organic_within_distance() < 0.0 ) {
		return
			"The "
			+ superposition_content_options_block::PO_INCLUDE_ORGANIC_WITHIN_DISTANCE
			+ " distance cannot be negative";
	}
	return none;
}