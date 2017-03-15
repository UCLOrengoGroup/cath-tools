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

using namespace cath;
using namespace cath::sup;

constexpr supn_regions_context superposition_content_spec::DEFAULT_REGIONS_CONTEXT;
constexpr double               superposition_content_spec::DEFAULT_DNA_MAX_DIST;
constexpr double               superposition_content_spec::DEFAULT_ORGANIC_MAX_DIST;

/// \brief TODOCUMENT
superposition_content_spec::superposition_content_spec(const supn_regions_context &arg_regions_context,                ///< The context to include when showing the specified region(s) of structure
                                                       const doub_opt             &arg_include_dna_within_distance,    ///< The maximum distance from the specified region(s) to DNA/RNA for the that to be included (or none to always exclude)
                                                       const doub_opt             &arg_include_organic_within_distance ///< The maximum distance from the specified region(s) to ligands for the that to be included (or none to always exclude)
                                                       ) : regions_context                { arg_regions_context                 },
                                                           include_dna_within_distance    { arg_include_dna_within_distance     },
                                                           include_organic_within_distance{ arg_include_organic_within_distance } {
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
superposition_content_spec & superposition_content_spec::set_regions_context(const supn_regions_context &arg_regions_context ///< The context to include when showing the specified region(s) of structure
                                                                             ) {
	regions_context = arg_regions_context;
	return *this;
}

/// \brief Setter for the maximum distance from the specified region(s) to DNA/RNA for the that to be included (or none to always exclude)
superposition_content_spec & superposition_content_spec::set_include_dna_within_distance(const doub_opt &arg_include_dna_within_distance ///< The maximum distance from the specified region(s) to DNA/RNA for the that to be included (or none to always exclude)
                                                                                         ) {
	include_dna_within_distance = arg_include_dna_within_distance;
	return *this;
}

/// \brief Setter for the maximum distance from the specified region(s) to ligands for the that to be included (or none to always exclude)
superposition_content_spec & superposition_content_spec::set_include_organic_within_distance(const doub_opt &arg_include_organic_within_distance ///< The maximum distance from the specified region(s) to ligands for the that to be included (or none to always exclude)
                                                                                             ) {
	include_organic_within_distance = arg_include_organic_within_distance;
	return *this;
}
