/// \file
/// \brief The broad_display_colour_spec class definitions

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

#include "broad_display_colour_spec.hpp"

#include <boost/lexical_cast.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/range/adaptor/filtered.hpp>
#include <boost/range/adaptor/map.hpp>
#include <boost/range/join.hpp>

#include "cath/common/algorithm/contains.hpp"
#include "cath/common/algorithm/copy_build.hpp"
#include "cath/common/algorithm/sort_uniq_build.hpp"
#include "cath/common/boost_addenda/range/indices.hpp"
#include "cath/common/exception/invalid_argument_exception.hpp"
#include "cath/display/viewer/viewer.hpp"

using namespace ::cath;
using namespace ::cath::common;
using namespace ::cath::file;

using ::boost::adaptors::filtered;
using ::boost::adaptors::map_keys;
using ::boost::adaptors::map_values;
using ::std::nullopt;
using ::std::ostream;
using ::std::ostringstream;
using ::std::string;

/// \brief Specify the base colour to use
broad_display_colour_spec & broad_display_colour_spec::colour_base(const display_colour &prm_colour,   ///< The colour to use
                                                                   const bool           &prm_overwrite ///< Whether to allow any previous base value to be overwritten without throwing
                                                                   ) {
	if ( ! prm_overwrite && base_clr ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Attempt to re-specify base colour"));
	}
	base_clr = prm_colour;
	return *this;
}

/// \brief Specify the colour to use for a specific structure
broad_display_colour_spec & broad_display_colour_spec::colour_pdb(const size_t         &prm_pdb_index, ///< The index of the structure for which the colour should be set
		                                                          const display_colour &prm_colour,    ///< The colour to use
		                                                          const bool           &prm_overwrite  ///< Whether to allow any previous base value to be overwritten without throwing
		                                                          ) {
	if ( ! prm_overwrite && contains ( clr_of_pdb, prm_pdb_index ) ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Attempt to re-specify colour of structure"));
	}
//	clr_of_pdb.insert( make_pair( prm_pdb_index, prm_colour ) );
	clr_of_pdb[ prm_pdb_index ] = prm_colour;
	return *this;
}

// /// \brief Specify a region and colour to apply to that region in a specific structure
// broad_display_colour_spec & broad_display_colour_spec::colour_region(const size_t         &prm_pdb_index, ///< The index of the structure for which the colour should be set
//                                                                      const region         &prm_region,    ///< The region to colour
//                                                                      const display_colour &prm_colour     ///< The colour to use
//                                                                      ) {
// 	clr_of_region[ prm_region ].emplace_back( prm_region, prm_colour );
// 	return *this;
// }

/// \brief Getter for the base colour
const display_colour_opt & broad_display_colour_spec::get_base_clr() const {
	return base_clr;
}

/// \brief Getter for colour by structure index
const size_display_colour_map & broad_display_colour_spec::get_clr_of_pdb() const {
	return clr_of_pdb;
}

// /// \brief Get the region/colour pairs by structure index
// const size_region_display_colour_pair_vec_map & broad_display_colour_spec::get_clr_of_regions() const {
// 	return clr_of_region;
// }

/// \brief Return whether the specified broad_display_colour_spec has a base colour
///
/// relates broad_display_colour_spec
bool cath::has_base_colour(const broad_display_colour_spec &prm_colour_spec ///< The broad_display_colour_spec to query
                           ) {
	return static_cast<bool>( prm_colour_spec.get_base_clr() );
}

/// \brief Get the base colour of the specified broad_display_colour_spec
///
/// \pre `has_base_colour( prm_colour_spec )`
///
/// relates broad_display_colour_spec
display_colour cath::get_base_colour(const broad_display_colour_spec &prm_colour_spec ///< The broad_display_colour_spec to query
                                     ) {
	return *prm_colour_spec.get_base_clr();
}

/// \brief Get the (optional) colour for the structure of the specified index
///
/// relates broad_display_colour_spec
display_colour_opt cath::get_clr_of_pdb_index(const broad_display_colour_spec &prm_colour_spec, ///< The broad_display_colour_spec to query
                                              const size_t                    &prm_entry        ///< The index of the structure required
                                              ) {
	const size_display_colour_map &colour_of_pdb = prm_colour_spec.get_clr_of_pdb();
	const bool                     has_entry     = contains( colour_of_pdb, prm_entry );
	return has_entry ? display_colour_opt( colour_of_pdb.at( prm_entry ) )
	                 : display_colour_opt( nullopt                       );
}

/// \brief Get a list of all colours used to colour whole structures
///
/// relates broad_display_colour_spec
display_colour_vec cath::get_pdb_colours(const broad_display_colour_spec &prm_colour_spec ///< The broad_display_colour_spec to query
                                         ) {
	return sort_uniq_build<display_colour_vec>(
		prm_colour_spec.get_clr_of_pdb() | map_values
	);
}

/// \brief Get a list of the indices of all PDBs of the specified colour
///
/// relates broad_display_colour_spec
size_vec cath::get_pdbs_of_colour(const broad_display_colour_spec &prm_colour_spec, ///< The broad_display_colour_spec to query
                                  const display_colour            &prm_colour       ///< The colour to query
                                  ) {
	return copy_build<size_vec>(
		prm_colour_spec.get_clr_of_pdb()
		| filtered( [&] (const size_display_colour_map_val &x) { return ( x.second == prm_colour ); } )
		| map_keys
	);
}

/// Write to specified stream the text to define all the specified colours with the specified viewer
///
/// \param prm_colours         The list of colours
/// \param prm_viewer          The viewer to specify how to render the instructions
/// \param prm_os              The ostream to which the instructions should be written
/// \param prm_colour_category The category of colouring (structure-only or structure-or-residue)
void cath::define_all_colours( const display_colour_vec &prm_colours,
                               const viewer &            prm_viewer,
                               ostream &                 prm_os,
                               const colour_category &   prm_colour_category ) {
	const size_t num_colours = prm_colours.size();
	for ( const size_t &colour_ctr : indices( num_colours ) ) {
		prm_viewer.define_colour(
		  prm_os, prm_colours[ colour_ctr ], generate_colour_name( colour_ctr, num_colours, prm_colour_category ) );
	}
}

/// \brief Write broad-level (ie not residue-specific) colouring instructions for the specified viewer to the specified ostream
///        using the specified broad_display_colour_spec and (viewer-cleaned) names
///        in the context of the specified list of colours
///
/// relates broad_display_colour_spec
void cath::detail::colour_base_and_pdbs_impl(const display_colour_vec        &prm_colours,                  ///< The list of colours
                                             const broad_display_colour_spec &prm_colour_spec,              ///< The broad_display_colour_spec to use
                                             const viewer                    &prm_viewer,                   ///< The viewer to specify how to render the instructions
                                             const pdb_list                  &/*prm_pdbs*/,                 ///< The key regions of the structures
                                             const str_vec                   &prm_cleaned_names_for_viewer, ///< The (viewer-cleaned) names for the structures to be coloured
                                             const colour_category           &prm_colour_category,          ///< The category of colouring (structure-only or structure-or-residue)
                                             ostream                         &prm_os                        ///< The ostream to which the instructions should be written
                                             ) {
	define_all_colours( prm_colours, prm_viewer, prm_os, prm_colour_category );

	const size_t num_colours = prm_colours.size();
	if ( has_base_colour( prm_colour_spec ) ) {
		prm_viewer.define_colour( prm_os, get_base_colour( prm_colour_spec ), base_colour_name() );
		prm_os << prm_viewer.get_colour_base_str( base_colour_name() );
	}

	for (const size_t colour_ctr : indices( num_colours ) ) {
		const size_vec pdb_indices = get_pdbs_of_colour( prm_colour_spec, prm_colours[ colour_ctr ] );
		for (const size_t &pdb_index : pdb_indices) {
			prm_os << prm_viewer.get_colour_pdb_str(
				generate_colour_name( colour_ctr, num_colours, prm_colour_category ),
				prm_cleaned_names_for_viewer[ pdb_index  ]
			);
		}
	}
}

/// \brief Write broad-level (ie not residue-specific) colouring instructions for the specified viewer to the specified ostream
///        using the specified broad_display_colour_spec and (viewer-cleaned) names
///
/// relates broad_display_colour_spec
void cath::colour_viewer_with_spec(const broad_display_colour_spec &prm_broad_spec,               ///< The specification of colours to use
                                   const viewer                    &prm_viewer,                   ///< The viewer to specify how to render the instructions
                                   const pdb_list                  &prm_pdbs,                     ///< The key regions of the structures
                                   const str_vec                   &prm_cleaned_names_for_viewer, ///< The (viewer-cleaned) names for the structures to be coloured
                                   ostream                         &prm_os                        ///< The ostream to which the instructions should be written
                                   ) {
	detail::colour_base_and_pdbs_impl(
		get_pdb_colours( prm_broad_spec ),
		prm_broad_spec,
		prm_viewer,
		prm_pdbs,
		prm_cleaned_names_for_viewer,
		colour_category::STRUC_ONLY,
		prm_os
	);
}
