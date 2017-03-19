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
#include <boost/range/irange.hpp>
#include <boost/range/join.hpp>

#include "common/algorithm/contains.hpp"
#include "common/algorithm/copy_build.hpp"
#include "common/algorithm/sort_uniq_build.hpp"
#include "common/size_t_literal.hpp"
#include "display/viewer/viewer.hpp"
#include "exception/invalid_argument_exception.hpp"

using namespace cath;
using namespace cath::common;
using namespace cath::file;

using boost::adaptors::filtered;
using boost::adaptors::map_keys;
using boost::adaptors::map_values;
using boost::irange;
using boost::lexical_cast;
using boost::none;
using boost::numeric_cast;
using boost::range::join;
using std::max;
using std::ostream;
using std::ostringstream;
using std::setfill;
using std::setw;
using std::string;

/// \brief Specify the base colour to use
broad_display_colour_spec & broad_display_colour_spec::colour_base(const display_colour &arg_colour,   ///< The colour to use
                                                                   const bool           &arg_overwrite ///< Whether to allow any previous base value to be overwritten without throwing
                                                                   ) {
	if ( ! arg_overwrite && base_clr ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Attempt to re-specify base colour"));
	}
	base_clr = arg_colour;
	return *this;
}

/// \brief Specify the colour to use for a specific structure
broad_display_colour_spec & broad_display_colour_spec::colour_pdb(const size_t         &arg_pdb_index, ///< The index of the structure for which the colour should be set
		                                                          const display_colour &arg_colour,    ///< The colour to use
		                                                          const bool           &arg_overwrite  ///< Whether to allow any previous base value to be overwritten without throwing
		                                                          ) {
	if ( ! arg_overwrite && contains ( clr_of_pdb, arg_pdb_index ) ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Attempt to re-specify colour of structure"));
	}
//	clr_of_pdb.insert( make_pair( arg_pdb_index, arg_colour ) );
	clr_of_pdb[ arg_pdb_index ] = arg_colour;
	return *this;
}

// /// \brief Specify a region and colour to apply to that region in a specific structure
// broad_display_colour_spec & broad_display_colour_spec::colour_region(const size_t         &arg_pdb_index, ///< The index of the structure for which the colour should be set
//                                                                      const region         &arg_region,    ///< The region to colour
//                                                                      const display_colour &arg_colour     ///< The colour to use
//                                                                      ) {
// 	clr_of_region[ arg_region ].emplace_back( arg_region, arg_colour );
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
bool cath::has_base_colour(const broad_display_colour_spec &arg_colour_spec ///< The broad_display_colour_spec to query
                           ) {
	return static_cast<bool>( arg_colour_spec.get_base_clr() );
}

/// \brief Get the base colour of the specified broad_display_colour_spec
///
/// \pre `has_base_colour( arg_colour_spec )`
///
/// relates broad_display_colour_spec
display_colour cath::get_base_colour(const broad_display_colour_spec &arg_colour_spec ///< The broad_display_colour_spec to query
                                     ) {
	return *arg_colour_spec.get_base_clr();
}

/// \brief Get the (optional) colour for the structure of the specified index
///
/// relates broad_display_colour_spec
display_colour_opt cath::get_clr_of_pdb_index(const broad_display_colour_spec &arg_colour_spec, ///< The broad_display_colour_spec to query
                                              const size_t                    &arg_entry        ///< The index of the structure required
                                              ) {
	const size_display_colour_map &colour_of_pdb = arg_colour_spec.get_clr_of_pdb();
	const bool                     has_entry     = contains( colour_of_pdb, arg_entry );
	return has_entry ? display_colour_opt( colour_of_pdb.at( arg_entry ) )
	                 : display_colour_opt( none                          );
}

/// \brief Get a list of all colours used to colour whole structures
///
/// relates broad_display_colour_spec
display_colour_vec cath::get_pdb_colours(const broad_display_colour_spec &arg_colour_spec ///< The broad_display_colour_spec to query
                                         ) {
	return sort_uniq_build<display_colour_vec>(
		arg_colour_spec.get_clr_of_pdb() | map_values
	);
}

/// \brief Get a list of the indices of all PDBs of the specified colour
///
/// relates broad_display_colour_spec
size_vec cath::get_pdbs_of_colour(const broad_display_colour_spec &arg_colour_spec, ///< The broad_display_colour_spec to query
                                  const display_colour            &arg_colour       ///< The colour to query
                                  ) {
	return copy_build<size_vec>(
		arg_colour_spec.get_clr_of_pdb()
		| filtered( [&] (const size_display_colour_map_val &x) { return ( x.second == arg_colour ); } )
		| map_keys
	);
}

/// \brief Write broad-level (ie not residue-specifiec) colouring instructions for the specified viewer to the specified ostream
///        using the specified broad_display_colour_spec and (viewer-cleaned) names
///        in the context of the specified list of colours
///
/// relates broad_display_colour_spec
void cath::detail::colour_base_and_pdbs_impl(const display_colour_vec        &arg_colours,                  ///< The list of colours
                                             const broad_display_colour_spec &arg_colour_spec,              ///< The broad_display_colour_spec to use
                                             const viewer                    &arg_viewer,                   ///< The viewer to specify how to render the instructions
                                             const pdb_list                  &/*arg_pdbs*/,                 ///< The key regions of the structures
                                             const str_vec                   &arg_cleaned_names_for_viewer, ///< The (viewer-cleaned) names for the structures to be coloured
                                             ostream                         &arg_os                        ///< The ostream to which the instructions should be written
                                             ) {
	const size_t num_colours = arg_colours.size();

	for (const size_t &colour_ctr : irange( 0_z, num_colours ) ) {
		arg_viewer.define_colour(
			arg_os,
			arg_colours[ colour_ctr ],
			generate_colour_name( colour_ctr, num_colours )
		);
	}

	if ( has_base_colour( arg_colour_spec ) ) {
		const string base_colour_name = "base_colour";
		arg_viewer.define_colour( arg_os, get_base_colour( arg_colour_spec ), base_colour_name );
		arg_os << arg_viewer.get_colour_base_str( base_colour_name );
	}

	for (const size_t colour_ctr : irange( 0_z, num_colours ) ) {
		const size_vec pdb_indices = get_pdbs_of_colour( arg_colour_spec, arg_colours[ colour_ctr ] );
		for (const size_t &pdb_index : pdb_indices) {
			arg_os << arg_viewer.get_colour_pdb_str(
				generate_colour_name( colour_ctr, num_colours ),
				arg_cleaned_names_for_viewer[ pdb_index  ]
			);
		}
	}
}

/// \brief Write broad-level (ie not residue-specifiec) colouring instructions for the specified viewer to the specified ostream
///        using the specified broad_display_colour_spec and (viewer-cleaned) names
///
/// relates broad_display_colour_spec
void cath::colour_viewer_with_spec(const broad_display_colour_spec &arg_broad_spec,               ///< The specification of colours to use
                                   const viewer                    &arg_viewer,                   ///< The viewer to specify how to render the instructions
                                   const pdb_list                  &arg_pdbs,                     ///< The key regions of the structures
                                   const str_vec                   &arg_cleaned_names_for_viewer, ///< The (viewer-cleaned) names for the structures to be coloured
                                   ostream                         &arg_os                        ///< The ostream to which the instructions should be written
                                   ) {
	detail::colour_base_and_pdbs_impl(
		get_pdb_colours( arg_broad_spec ),
		arg_broad_spec,
		arg_viewer,
		arg_pdbs,
		arg_cleaned_names_for_viewer,
		arg_os
	);
}
