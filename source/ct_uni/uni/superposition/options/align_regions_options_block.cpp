/// \file
/// \brief The align_regions_options_block class definitions

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

#include "align_regions_options_block.hpp"

#include <boost/optional.hpp>

#include "chopping/domain/domain.hpp"
#include "chopping/region/region.hpp"
#include "common/clone/make_uptr_clone.hpp"

using namespace cath;
using namespace cath::chop;
using namespace cath::common;
using namespace cath::opts;

using boost::none;
using boost::program_options::options_description;
using boost::program_options::value;
using boost::program_options::variables_map;
using std::string;
using std::unique_ptr;

/// \brief The option name for the regions of the alignment (either to align or that have been aligned)
const string align_regions_options_block::PO_ALN_REGIONS { "align-regions"   };

/// \brief A standard do_clone method
unique_ptr<options_block> align_regions_options_block::do_clone() const {
	return { make_uptr_clone( *this ) };
}

/// \brief Define this block's name (used as a header for the block in the usage)
string align_regions_options_block::do_get_block_name() const {
	return "Regions";
}

/// \brief Add this block's options to the provided options_description
void align_regions_options_block::do_add_visible_options_to_description(options_description &prm_desc,           ///< The options_description to which the options are added
                                                                        const size_t        &/*prm_line_length*/ ///< The line length to be used when outputting the description (not very clearly documented in Boost)
                                                                        ) {
	const string regions_varname { "<regions>" };

	const auto align_regions_notifier = [&] (const domain_vec &x) { align_domains = x; };

	prm_desc.add_options()
		(
			( PO_ALN_REGIONS ).c_str(),
			value<domain_vec>()
				->value_name   ( regions_varname        )
				->notifier     ( align_regions_notifier ),
			( "Handle region(s) " + regions_varname + " as the alignment part of the structure.\n"
				+ "May be specified multiple times, in correspondence with the structures.\n"
				+ "Format is: D[5inwB02]251-348:B,408-416A:B\n"
				+ "(Put " + regions_varname + R"( in quotes to prevent the square brackets confusing your shell ("No match")))").c_str()
		);
}

/// \brief Generate a description of any problem that makes the specified align_regions_options_block invalid
///        or none otherwise
str_opt align_regions_options_block::do_invalid_string(const variables_map &/*prm_variables_map*/ ///< The variables map, which options_blocks can use to determine which options were specified, defaulted etc
                                                       ) const {
	return none;
}

/// \brief Return all options names for this block
str_vec align_regions_options_block::do_get_all_options_names() const {
	return {
		align_regions_options_block::PO_ALN_REGIONS,
	};
}

/// \brief Getter for the align regions that the align_regions_options_block configures
const domain_vec & align_regions_options_block::get_align_domains() const {
	return align_domains;
}
