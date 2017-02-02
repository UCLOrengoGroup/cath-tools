/// \file
/// \brief The crh_html_options_block class definitions

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

#include "crh_html_options_block.hpp"

#include <boost/algorithm/cxx11/any_of.hpp>

#include "common/clone/make_uptr_clone.hpp"
#include "common/program_options/prog_opt_num_range.hpp"

#include <limits>

using namespace cath;
using namespace cath::common;
using namespace cath::opts;
using namespace cath::rslv;

using boost::algorithm::any_of;
using boost::none;
using boost::program_options::bool_switch;
using boost::program_options::options_description;
using boost::program_options::value;
using boost::program_options::variables_map;
using std::numeric_limits;
using std::string;
using std::unique_ptr;

/// \brief  The option name for whether to restrict HTML output to the contents of the body tag
const string crh_html_options_block::PO_RESTRICT_HTML_WITHIN_BODY { "restrict-html-within-body"  };

/// \brief  The option name for the maximum number of non-solution hits to display in the HTML
const string crh_html_options_block::PO_MAX_NUM_NON_SOLN_HITS     { "html-max-num-non-soln-hits" };

/// \brief  The option name for whether to exclude hits rejected by the score filters from the HTML
const string crh_html_options_block::PO_EXCLUDE_REJECTED_HITS     { "html-exclude-rejected-hits" };

/// \brief A list of all the option names used in this options_block
const str_vec crh_html_options_block::ALL_BLOCK_POS = { crh_html_options_block::PO_RESTRICT_HTML_WITHIN_BODY,
                                                        crh_html_options_block::PO_MAX_NUM_NON_SOLN_HITS,
                                                        crh_html_options_block::PO_EXCLUDE_REJECTED_HITS };

/// \brief A standard do_clone method
unique_ptr<options_block> crh_html_options_block::do_clone() const {
	return { make_uptr_clone( *this ) };
}

/// \brief Define this block's name (used as a header for the block in the usage)
string crh_html_options_block::do_get_block_name() const {
	return "HTML";
}

/// \brief Add this block's options to the provided options_description
void crh_html_options_block::do_add_visible_options_to_description(options_description &arg_desc ///< The options_description to which the options are added
                                                                    ) {
	const string num_varname { "<num>" };

	const auto set_restrict_html_within_body_notifier = [&] (const bool   &x) { the_spec.set_restrict_html_within_body( x ); };
	const auto max_num_non_soln_hits_notifier         = [&] (const size_t &x) { the_spec.set_max_num_non_soln_hits    ( x ); };
	const auto exclude_rejected_hits_notifier         = [&] (const bool   &x) { the_spec.set_exclude_rejected_hits    ( x ); };

	arg_desc.add_options()
		(
			( PO_RESTRICT_HTML_WITHIN_BODY ).c_str(),
			bool_switch()
				->notifier     ( set_restrict_html_within_body_notifier           )
				->default_value( crh_html_spec::DEFAULT_RESTRICT_HTML_WITHIN_BODY ),
			"Restrict HTML output to the contents of the body tag.\n"
				"The contents should be included inside a body tag of class crh-body"
		)
		(
			( PO_MAX_NUM_NON_SOLN_HITS ).c_str(),
			value< prog_opt_num_range<size_t, 0, numeric_limits<uint32_t>::max(), int64_t> >()
				->value_name   ( num_varname                                      )
				->notifier     ( max_num_non_soln_hits_notifier                   )
				->default_value( crh_html_spec::DEFAULT_MAX_NUM_NON_SOLN_HITS     ),
			( "Only display up to " + num_varname + " non-solution hits in the HTML" ).c_str()
		)
		(
			( PO_EXCLUDE_REJECTED_HITS ).c_str(),
			bool_switch()
				->notifier     ( exclude_rejected_hits_notifier                   )
				->default_value( crh_html_spec::DEFAULT_EXCLUDE_REJECTED_HITS     ),
			"Exclude hits rejected by the score filters from the HTML"
		);

	static_assert( ! crh_html_spec::DEFAULT_RESTRICT_HTML_WITHIN_BODY,
		"If crh_output_spec::DEFAULT_RESTRICT_HTML_WITHIN_BODY isn't false, it might mess up the bool switch in here" );
	static_assert( ! crh_html_spec::DEFAULT_EXCLUDE_REJECTED_HITS,
		"If crh_output_spec::DEFAULT_EXCLUDE_REJECTED_HITS isn't false, it might mess up the bool switch in here" );
}

/// \brief Generate a description of any problem that makes the specified crh_html_options_block invalid
///        or none otherwise
str_opt crh_html_options_block::do_invalid_string(const variables_map &/*arg_variables_map*/ ///< The variables map, which options_blocks can use to determine which options were specified, defaulted etc
                                                  ) const {
	return none;
}

/// \brief Getter for the crh_html_spec that the crh_html_options_block configures
const crh_html_spec & crh_html_options_block::get_crh_html_spec() const {
	return the_spec;
}

/// \brief Return whether the specified variables_map indicates any of the crh_html_options_block options have been specified
bool cath::rslv::has_specified_crh_html_options(const variables_map &arg_vm ///< The variables_map to query
                                                ) {
	return any_of(
		crh_html_options_block::ALL_BLOCK_POS,
		[&] (const string &x) { return ( arg_vm.count( x ) && ! arg_vm[ x ].defaulted() ); }
	);
}
