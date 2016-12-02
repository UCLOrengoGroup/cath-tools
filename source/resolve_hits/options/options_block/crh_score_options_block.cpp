/// \file
/// \brief The crh_score_options_block class definitions

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

#include "crh_score_options_block.hpp"

#include <boost/format.hpp>

#include "common/clone/make_uptr_clone.hpp"
#include "common/program_options/prog_opt_num_range.hpp"

using namespace cath::common;
using namespace cath::opts;
using namespace cath::rslv;
using namespace cath;

using boost::format;
using boost::none;
using boost::program_options::bool_switch;
using boost::program_options::options_description;
using boost::program_options::value;
using boost::program_options::variables_map;
using std::string;
using std::unique_ptr;

/// \brief The option name for the degree to which long domains are preferred
const string crh_score_options_block::PO_LONG_DOMAINS_PREFERENCE { "long-domains-preference" };

/// \brief The option name for the degree to which high scores are preferred
const string crh_score_options_block::PO_HIGH_SCORES_PREFERENCE  { "high-scores-preference"  };

/// \brief The option name for whether to apply rules specific to CATH-Gene3D
const string crh_score_options_block::PO_APPLY_CATH_RULES        { "apply-cath-rules"        };

/// \brief A standard do_clone method
unique_ptr<options_block> crh_score_options_block::do_clone() const {
	return { make_uptr_clone( *this ) };
}

/// \brief Define this block's name (used as a header for the block in the usage)
string crh_score_options_block::do_get_block_name() const {
	return "Hit preference";
}

/// \brief Add this block's options to the provided options_description
void crh_score_options_block::do_add_visible_options_to_description(options_description &arg_desc ///< The options_description to which the options are added
                                                                    ) {
	const string val_varname { "<val>" };

	const auto long_domains_preference_notifier = [&] (const resscr_t &x) { the_spec.set_long_domains_preference ( x ); };
	const auto high_scores_preference_notifier  = [&] (const resscr_t &x) { the_spec.set_high_scores_preference  ( x ); };
	const auto apply_cath_rules_notifier        = [&] (const bool     &x) { the_spec.set_apply_cath_rules        ( x ); };

	arg_desc.add_options()
		(
			( PO_LONG_DOMAINS_PREFERENCE ).c_str(),
			value< prog_opt_num_range<resscr_t, -100, 100 > >()
				->value_name   ( val_varname                                     )
				->notifier     ( long_domains_preference_notifier                )
				->default_value( crh_score_spec::DEFAULT_LONG_DOMAINS_PREFERENCE ),
			( "Prefer longer hits to degree " + val_varname
				+ "\n(" + val_varname + " may be negative to prefer shorter; 0 leaves scores unaffected)").c_str()
		)
		(
			( PO_HIGH_SCORES_PREFERENCE ).c_str(),
			value< prog_opt_num_range<resscr_t, -100, 100 > >()
				->value_name   ( val_varname                                     )
				->notifier     ( high_scores_preference_notifier                 )
				->default_value(
					crh_score_spec::DEFAULT_HIGH_SCORES_PREFERENCE,
					( format( "%.4g" ) % crh_score_spec::DEFAULT_HIGH_SCORES_PREFERENCE ).str()
				),
			( "Prefer higher scores to degree " + val_varname
				+ "\n(" + val_varname + " may be negative to reduce preference for higher scores; 0 leaves scores unaffected)" ).c_str()
		)
		(
			( PO_APPLY_CATH_RULES ).c_str(),
			bool_switch()
				->notifier     ( apply_cath_rules_notifier                       )
				->default_value( crh_score_spec::DEFAULT_APPLY_CATH_RULES        ),
			"Apply rules specific to CATH-Gene3D during the parsing and processing"
		);

	static_assert( ! crh_score_spec::DEFAULT_APPLY_CATH_RULES, "If crh_score_spec::DEFAULT_APPLY_CATH_RULES isn't false, it might mess up the bool switch in here" );
}

/// \brief Generate a description of any problem that makes the specified crh_score_options_block invalid
///        or none otherwise
str_opt crh_score_options_block::do_invalid_string(const variables_map &/*arg_variables_map*/ ///< The variables map, which options_blocks can use to determine which options were specified, defaulted etc
                                                   ) const {
	return none;
}

/// \brief Getter for the crh_score_spec that the crh_score_options_block configures
const crh_score_spec & crh_score_options_block::get_crh_score_spec() const {
	return the_spec;
}
