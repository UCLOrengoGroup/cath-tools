/// \file
/// \brief The cath_hit_resolver class definitions

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

#include "cath_hit_resolver.hpp"

#include "common/file/ofstream_list.hpp"
#include "common/file/open_fstream.hpp"
#include "common/logger.hpp"
#include "exception/not_implemented_exception.hpp"
#include "exception/out_of_range_exception.hpp"
#include "resolve_hits/calc_hit_list.hpp"
#include "resolve_hits/file/hmmer_hmmsearch_domtblout.hpp"
#include "resolve_hits/file/hmmer_hmmsearch_out.hpp"
#include "resolve_hits/html_output/resolve_hits_html_outputter.hpp"
#include "resolve_hits/options/crh_options.hpp"
#include "resolve_hits/options/spec/crh_score_spec.hpp"
#include "resolve_hits/options/spec/crh_spec.hpp"
#include "resolve_hits/read_and_process_hits/read_and_process_mgr.hpp"

#include <fstream>

using namespace cath::common;
using namespace cath::opts;
using namespace cath::rslv;
using namespace std::literals::string_literals;

using std::istream;
using std::ifstream;
using std::ostream;
using std::ofstream;

/// \brief Perform resolve-hits according to the specified arguments strings with the specified i/o streams
void cath::rslv::perform_resolve_hits(const str_vec &args,        ///< The arguments strings specifying the resolve-hits action to perform
                                      istream       &arg_istream, ///< The input stream
                                      ostream       &arg_stdout   ///< The output stream
                                      ) {
	perform_resolve_hits(
		make_and_parse_options<crh_options>( args ),
		arg_istream,
		arg_stdout
	);
}

/// \brief Perform resolve-hits according to the specified crh_options with the specified i/o streams
void cath::rslv::perform_resolve_hits(const crh_options &arg_opts,    ///< The crh_options specifying the resolve-hits action to perform
                                      istream           &arg_istream, ///< The input stream
                                      ostream           &arg_stdout   ///< The output stream
                                      ) {
	// If the options are invalid or specify to do_nothing, then just return
	const auto &error_or_help_string = arg_opts.get_error_or_help_string();
	if ( error_or_help_string ) {
		arg_stdout << *error_or_help_string << "\n";
		return;
	}

	perform_resolve_hits(
		arg_opts.get_crh_spec(),
		arg_istream,
		arg_stdout
	);
}

/// \brief Perform resolve-hits according to the specified crh_spec with the specified i/o streams
void cath::rslv::perform_resolve_hits(const crh_spec &arg_crh_spec, ///< The crh_input_spec specifying the resolve-hits action to perform
                                      istream        &arg_istream,  ///< The input stream
                                      ostream        &arg_stdout    ///< The output stream
                                      ) {
	const auto &in_spec         = arg_crh_spec.get_input_spec();
	const auto &out_spec        = arg_crh_spec.get_output_spec();
	const auto &score_spec      = arg_crh_spec.get_score_spec();
	const auto &css_file_opt    = out_spec.get_export_css_file();
	const auto &input_file_opt  = in_spec.get_input_file();
	const auto &read_from_stdin = in_spec.get_read_from_stdin();

	// If CSS requested, export it
	if ( css_file_opt ) {
		ofstream_list oftreams{ arg_stdout };
		const auto ostream_refs = oftreams.open_ofstreams( { *css_file_opt } );
		if ( ostream_refs.size() != 1 ) {
			BOOST_THROW_EXCEPTION(out_of_range_exception("argh"));
		}
		ostream_refs.front().get() << resolve_hits_html_outputter::css_string();
		oftreams.close_all();
	}

	// If no input specified, stop here
	if ( ! input_file_opt && ! read_from_stdin ) {
		return;
	}

	// Organise the input stream
	ifstream input_file_stream;
	if ( input_file_opt ) {
		if ( ! exists( *input_file_opt ) ) {
			logger::log_and_exit(
				logger::return_code::NO_SUCH_FILE,
				"No such resolve-hits input data file \"" + input_file_opt->string() + "\""
			);
		}
		open_ifstream( input_file_stream, *input_file_opt );
	}
	istream &the_istream_ref = ( read_from_stdin ? arg_istream : input_file_stream );

	// Prepare a read_and_process_mgr object
	ofstream_list oftreams{ arg_stdout };
	read_and_process_mgr the_read_and_process_mgr = make_read_and_process_mgr(
		oftreams,
		arg_crh_spec
	);

	try {
		switch( in_spec.get_input_format() ) {
			case ( hits_input_format_tag::HMMER_DOMTBLOUT ) : {
				parse_domain_hits_table(
					the_read_and_process_mgr,
					the_istream_ref,
					score_spec.get_apply_cath_rules()
				);
				break;
			}
			case ( hits_input_format_tag::HMMSEARCH_OUT ) : {
				parse_hmmsearch_out(
					the_read_and_process_mgr,
					the_istream_ref,
					score_spec.get_apply_cath_rules(),
					in_spec.get_min_gap_length(),
					out_spec.get_output_hmmsearch_aln()
				);
				break;
			}
			case ( hits_input_format_tag::RAW_WITH_SCORES ) : {
				read_hit_list_from_istream(
					the_read_and_process_mgr,
					the_istream_ref,
					hit_score_type::CRH_SCORE
				);
				break;
			}
			case ( hits_input_format_tag::RAW_WITH_EVALUES ) : {
				read_hit_list_from_istream(
					the_read_and_process_mgr,
					the_istream_ref,
					hit_score_type::FULL_EVALUE
				);
				break;
			}
			default : {
				BOOST_THROW_EXCEPTION(out_of_range_exception("Value of hits_input_format_tag not recognised"));
			}
		}
	}
	catch (const std::exception &arg_exception) {
		logger::log_and_exit(
			logger::return_code::MALFORMED_RESOLVE_HITS_INFILE,
			"Unable to parse/process resolve-hits input data file \""
				+ ( input_file_opt ? input_file_opt->string() : "<stdin>"s )
				+ "\" of format "
				+ to_string( in_spec.get_input_format() )
				+ ". Error was:\n"
				+ arg_exception.what()
		);
	}

	// Close any open file streams
	oftreams.close_all();
	if ( input_file_opt ) {
		input_file_stream.close();
	}
}
