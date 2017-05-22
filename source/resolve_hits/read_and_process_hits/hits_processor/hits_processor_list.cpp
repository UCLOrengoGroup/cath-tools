/// \file
/// \brief The hits_processor_list class definitions

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

#include "hits_processor_list.hpp"

#include <boost/algorithm/cxx11/any_of.hpp>

#include "resolve_hits/options/spec/crh_single_output_spec.hpp"
#include "resolve_hits/read_and_process_hits/hits_processor/summarise_hits_processor.hpp"
#include "resolve_hits/read_and_process_hits/hits_processor/write_html_hits_processor.hpp"
#include "resolve_hits/read_and_process_hits/hits_processor/write_json_hits_processor.hpp"
#include "resolve_hits/read_and_process_hits/hits_processor/write_results_hits_processor.hpp"

using namespace cath::common;
using namespace cath::rslv;
using namespace cath::rslv::detail;

using boost::algorithm::any_of;
using std::initializer_list;
using std::make_unique;
using std::ostream;
using std::string;
using std::unique_ptr;

/// \brief hits_processor_list ctor from crh_score_spec and crh_segment_spec
hits_processor_list::hits_processor_list(const crh_score_spec                   &arg_crh_score_spec,   ///< The crh_segment_spec with which to construct this hits_processor_list
                                         const crh_segment_spec                 &arg_crh_segment_spec, ///< The crh_segment_spec with which to construct this hits_processor_list
                                         initializer_list<hits_processor_clptr>  arg_hits_processors   ///< Any hits_processor_clprts with which this hits_processor_list should be initialized
                                         ) : processors          { arg_hits_processors  },
                                             the_score_spec      { arg_crh_score_spec   },
                                             the_segment_spec    { arg_crh_segment_spec } {
}

/// \brief Getter for the score spec to apply to incoming hits
const crh_score_spec & hits_processor_list::get_score_spec() const {
	return the_score_spec;
}

/// \brief Getter for the segment spec to apply to incoming hits
const crh_segment_spec & hits_processor_list::get_segment_spec() const {
	return the_segment_spec;
}

/// \brief Add a processor to the list
hits_processor_list & hits_processor_list::add_processor(const hits_processor &arg_hits_processor ///< The processor to add
                                                         ) {
	processors.push_back( arg_hits_processor.clone() );
	return *this;
}

/// \brief Add a processor to the list
hits_processor_list & hits_processor_list::add_processor(hits_processor_uptr arg_hits_processor_uptr ///< (A pointer to) the processor to add
                                                         ) {
	processors.push_back( std::move( arg_hits_processor_uptr ) );
	return *this;
}

/// \brief Add a processor to the list
hits_processor_list & hits_processor_list::add_processor(hits_processor_clptr arg_hits_processor_uptr ///< (A pointer to) the processor to add
                                                         ) {
	processors.push_back( std::move( arg_hits_processor_uptr ) );
	return *this;
}

/// \brief Return whether the list of hits_processors is empty
bool hits_processor_list::empty() const {
	return processors.empty();
}

/// \brief Return the number ofhits_processors in the list
size_t hits_processor_list::size() const {
	return processors.size();
}

/// \brief Access the hits_processor at the specified index
///
/// This isn't bounds-checked (unless in a build with a bounds-checked STL)
const hits_processor & hits_processor_list::operator[](const size_t &arg_index ///< The index of the hits_processor to access
                                                       ) const {
	return *( processors[ arg_index ] );
}

/// \brief Return whether any of the hits_processors in the list want to hear about hits that fail the score_filter
bool hits_processor_list::wants_hits_that_fail_score_filter() const {
	return any_of( *this, [] (const hits_processor &x) { return x.wants_hits_that_fail_score_filter(); } );
}

/// \brief Standard const begin() method, as part of making this a range over hits_processors
///
/// Note that this pipes through boost::indirected_range so the range is
/// over elements of type `const hits_processor &` rather than `clone_ptr<const hits_processor>`
auto hits_processor_list::begin() const -> const_iterator {
	return common::cbegin( processors | boost::adaptors::indirected );
}

/// \brief Standard const end() method, as part of making this a range over hits_processors
///
/// Note that this pipes through boost::indirected_range so the range is
/// over elements of type `const hits_processor &` rather than `clone_ptr<const hits_processor>`
auto hits_processor_list::end() const -> const_iterator {
	return common::cend( processors | boost::adaptors::indirected );
}

/// \brief Make the hits_processor implied by the specified spec object and the specified ostream
///
/// \relates hits_processor_list
hits_processor_list cath::rslv::detail::make_hits_processors(ostream                      &arg_ostream,            ///< The ostream to which the new hits_processor should write
                                                             const crh_single_output_spec &arg_single_output_spec, ///< The crh_single_output_spec defining the type of hits_processor to make
                                                             const crh_score_spec         &arg_score_spec,         ///< The crh_score_spec how to handle scores
                                                             const crh_segment_spec       &arg_segment_spec,       ///< The crh_segment_spec how to handle segments
                                                             const crh_html_spec          &arg_html_spec           ///< The crh_html_spec defining how to render any HTML
                                                             ) {
	const auto &bound_out = arg_single_output_spec.get_boundary_output();
	hits_processor_list the_list{ arg_score_spec, arg_segment_spec };
	the_list.add_processor( [&] () -> unique_ptr<hits_processor> {
		switch ( get_out_format( arg_single_output_spec ) ) {
			case ( crh_out_format::HTML     ) : { return make_unique< write_html_hits_processor    >( arg_ostream, arg_html_spec ); }
			case ( crh_out_format::SUMMARY  ) : { return make_unique< summarise_hits_processor     >( arg_ostream                ); }
			case ( crh_out_format::STANDARD ) : { return make_unique< write_results_hits_processor >( arg_ostream, bound_out     ); }
			case ( crh_out_format::JSON     ) : { return make_unique< write_json_hits_processor    >( arg_ostream                ); }
		}
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Value of crh_out_format not recognised whilst converting to_string()"));
	} () );
	return the_list;
}
