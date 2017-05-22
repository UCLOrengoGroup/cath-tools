/// \file
/// \brief The crh_spec class definitions

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

#include "crh_spec.hpp"

using namespace cath::rslv;

/// \brief Ctor from the individual specs
crh_spec::crh_spec(crh_input_spec         arg_input_spec,         ///< The input spec
                   crh_segment_spec       arg_segment_spec,       ///< The segment spec
                   crh_score_spec         arg_score_spec,         ///< The score spec
                   crh_filter_spec        arg_filter_spec,        ///< The filter spec
                   crh_single_output_spec arg_single_output_spec, ///< The output spec
                   crh_html_spec          arg_html_spec           ///< The html spec
                   ) : the_input_spec         { std::move( arg_input_spec         ) },
                       the_segment_spec       { std::move( arg_segment_spec       ) },
                       the_score_spec         { std::move( arg_score_spec         ) },
                       the_filter_spec        { std::move( arg_filter_spec        ) },
                       the_single_output_spec { std::move( arg_single_output_spec ) },
                       the_html_spec          { std::move( arg_html_spec          ) } {
}

/// \brief Non-const overload of getter for input_spec
crh_input_spec & crh_spec::get_input_spec() {
	return the_input_spec;
}

/// \brief Const overload of getter for input_spec
const crh_input_spec & crh_spec::get_input_spec() const {
	return the_input_spec;
}

/// \brief Non-const overload of getter for segment_spec
crh_segment_spec & crh_spec::get_segment_spec() {
	return the_segment_spec;
}

/// \brief Const overload of getter for segment_spec
const crh_segment_spec & crh_spec::get_segment_spec() const {
	return the_segment_spec;
}

/// \brief Non-const overload of getter for score_spec
crh_score_spec & crh_spec::get_score_spec() {
	return the_score_spec;
}

/// \brief Const overload of getter for score_spec
const crh_score_spec & crh_spec::get_score_spec() const {
	return the_score_spec;
}

/// \brief Non-const overload of getter for filter_spec
crh_filter_spec & crh_spec::get_filter_spec() {
	return the_filter_spec;
}

/// \brief Const overload of getter for filter_spec
const crh_filter_spec & crh_spec::get_filter_spec() const {
	return the_filter_spec;
}

/// \brief Non-const overload of getter for output_spec
crh_single_output_spec & crh_spec::get_single_output_spec() {
	return the_single_output_spec;
}

/// \brief Const overload of getter for output_spec
const crh_single_output_spec & crh_spec::get_single_output_spec() const {
	return the_single_output_spec;
}

/// \brief Non-const overload of getter for html_spec
crh_html_spec & crh_spec::get_html_spec() {
	return the_html_spec;
}

/// \brief Const overload of getter for html_spec
const crh_html_spec & crh_spec::get_html_spec() const {
	return the_html_spec;
}

/// \brief Setter for the input spec
crh_spec & crh_spec::set_input_spec(const crh_input_spec &arg_input_spec ///< The input spec
                                    ) {
	get_input_spec() = arg_input_spec;
	return *this;
}

/// \brief Setter for the segment spec
crh_spec & crh_spec::set_segment_spec(const crh_segment_spec &arg_segment_spec ///< The segment spec
                                      ) {
	get_segment_spec() = arg_segment_spec;
	return *this;
}

/// \brief Setter for the score spec
crh_spec & crh_spec::set_score_spec(const crh_score_spec &arg_score_spec ///< The score spec
                                    ) {
	get_score_spec() = arg_score_spec;
	return *this;
}

/// \brief Setter for the filter spec
crh_spec & crh_spec::set_filter_spec(const crh_filter_spec &arg_filter_spec ///< The filter spec
                                     ) {
	get_filter_spec() = arg_filter_spec;
	return *this;
}

/// \brief Setter for the output spec
crh_spec & crh_spec::set_single_output_spec(const crh_single_output_spec &arg_single_output_spec ///< The output spec
                                            ) {
	get_single_output_spec() = arg_single_output_spec;
	return *this;
}

/// \brief Setter for the html spec
crh_spec & crh_spec::set_html_spec(const crh_html_spec &arg_html_spec ///< The html spec
                                   ) {
	get_html_spec() = arg_html_spec;
	return *this;
}
