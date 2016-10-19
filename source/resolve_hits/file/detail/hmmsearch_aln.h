/// \file
/// \brief The hmmsearch_aln class header

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

#ifndef HMMSEARCH_ALN_H_INCLUDED
#define HMMSEARCH_ALN_H_INCLUDED

#include <boost/optional.hpp>
#include <boost/range/combine.hpp>
#include <boost/utility/string_ref.hpp>

#include "common/boost_addenda/make_string_ref.h"
#include "common/string/string_parse_tools.h"
#include "exception/runtime_error_exception.h"
#include "resolve_hits/hit_seg.h"
#include "resolve_hits/res_arrow.h"
#include "resolve_hits/resolve_hits_type_aliases.h"

#include <string>

namespace cath {
	namespace rslv {
		namespace detail {

			/// \brief Represent whether a certain alignment position (from hmmsearch output) contains entries from
			///        the first sequence, the second or both
			enum class aln_presence {
				A,   ///< Contains only an entry from the first  sequence
				B,   ///< Contains only an entry from the second sequence
				BOTH ///< Contains entries from both sequences
			};

			/// \brief Type alias for a pair of aln_presence and a residue index type
			///
			/// This is used to represent a stretch of the alignment that has the same
			/// aln_presence; the second field indicates the length of the stretch
			using aln_stretch     = std::pair<aln_presence, residx_t>;

			/// \brief Type alias for a vector of aln_stretch types
			using aln_stretch_vec = std::vector<aln_stretch>;

			/// \brief Represent an alignment parsed from an hmmsearch output
			///
			/// Parsing the alignment is made more difficult because:
			///  * the alignment can be split over multiple blocks
			///  * the alignment doesn't necessarily start at 0 in either sequence
			///
			/// For example:
			///
			/// ~~~~~
			///                                        xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx RF
			///                    1aqkL01_round_1   3 ppsvsvspgesvtlsCtasssnigssyylhWyqqkpgqapklliysasnrasgvpdrfsgsksgtsatLtisslqae 79 
			///                                        + sv+v++g+++ l+C +s+++     ++ W++   +         ++n++ +v+   +++ +g sa+L+i+ + a 
			///   404421f6d00e5b5f33713585fb60b9bf 811 SRSVMVKKGDTALLQCAVSGDK---PINIVWMRSGKNTLN-----PSTNYKISVK--QEATPDGVSAELQIRTVDAT 877
			///                                        5799***************986...46899*998666544.....3444443333..34455677789********* PP
			///
			///                                        xxxxxxxxx RF
			///                    1aqkL01_round_1  80 DeAvYyCas 88 
			///                                        D++ Y+C +
			///   404421f6d00e5b5f33713585fb60b9bf 878 DSGPYFCRA 886
			///                                        *******76 PP
			/// ~~~~~
			class hmmsearch_aln final {
			private:
				/// \brief The ID of the first protein
				std::string id_a;

				/// \brief The ID of the second protein
				std::string id_b;

				/// \brief The (optional) index of the start of the alignment in the first sequence
				residx_opt start_a;

				/// \brief The (optional) index of the start of the alignment in the second sequence
				residx_opt start_b;

				/// \brief The string of the aligned first sequence
				std::string aln_seq_a;

				/// \brief A list of the stretches (each being an aln_stretch and associated length)
				aln_stretch_vec stretches;

				/// \brief A list of the segments
				hit_seg_vec segments;

				/// \brief The field offset to the ID in an alignment line
				static constexpr size_t LINE_ID_OFFSET    = 0;

				/// \brief The field offset to the start in an alignment line
				static constexpr size_t LINE_START_OFFSET = 1;

				/// \brief The field offset to the alignment in an alignment line
				static constexpr size_t LINE_ALN_OFFSET   = 2;

			public:
				void reset();

				void add_a(const std::string &);

				void add_b(const std::string &);

				std::pair<std::string &, hit_seg_vec &> process_aln(const residx_t &);

				bool empty() const;
			};

			/// \brief Reset this hmmsearch_aln to be reused
			inline void hmmsearch_aln::reset() {
				id_a.clear();
				id_b.clear();
				start_a = boost::none;
				start_b = boost::none;
				aln_seq_a.clear();
				stretches.clear();
			}

			/// \brief Add the data from a first alignment line
			///
			/// This should always be called before add_b() and after reset() or add_b() 
			inline void hmmsearch_aln::add_a(const std::string &arg_line ///< The first alignment line (eg "                   1aqkL01_round_1   3 ppsvsvspgesvtlsCtasssnigssyylhWyqqkpgqapklliysasnrasgvpdrfsgsksgtsatLtisslqae 79 ")
			                                 ) {
				const auto id_itrs    = common::find_field_itrs( arg_line, LINE_ID_OFFSET                                              );
				const auto start_itrs = common::find_field_itrs( arg_line, LINE_START_OFFSET, 1 + LINE_ID_OFFSET,    id_itrs.second    );
				const auto aln_itrs   = common::find_field_itrs( arg_line, LINE_ALN_OFFSET,   1 + LINE_START_OFFSET, start_itrs.second );

				if ( id_a.empty() ) {
					id_a.assign( id_itrs.first, id_itrs.second );
				}
				if ( ! start_a ) {
					start_a = common::parse_uint_from_field( start_itrs.first, start_itrs.second );
				}
				aln_seq_a.assign( aln_itrs.first, aln_itrs.second );
			}

			/// \brief Add the data from a second alignment line
			///
			/// This should always be called after add_a() and before add_a() or process_aln()
			inline void hmmsearch_aln::add_b(const std::string &arg_line ///< The second alignment line (eg "  404421f6d00e5b5f33713585fb60b9bf 811 SRSVMVKKGDTALLQCAVSGDK---PINIVWMRSGKNTLN-----PSTNYKISVK--QEATPDGVSAELQIRTVDAT 877")
			                                 ) {
				const auto id_itrs    = common::find_field_itrs( arg_line, LINE_ID_OFFSET                                              );
				const auto start_itrs = common::find_field_itrs( arg_line, LINE_START_OFFSET, 1 + LINE_ID_OFFSET,    id_itrs.second    );
				const auto aln_itrs   = common::find_field_itrs( arg_line, LINE_ALN_OFFSET,   1 + LINE_START_OFFSET, start_itrs.second );

				if ( id_b.empty() ) {
					id_b.assign( id_itrs.first, id_itrs.second );
				}
				if ( ! start_b ) {
					start_b = common::parse_uint_from_field( start_itrs.first, start_itrs.second );
				}
				assert( start_a );

				const auto the_sequence = common::make_string_ref( aln_itrs.first, aln_itrs.second );
				for (const boost::tuple<const char &, const char &> &x : boost::range::combine( aln_seq_a, the_sequence ) ) {
					const char &char_a = x.get<0>();
					const char &char_b = x.get<1>();

					assert( char_b != '-' || char_a != '.' );
					if ( char_b == '-' && char_a == '.' ) {
						BOOST_THROW_EXCEPTION(common::runtime_error_exception("hmmsearch aligment files shouldn't contain alignment positions that have gaps on both sides"));
					}

					const aln_presence presence = ( char_a == '.' ) ? aln_presence::B :
					                              ( char_b == '-' ) ? aln_presence::A :
					                                                  aln_presence::BOTH;

					if ( stretches.empty() || stretches.back().first != presence ) {
						stretches.emplace_back( presence, 1 );
					}
					else {
						++( stretches.back().second );
					}
				}
			}

			/// \brief Process the alignment that has been loaded into this hmmsearch_aln
			///
			/// This should always be called after add_b()
			inline std::pair<std::string &, hit_seg_vec &> hmmsearch_aln::process_aln(const residx_t &arg_min_gap_length ///< The minimum length for a gap to be considered a gap
			                                                                          ) {
				segments.clear();
				res_arrow arrow_b = arrow_before_res( *start_b );
				for (const aln_stretch &stretch : stretches) {
					const auto &presence = stretch.first;
					const auto &length   = stretch.second;

					if ( presence == aln_presence::BOTH || presence == aln_presence::B ) {
						const res_arrow  start = arrow_b;
						arrow_b += length;
						const res_arrow &stop  = arrow_b;

						if ( presence == aln_presence::BOTH || length < arg_min_gap_length ) {
							if ( segments.empty() || segments.back().get_stop_arrow() != start ) {
								segments.emplace_back( start, stop );
							}
							else {
								segments.back().set_stop_arrow( stop );
							}
						}
					}
				}
				return { id_a, segments };
			}

			/// \brief Return whether the stretches are empty (ie none has been parsed in)
			inline bool hmmsearch_aln::empty() const {
				return stretches.empty();
			}

		}
	}
}

#endif
