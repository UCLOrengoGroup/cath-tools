/// \file
/// \brief The hmmer_aln class header

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

#ifndef _CATH_TOOLS_SOURCE_RESOLVE_HITS_FILE_DETAIL_HMMER_ALN_HPP
#define _CATH_TOOLS_SOURCE_RESOLVE_HITS_FILE_DETAIL_HMMER_ALN_HPP

#include <boost/optional.hpp>
#include <boost/range/combine.hpp>
#include <boost/utility/string_ref.hpp>

#include "common/boost_addenda/make_string_ref.hpp"
#include "common/exception/runtime_error_exception.hpp"
#include "common/string/string_parse_tools.hpp"
#include "resolve_hits/file/alnd_rgn.hpp"
#include "resolve_hits/resolve_hits_type_aliases.hpp"
#include "seq/seq_arrow.hpp"
#include "seq/seq_seg.hpp"

#include <string>

namespace cath {
	namespace rslv {
		namespace detail {

			/// \brief Represent whether a certain alignment position (from hmmsearch/hmmscan output) contains entries from
			///        the first sequence, the second or both
			enum class aln_presence : char {
				A,   ///< Contains only an entry from the first  sequence
				B,   ///< Contains only an entry from the second sequence
				BOTH ///< Contains entries from both sequences
			};

			/// \brief Type alias for a pair of aln_presence and a residue index type
			///
			/// This is used to represent a stretch of the alignment that has the same
			/// aln_presence; the second field indicates the length of the stretch
			using aln_stretch     = std::pair<aln_presence, seq::residx_t>;

			/// \brief Type alias for a vector of aln_stretch types
			using aln_stretch_vec = std::vector<aln_stretch>;

			/// \brief Represent an alignment parsed from an hmmsearch/hmmscan output
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
			class hmmer_aln final {
			private:
				/// \brief The ID of the first protein
				///
				/// An empty string means "not yet been initialised" (because IDs must always have length >= 1)
				std::string id_a;

				/// \brief The ID of the second protein
				///
				/// An empty string means "not yet been initialised" (because IDs must always have length >= 1)
				std::string id_b;

				/// \brief The (optional) index of the start of the alignment in the first sequence
				///
				/// This is optional to allow it be initialised as the data is read in
				seq::residx_opt start_a;

				/// \brief The (optional) index of the start of the alignment in the second sequence
				///
				/// This is optional to allow it be initialised as the data is read in
				seq::residx_opt start_b;

				/// \brief The string of the aligned first sequence
				///
				/// This is stored during add_a() so that it's available during add_b()
				std::string aln_seq_a;

				/// \brief A list of the stretches (each being an aln_stretch and associated length)
				aln_stretch_vec stretches;

				/// \brief A list of the segments
				///
				/// This is only used in process_aln but is stored here to allow the
				/// vector's allocated memory to be reused between results
				seq::seq_seg_vec segments;

				/// \brief A list of the aligned regions
				///
				/// This is only used in process_aln but is stored here to allow the
				/// vector's allocated memory to be reused between results
				alnd_rgn_vec aligned_regions;

				/// \brief The field offset to the ID in an alignment line
				static constexpr size_t LINE_ID_OFFSET    = 0;

				/// \brief The field offset to the start in an alignment line
				static constexpr size_t LINE_START_OFFSET = 1;

				/// \brief The field offset to the alignment in an alignment line
				static constexpr size_t LINE_ALN_OFFSET   = 2;

				void update_segments(const seq::residx_t &);
				void update_aligned_regions(const bool &);

			public:
				void reset();

				void add_a(const std::string &);

				void add_b(const std::string &);

				std::tuple<std::string &, seq::seq_seg_vec &, alnd_rgn_vec &> process_aln(const seq::residx_t &,
				                                                                          const bool &);

				bool empty() const;
			};

			/// \brief Update the segments based on the parsed information
			inline void hmmer_aln::update_segments(const seq::residx_t &arg_min_gap_length ///< The minimum length for a missing stretch to be considered a gap
			                                       ) {
				segments.clear();
				seq::seq_arrow arrow_b = seq::arrow_before_res( *start_b );
				for (const aln_stretch &stretch : stretches) {
					const auto &presence = stretch.first;
					const auto &length   = stretch.second;

					if ( presence == aln_presence::BOTH || presence == aln_presence::B ) {
						const seq::seq_arrow  start = arrow_b;
						arrow_b += length;
						const seq::seq_arrow &stop  = arrow_b;

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
			}

			/// \brief Update the aligned regions based on the parsed information
			inline void hmmer_aln::update_aligned_regions(const bool &arg_parse_hmmer_aln ///< Whether to actually bother parsing this HMMER alignment information
			                                              ) {
				aligned_regions.clear();
				if ( arg_parse_hmmer_aln ) {
					seq::seq_arrow arrow_a = seq::arrow_before_res( *start_a );
					seq::seq_arrow arrow_b = seq::arrow_before_res( *start_b );
					for (const aln_stretch &stretch : stretches) {
						const auto &presence = stretch.first;
						const auto &length   = stretch.second;

						if ( presence == aln_presence::BOTH ) {
							aligned_regions.emplace_back( arrow_a, arrow_b, length );
						}

						if ( presence == aln_presence::BOTH || presence == aln_presence::A ) {
							arrow_a += length;
						}

						if ( presence == aln_presence::BOTH || presence == aln_presence::B ) {
							arrow_b += length;
						}
					}
				}
			}

			/// \brief Reset this hmmer_aln to be reused
			inline void hmmer_aln::reset() {
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
			inline void hmmer_aln::add_a(const std::string &arg_line ///< The first alignment line (eg "                   1aqkL01_round_1   3 ppsvsvspgesvtlsCtasssnigssyylhWyqqkpgqapklliysasnrasgvpdrfsgsksgtsatLtisslqae 79 ")
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
			inline void hmmer_aln::add_b(const std::string &arg_line ///< The second alignment line (eg "  404421f6d00e5b5f33713585fb60b9bf 811 SRSVMVKKGDTALLQCAVSGDK---PINIVWMRSGKNTLN-----PSTNYKISVK--QEATPDGVSAELQIRTVDAT 877")
			                             ) {
				constexpr char GAP_CHAR_A = '.';
				constexpr char GAP_CHAR_B = '-';

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

					assert( char_b != GAP_CHAR_B || char_a != GAP_CHAR_A );
					if ( char_b == GAP_CHAR_B && char_a == GAP_CHAR_A ) {
						BOOST_THROW_EXCEPTION(common::runtime_error_exception("HMMER alignment files shouldn't contain alignment positions that have gaps on both sides"));
					}

					const aln_presence presence = ( char_a == GAP_CHAR_A ) ? aln_presence::B :
					                              ( char_b == GAP_CHAR_B ) ? aln_presence::A :
					                                                         aln_presence::BOTH;

					if ( stretches.empty() || stretches.back().first != presence ) {
						stretches.emplace_back( presence, 1 );
					}
					else {
						++( stretches.back().second );
					}
				}
			}

			/// \brief Process the alignment that has been loaded into this hmmer_aln
			///
			/// This should always be called after add_b()
			inline std::tuple<std::string &, seq::seq_seg_vec &, alnd_rgn_vec &> hmmer_aln::process_aln(const seq::residx_t &arg_min_gap_length, ///< The minimum length for a gap to be considered a gap
			                                                                                            const bool          &arg_parse_hmmer_aln ///< Whether to parse the HMMER alignment information for outputting later
			                                                                                            ) {
				update_segments       ( arg_min_gap_length      );
				update_aligned_regions( arg_parse_hmmer_aln );
				return std::tie( id_a, segments, aligned_regions );
			}

			/// \brief Return whether the stretches are empty (ie none has been parsed in)
			inline bool hmmer_aln::empty() const {
				return stretches.empty();
			}

		} // namespace detail
	} // namespace rslv
} // namespace cath

#endif
