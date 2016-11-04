/// \file
/// \brief The trim_spec class header

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

#ifndef _CATH_TOOLS_SOURCE_RESOLVE_HITS_TRIM_TRIM_SPEC_H
#define _CATH_TOOLS_SOURCE_RESOLVE_HITS_TRIM_TRIM_SPEC_H

#include <boost/any.hpp>
#include <boost/optional/optional_fwd.hpp>

#include "common/debug_numeric_cast.h"
#include "common/type_aliases.h"
#include "resolve_hits/hit_seg.h"
#include "resolve_hits/resolve_hits_type_aliases.h"

#include <stdexcept>

namespace cath {
	namespace rslv {

		/// \brief Specify a trimming policy that can be applied to trim the edges off hits' segments
		///        during calculation so that they can overlap a bit
		///
		/// The policy is specified by stating a segment length and the total amount of trimming (split
		/// between the start and end) to be performed.
		///
		/// For longer segments, the total trim stays at the value specified here.
		/// For shorter segments, it decreases linearly down to a total trim of 0 for a segment of length 1
		class trim_spec final {
		private:
			/// \brief The length of a segment on which the trimming is to be defined
			residx_t full_length = DEFAULT_FULL_LENGTH;

			/// \brief The total amount of trimming (split between start and end) to be performed on a segment of the specified length
			residx_t total_trimming = DEFAULT_TOTAL_TRIMMING;

		public:
			/// \brief The default value for the length of a segment on which the trimming is to be defined
			static constexpr residx_t DEFAULT_FULL_LENGTH = 1;

			/// \brief The default value for the total amount of trimming (split between start and end) to be performed on a segment of the specified length
			static constexpr residx_t DEFAULT_TOTAL_TRIMMING = 0;

			constexpr trim_spec() noexcept = default;

			constexpr trim_spec(const residx_t &,
			                    const residx_t &);

			constexpr const residx_t & get_full_length() const noexcept;
			constexpr const residx_t & get_total_trimming() const noexcept;
		};
		

		std::string to_options_string(const trim_spec &);

		std::string to_string(const trim_spec &);

		std::ostream & operator<<(std::ostream &,
		                          const trim_spec &);

		trim_spec parse_trim_spec(const std::string &);

		std::istream & operator>>(std::istream &,
		                          trim_spec &);

		std::string to_possibly_trimmed_simple_string(const hit_seg &,
		                                              const boost::optional<trim_spec> &);

		std::string get_segments_string(const hit_seg_vec &,
		                                const boost::optional<trim_spec> &);

		void validate(boost::any &,
		              const str_vec &,
		              trim_spec*,
		              int);

		/// \brief Ctor for a trim_spec from the full_length and total_trimming
		constexpr trim_spec::trim_spec(const residx_t &arg_full_length,   ///< The length of a segment on which the trimming is to be defined
		                               const residx_t &arg_total_trimming ///< The total amount of trimming (split between start and end) to be performed on a segment of the specified length
		                               ) : full_length   { arg_full_length },
		                                   total_trimming{
		                                   	( arg_total_trimming < arg_full_length )
		                                   	? arg_total_trimming
		                                   	: throw "trim_spec total_trimming must be non-negative and less than the full_length"
		                                   } {
			static_assert( std::is_unsigned< std::decay_t< decltype( arg_full_length    ) > >::value, "If arg_full_length    isn't of an unsigned type, then add an explicit >= 0 check" );
			static_assert( std::is_unsigned< std::decay_t< decltype( arg_total_trimming ) > >::value, "If arg_total_trimming isn't of an unsigned type, then add an explicit >= 0 check" );
		}

		/// \brief Getter for the length of a segment on which the trimming is to be defined
		constexpr const residx_t & trim_spec::get_full_length() const noexcept {
			return full_length;
		}

		/// \brief Getter for the total amount of trimming (split between start and end) to be performed on a segment of the specified length
		constexpr const residx_t & trim_spec::get_total_trimming() const noexcept {
			return total_trimming;
		}

		/// \brief Make a trim_spec that implies no trimming
		constexpr trim_spec make_no_trim_trim_spec() {
			return { 1, 0 };
		}

		/// \brief Get the total trimming resulting from applying the specified trim_spec on a segment of the specified length
		///
		/// \relates trim_spec
		constexpr residx_t total_trimming_of_length(const trim_spec &arg_trim_spec, ///< The trim_spec to apply
		                                            const residx_t  &arg_length     ///< The length of the segment in question
		                                            ) {
			return
				( arg_length == 0 )
					? 0
					: ( arg_length >= arg_trim_spec.get_full_length() )
						? arg_trim_spec.get_total_trimming()
						: ( ( ( arg_length - 1 ) * arg_trim_spec.get_total_trimming() ) / ( arg_trim_spec.get_full_length() - 1 ) );
		}

		/// \brief Get the length-after-trimming resulting from applying the specified trim_spec on a segment of the specified length
		///
		/// \relates trim_spec
		constexpr residx_t length_after_trimming(const trim_spec &arg_trim_spec, ///< The trim_spec to apply
		                                         const residx_t  &arg_length     ///< The length of the segment in question
		                                         ) {
			return arg_length - total_trimming_of_length( arg_trim_spec, arg_length );
		}

		/// \brief Get the trimming to be applied at the start resulting from applying the specified trim_spec on a segment of the specified length
		///
		/// \relates trim_spec
		constexpr residx_t start_trimming_of_length(const trim_spec &arg_trim_spec, ///< The trim_spec to apply
		                                            const residx_t  &arg_length     ///< The length of the segment in question
		                                            ) {
			return total_trimming_of_length( arg_trim_spec, arg_length ) / 2;
		}

		/// \brief Get the trimming to be applied at the stop resulting from applying the specified trim_spec on a segment of the specified length
		///
		/// \relates trim_spec
		constexpr residx_t stop_trimming_of_length(const trim_spec &arg_trim_spec, ///< The trim_spec to apply
		                                           const residx_t  &arg_length     ///< The length of the segment in question
		                                           ) {
			return
				( total_trimming_of_length( arg_trim_spec, arg_length ) / 2 )
				+
				( total_trimming_of_length( arg_trim_spec, arg_length ) % 2 );
		}

		/// \brief Trim the specified hit_seg according to the specified trim_spec
		///
		/// \todo Come relaxed constexpr, make this constexpr
		///
		/// \relates trim_spec
		inline void trim_hit_seg(hit_seg         &arg_hit_seg,    ///< The hit_seg to trim
		                         const trim_spec &arg_trim_spec ///< The trim_spec to apply
		                         ) {
			const auto length = get_length( arg_hit_seg );
			arg_hit_seg.set_start_arrow( arg_hit_seg.get_start_arrow() + start_trimming_of_length( arg_trim_spec, length ) );
			arg_hit_seg.set_stop_arrow ( arg_hit_seg.get_stop_arrow () - stop_trimming_of_length ( arg_trim_spec, length ) );
		}

		/// \brief Return a copy of the specified hit_seg, trimmed according to the specified trim_spec
		///
		/// \todo Come relaxed constexpr, make this constexpr
		///
		/// \relates trim_spec
		inline hit_seg trim_hit_seg_copy(hit_seg          arg_hit_seg,  ///< The hit_seg to trim
		                                 const trim_spec &arg_trim_spec ///< The trim_spec to apply
		                                 ) {
			trim_hit_seg( arg_hit_seg, arg_trim_spec );
			return arg_hit_seg;
		}

		/// \brief Get the trimmed start/stop resulting from applying the specified trim_spec on a segment with the specified start/stop
		///
		/// \relates trim_spec
		constexpr residx_residx_pair trim_copy_start_stop(const trim_spec &arg_trim_spec, ///< The trim_spec to apply
		                                                  const residx_t  &arg_start,     ///< The start residue index of the segment in question
		                                                  const residx_t  &arg_stop       ///< The stop  residue index of the segment in question
		                                                  ) {
			return ( arg_stop < arg_start )
				? throw std::invalid_argument("stop must not come before start")
				: residx_residx_pair{
					arg_start + start_trimming_of_length( arg_trim_spec, arg_stop + 1 - arg_start ),
					arg_stop  - stop_trimming_of_length ( arg_trim_spec, arg_stop + 1 - arg_start )
				};
		}

	} // namespace rslv
} // namespace cath

#endif
