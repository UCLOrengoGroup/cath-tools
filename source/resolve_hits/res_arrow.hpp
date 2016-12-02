/// \file
/// \brief The res_arrow class header

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

#ifndef _CATH_TOOLS_SOURCE_RESOLVE_HITS_RES_ARROW_H
#define _CATH_TOOLS_SOURCE_RESOLVE_HITS_RES_ARROW_H

#include <boost/operators.hpp>

#include "exception/invalid_argument_exception.hpp"
#include "resolve_hits/resolve_hits_type_aliases.hpp"

#include <iosfwd>

namespace cath {
	namespace rslv {

		/// \brief Represent the break-point between two neighouring residues
		///
		/// This makes many calculations a bit simpler to reason about than
		/// using raw residue indices.
		///
		/// All conversions between residx_t (unsigned int) and res_arrow,
		/// should be done through the explicit functions:
		///  * `arrow_before_res(const residx_t &)`
		///  * `arrow_after_res (const residx_t &)`
		///  * `res_arrow::res_before()`
		///  * `res_arrow::res_after ()`
		class res_arrow final : private boost::totally_ordered<res_arrow> {
		private:
			friend constexpr res_arrow arrow_before_res(const residx_t &);
			friend constexpr res_arrow arrow_after_res (const residx_t &);

			/// \brief The number used to represent the arrow's position (where 0 means "before the first residue")
			resarw_t arrow;

			explicit constexpr res_arrow(const resarw_t &);

		public:
			constexpr residx_t res_before() const;
			constexpr const residx_t & res_after () const;

			constexpr const resarw_t & get_index() const;

			res_arrow & operator+=(const resarw_t &);
			res_arrow & operator-=(const resarw_t &);
		};

		/// \brief Return whether the first specified arrow appears earlier in the sequence than the second
		///
		/// \relates res_arrow
		constexpr bool operator<(const res_arrow &arg_arrow_a, ///< The first  res_arrow to compare
		                         const res_arrow &arg_arrow_b  ///< The second res_arrow to compare
		                         ) {
			return arg_arrow_a.get_index() < arg_arrow_b.get_index();
		}

		/// \brief Return whether the two specified arrows are identical
		///
		/// \relates res_arrow
		constexpr bool operator==(const res_arrow &arg_arrow_a, ///< The first  res_arrow to compare
		                          const res_arrow &arg_arrow_b  ///< The second res_arrow to compare
		                          ) {
			return arg_arrow_a.get_index() == arg_arrow_b.get_index();
		}

		std::string to_string(const res_arrow &);
		std::ostream & operator<<(std::ostream &,
		                          const res_arrow &);

		constexpr res_arrow arrow_before_res(const residx_t &);
		constexpr res_arrow arrow_after_res (const residx_t &);
		constexpr res_arrow start_arrow();

		/// \brief Ctor from a resarw_t index
		///
		/// This is private and should only be called by
		/// `arrow_before_res(const residx_t &)` and `arrow_after_res (const residx_t &)`
		inline constexpr res_arrow::res_arrow(const resarw_t &arg_res_arrow_index ///< The index of the arrow
		                                      ) : arrow ( arg_res_arrow_index ) {
		}

		/// \brief Get the index of the residue immediately before the arrow
		inline constexpr residx_t res_arrow::res_before() const {
#ifndef NDEBUG
			return ( arrow > 0 ) ? ( arrow - 1 ) : throw( "Cannot return a residue before an arrow that precedes the first residue" );
#else
			return                 ( arrow - 1 );
#endif
		}

		/// \brief Get the index of the residue immediately after the arrow
		inline constexpr const residx_t & res_arrow::res_after () const {
			return arrow;
		}

		/// \brief Getter for the interal index representation
		///
		/// You probably don't want to call this. Consider res_before() and res_after() instead.
		inline constexpr const resarw_t & res_arrow::get_index() const {
			return arrow;
		}

		/// \brief Increment the arrow by the specified offset
		///
		/// \todo Come relaxed constexpr in all supported compilers, make this constexpr
		inline res_arrow & res_arrow::operator+=(const resarw_t &arg_offset ///< The offset by which to increment the res_arrow
		                                         ) {
			arrow += arg_offset;
			return *this;
		}

		/// \brief Decrement the arrow by the specified offset
		///
		/// \todo Come relaxed constexpr in all supported compilers, make this constexpr
		inline res_arrow & res_arrow::operator-=(const resarw_t &arg_offset ///< The offset by which to decrement the res_arrow
		                                         ) {
#ifndef NDEBUG
			if ( arg_offset > arrow ) {
				BOOST_THROW_EXCEPTION(common::invalid_argument_exception("Cannot decrement res_arrow below 0"));
			}
#endif
			arrow -= arg_offset;
			return *this;
		}

		/// \brief Return the result of decrementing the specified arrow by some specified offset
		///
		/// \todo Come constexpr Boost.Operators and relaxed constexpr support in all supported compilers,
		///       implement this using boost::additive<>
		///
		/// \relates res_arrow
		inline constexpr res_arrow operator-(const res_arrow &arg_res_arrow, ///< The res_arrow whose copy should be decremented and returned
		                                     const residx_t  &arg_offset     ///< The offset by which to decrement the the copy of the res_arrow
		                                     ) {


			return
#ifndef NDEBUG
				( arg_res_arrow.res_after() < arg_offset )
				?
					throw std::invalid_argument("Cannot decrement res_arrow beyond 0")
				:
#endif
					arrow_before_res( arg_res_arrow.res_after() - arg_offset );
		}

		/// \brief Return the result of incrementing the specified arrow by some specified offset
		///
		/// \todo Come constexpr Boost.Operators and relaxed constexpr support in all supported compilers,
		///       implement this using boost::additive<>
		///
		/// \relates res_arrow
		inline constexpr res_arrow operator+(const res_arrow &arg_res_arrow, ///< The res_arrow whose copy should be incremented and returned
		                                     const residx_t  &arg_offset     ///< The offset by which to increment the the copy of the res_arrow
		                                     ) {
			return arrow_before_res( arg_res_arrow.res_after() + arg_offset );
		}

		/// \brief Return the result of subtracting one res_arrow from another
		///
		/// \relates res_arrow
		inline constexpr residx_t operator-(const res_arrow &arg_lhs, ///< The res_arrow from which to subtract
		                                    const res_arrow &arg_rhs  ///< The res_arrow to subtract
		                                    ) {
			return arg_lhs.get_index() - arg_rhs.get_index();
		}

		/// \brief Get the arrow immedately before the residue with the specified index
		///
		/// \relates res_arrow
		inline constexpr res_arrow arrow_before_res(const residx_t &arg_res_index ///< The index of the residue immediately after the arrow to return
		                                            ) {
			/// \todo Come C++17, if Herb Sutter has gotten his way (n4029), just use braced list here
			return res_arrow{ arg_res_index     };
		}

		/// \brief Get the arrow immediately after the residue with the specified index
		///
		/// \relates res_arrow
		inline constexpr res_arrow arrow_after_res(const residx_t &arg_res_index ///< The index of the residue immediately before the arrow to return
		                                           ) {
			/// \todo Come C++17, if Herb Sutter has gotten his way (n4029), just use braced list here
			return res_arrow{ arg_res_index + 1 };
		}

		/// \brief Get the start arrow (ie the arrow before the residue with index 0)
		///
		/// \relates res_arrow
		inline constexpr res_arrow start_arrow() {
			return arrow_before_res( 0 );
		}



	} // namespace rslv
} // namespace cath

#endif
