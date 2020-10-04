/// \file
/// \brief The seq_arrow class header

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

#ifndef _CATH_TOOLS_SOURCE_SEQ_SEQ_ARROW_HPP
#define _CATH_TOOLS_SOURCE_SEQ_SEQ_ARROW_HPP

#include <boost/operators.hpp>

#include "common/exception/invalid_argument_exception.hpp"
#include "seq/seq_type_aliases.hpp"

#include <iosfwd>

namespace cath {
	namespace seq {

		/// \brief Represent the break-point between two neighbouring residues
		///
		/// This makes many calculations a bit simpler to reason about than
		/// using raw residue indices.
		///
		/// All conversions between residx_t (unsigned int) and seq_arrow,
		/// should be done through the explicit functions:
		///  * `arrow_before_res(const residx_t &)`
		///  * `arrow_after_res (const residx_t &)`
		///  * `seq_arrow::res_before()`
		///  * `seq_arrow::res_after ()`
		class seq_arrow final : private boost::totally_ordered<seq_arrow> {
		private:
			friend constexpr seq_arrow arrow_before_res(const residx_t &);
			friend constexpr seq_arrow arrow_after_res (const residx_t &);

			/// \brief The number used to represent the arrow's position (where 0 means "before the first residue")
			resarw_t arrow;

			explicit constexpr seq_arrow(const resarw_t &);

		public:
			constexpr residx_t res_before() const;
			constexpr const residx_t & res_after() const;

			constexpr const resarw_t & get_index() const;

			seq_arrow & operator+=(const resarw_t &);
			seq_arrow & operator-=(const resarw_t &);
		};

		/// \brief Return whether the first specified arrow appears earlier in the sequence than the second
		///
		/// \relates seq_arrow
		constexpr bool operator<(const seq_arrow &prm_arrow_a, ///< The first  seq_arrow to compare
		                         const seq_arrow &prm_arrow_b  ///< The second seq_arrow to compare
		                         ) {
			return prm_arrow_a.get_index() < prm_arrow_b.get_index();
		}

		/// \brief Return whether the two specified arrows are identical
		///
		/// \relates seq_arrow
		constexpr bool operator==(const seq_arrow &prm_arrow_a, ///< The first  seq_arrow to compare
		                          const seq_arrow &prm_arrow_b  ///< The second seq_arrow to compare
		                          ) {
			return prm_arrow_a.get_index() == prm_arrow_b.get_index();
		}

		std::string to_string(const seq_arrow &);
		std::ostream & operator<<(std::ostream &,
		                          const seq_arrow &);

		constexpr seq_arrow arrow_before_res(const residx_t &);
		constexpr seq_arrow arrow_after_res (const residx_t &);
		constexpr seq_arrow start_arrow();

		/// \brief Ctor from a resarw_t index
		///
		/// This is private and should only be called by
		/// `arrow_before_res(const residx_t &)` and `arrow_after_res (const residx_t &)`
		inline constexpr seq_arrow::seq_arrow(const resarw_t &prm_res_arrow_index ///< The index of the arrow
		                                      ) : arrow ( prm_res_arrow_index ) {
		}

		/// \brief Get the index of the residue immediately before the arrow
		inline constexpr residx_t seq_arrow::res_before() const {
#ifndef NDEBUG
			return ( arrow > 0 ) ? ( arrow - 1 ) : throw( "Cannot return a residue before an arrow that precedes the first residue" );
#else
			return                 ( arrow - 1 );
#endif
		}

		/// \brief Get the index of the residue immediately after the arrow
		inline constexpr const residx_t & seq_arrow::res_after () const {
			return arrow;
		}

		/// \brief Getter for the internal index representation
		///
		/// You probably don't want to call this. Consider res_before() and res_after() instead.
		inline constexpr const resarw_t & seq_arrow::get_index() const {
			return arrow;
		}

		/// \brief Increment the arrow by the specified offset
		///
		/// \todo Come relaxed constexpr in all supported compilers, make this constexpr
		inline seq_arrow & seq_arrow::operator+=(const resarw_t &prm_offset ///< The offset by which to increment the seq_arrow
		                                         ) {
			arrow += prm_offset;
			return *this;
		}

		/// \brief Decrement the arrow by the specified offset
		///
		/// \todo Come relaxed constexpr in all supported compilers, make this constexpr
		inline seq_arrow & seq_arrow::operator-=(const resarw_t &prm_offset ///< The offset by which to decrement the seq_arrow
		                                         ) {
#ifndef NDEBUG
			if ( prm_offset > arrow ) {
				BOOST_THROW_EXCEPTION(common::invalid_argument_exception("Cannot decrement seq_arrow below 0"));
			}
#endif
			arrow -= prm_offset;
			return *this;
		}

		/// \brief Return the result of decrementing the specified arrow by some specified offset
		///
		/// \todo Come constexpr Boost.Operators and relaxed constexpr support in all supported compilers,
		///       implement this using boost::additive<>
		///
		/// \relates seq_arrow
		inline constexpr seq_arrow operator-(const seq_arrow &prm_res_arrow, ///< The seq_arrow whose copy should be decremented and returned
		                                     const residx_t  &prm_offset     ///< The offset by which to decrement the copy of the seq_arrow
		                                     ) {


			return
#ifndef NDEBUG
				( prm_res_arrow.res_after() < prm_offset )
				?
					throw std::invalid_argument("Cannot decrement seq_arrow beyond 0")
				:
#endif
					arrow_before_res( prm_res_arrow.res_after() - prm_offset );
		}

		/// \brief Return the result of incrementing the specified arrow by some specified offset
		///
		/// \todo Come constexpr Boost.Operators and relaxed constexpr support in all supported compilers,
		///       implement this using boost::additive<>
		///
		/// \relates seq_arrow
		inline constexpr seq_arrow operator+(const seq_arrow &prm_res_arrow, ///< The seq_arrow whose copy should be incremented and returned
		                                     const residx_t  &prm_offset     ///< The offset by which to increment the copy of the seq_arrow
		                                     ) {
			return arrow_before_res( prm_res_arrow.res_after() + prm_offset );
		}

		/// \brief Return the result of subtracting one seq_arrow from another
		///
		/// \relates seq_arrow
		inline constexpr residx_t operator-(const seq_arrow &prm_lhs, ///< The seq_arrow from which to subtract
		                                    const seq_arrow &prm_rhs  ///< The seq_arrow to subtract
		                                    ) {
			return prm_lhs.get_index() - prm_rhs.get_index();
		}

		/// \brief Get the arrow immediately before the residue with the specified index
		///
		/// \relates seq_arrow
		inline constexpr seq_arrow arrow_before_res(const residx_t &prm_res_index ///< The index of the residue immediately after the arrow to return
		                                            ) {
			/// \todo Come C++17, if Herb Sutter has gotten his way (n4029), just use braced list here
			return seq_arrow{ prm_res_index     };
		}

		/// \brief Get the arrow immediately after the residue with the specified index
		///
		/// \relates seq_arrow
		inline constexpr seq_arrow arrow_after_res(const residx_t &prm_res_index ///< The index of the residue immediately before the arrow to return
		                                           ) {
			/// \todo Come C++17, if Herb Sutter has gotten his way (n4029), just use braced list here
			return seq_arrow{ prm_res_index + 1 };
		}

		/// \brief Get the start arrow (ie the arrow before the residue with index 0)
		///
		/// \relates seq_arrow
		inline constexpr seq_arrow start_arrow() {
			return arrow_before_res( 0 );
		}



	} // namespace seq
} // namespace cath

#endif
