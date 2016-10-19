/// \file
/// \brief The limit_itr class header

/// \copyright
/// Tony Lewis's Common C++ Library Code (here imported into the CATH Tools project and then tweaked, eg namespaced in cath)
/// Copyright (C) 2007, Tony Lewis
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

#ifndef _CATH_TOOLS_SOURCE_COMMON_BOOST_ADDENDA_RANGE_ADAPTOR_ITERATOR_LIMIT_ITR_H
#define _CATH_TOOLS_SOURCE_COMMON_BOOST_ADDENDA_RANGE_ADAPTOR_ITERATOR_LIMIT_ITR_H

#include <boost/iterator/iterator_adaptor.hpp>
#include <boost/numeric/conversion/cast.hpp>

#include "common/boost_addenda/range/range_concept_type_aliases.h"

//#include <ios>

namespace cath { namespace common { template <class RNG> class limit_itr; } }

namespace cath {
	namespace common {
		namespace detail {
//			/// \brief Convenience type alias for a pair of references to the elements in a range of type RNG
//			template <typename RNG>
//			using range_limit_ref_pair = std::pair<range_reference_t<RNG> &,
//			                                       range_reference_t<RNG> &>;

			/// \brief Convenience type alias for the iterator_adaptor through which limit_itr is implemented
			template <typename RNG>
			using limit_itr_impl = boost::iterator_adaptor<limit_itr<RNG>,        // Derived
			                                               range_iterator_t<RNG>, // Base
			                                               range_value_t<RNG>,    // Value
			                                               range_category_t<RNG>, // CategoryOrTraversal
			                                               range_reference_t<RNG> // Reference
			                                               >;
		} // namespace detail

		/// \brief Iterator for wrapping a range and limiting it to, at most, n elements
		///        (but without requiring random_access as sliced does)
		///
		/// This iterator stores the underlying end iterator and uses this to stop safely
		/// if the underlying range has fewer than n elements.
		///
		/// This make it useful for, eg, outputting the first few elements of a range of unknown length
		/// via the limited adaptor:
		///
		/// ~~~~~.cpp
		/// cerr << "First few elements are : " << join( my_numbers | limited( 5 ) | lexical_cast<string>(), ", " ) << endl;
		/// ~~~~~
		///
		/// The underlying range shouldn't be modified during the iterator's lifetime.
		///
		/// For a const_iterator, use a const RNG type.
		///
		/// Invariants:
		///  * The base() iterator should always point to some member of the original range
		///  * The client must preserve the validity of the original range and its one-past-end iterator throughout the limit_itr's lifetime
		template <typename RNG>
		class limit_itr final : public detail::limit_itr_impl<RNG> {
		private:
			friend class boost::iterator_core_access;

			/// \brief TODOCUMENT
			using super              = detail::limit_itr_impl<RNG>;

			/// \brief TODOCUMENT
			using base_iterator_type = range_iterator_t<RNG>;

//			/// \brief TODOCUMENT
//			using limit_pair_type    = detail::range_limit_ref_pair<RNG>;

		public:
			/// \brief TODOCUMENT
			using difference_type    = typename super::difference_type;

//			/// \brief TODOCUMENT
//			using reference_type     = range_reference_t<RNG>;

		private:
			/// \brief TODOCUMENT
			difference_type offset = 0;

			/// \brief TODOCUMENT
			base_iterator_type end_itr;

			/// \brief TODOCUMENT
			difference_type max_num_elements;

			void advance(const difference_type &);
			void decrement();
			void increment();

			bool is_at_end() const;
			difference_type distance_to_end() const;

			template <typename OTHER_RNG>
			bool equal(const limit_itr<OTHER_RNG> &) const;

			template <class OTHER_RNG>
			difference_type distance_to(const limit_itr<OTHER_RNG> &) const;

		public:
			limit_itr(const base_iterator_type &,
			          const base_iterator_type &,
			          const size_t &);
			limit_itr(RNG &,
			          const size_t &);
		};

		/// \brief TODOCUMENT
		template <class RNG>
		void limit_itr<RNG>::advance(const difference_type &arg_offset ///< TODOCUMENT
		                             ) {
			( this->base_reference() ) += arg_offset;
			offset += arg_offset;
		}

		/// \brief TODOCUMENT
		template <class RNG>
		void limit_itr<RNG>::decrement() {
			--( this->base_reference() );
			--offset;
		}

		/// \brief TODOCUMENT
		template <class RNG>
		void limit_itr<RNG>::increment() {
			++( this->base_reference() );
			++offset;
		}

		/// \brief TODOCUMENT
		template <class RNG>
		bool limit_itr<RNG>::is_at_end() const {
			return
				( offset >= max_num_elements )
				||
				( this->base() == end_itr );
		}

		/// \brief TODOCUMENT
		template <class RNG>
		typename limit_itr<RNG>::difference_type limit_itr<RNG>::distance_to_end() const {
			return std::min(
				end_itr          - this->base(),
				max_num_elements - offset
			);
		}

		/// \brief TODOCUMENT
		template <class RNG>
		template <typename OTHER_RNG>
		bool limit_itr<RNG>::equal(const limit_itr<OTHER_RNG> &arg_limit_itr ///< TODOCUMENT
		                           ) const {
			if ( end_itr == arg_limit_itr.end_itr && is_at_end() && arg_limit_itr.is_at_end() ) {
				return true;
			}
			if ( ! is_at_end() && ! arg_limit_itr.is_at_end() ) {
				if ( this->base() == arg_limit_itr.base() ) {
					return true;
				}
			}
			return false;
		}

		/// \brief TODOCUMENT
		template <class RNG>
		template <class OTHER_RNG>
		typename limit_itr<RNG>::difference_type limit_itr<RNG>::distance_to(const limit_itr<OTHER_RNG> &arg_other_itr ///< TODOCUMENT
		                                                                     ) const {
			if ( is_at_end() || arg_other_itr.is_at_end() ) {
				if ( is_at_end() && arg_other_itr.is_at_end() ) {
					return 0;
				}
				return arg_other_itr.is_at_end() ? distance_to_end()
				                                 : - arg_other_itr.distance_to_end();
			}
			return ( arg_other_itr.base() - this->base() );
		}


		/// \brief Ctor from iterators
		template <class RNG>
		limit_itr<RNG>::limit_itr(const base_iterator_type &arg_begin,           ///< Begin iterator for the original range
		                          const base_iterator_type &arg_end,             ///< End   iterator for the original range
		                          const size_t             &arg_max_num_elements ///< TODOCUMENT
		                          ) : super           ( arg_begin           ),
		                              end_itr         ( arg_end             ),
		                              max_num_elements( boost::numeric_cast<difference_type>( arg_max_num_elements ) ) {
		}

		/// \brief Ctor from a range
		template <class RNG>
		limit_itr<RNG>::limit_itr(RNG          &arg_range,           ///< The range over which this limit_itr should act
		                          const size_t &arg_max_num_elements ///< TODOCUMENT
		                          ) : limit_itr(
		                              	std::begin( arg_range ),
		                              	std::end  ( arg_range ),
		                              	arg_max_num_elements
		                              ) {
		}

	} // namespace common
} // namespace cath

#endif
