/// \file
/// \brief The equal_group_itr class header

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

#ifndef CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_BOOST_ADDENDA_RANGE_ADAPTOR_ITERATOR_EQUAL_GROUP_ITR_HPP
#define CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_BOOST_ADDENDA_RANGE_ADAPTOR_ITERATOR_EQUAL_GROUP_ITR_HPP

#include <boost/range/sub_range.hpp>

#include "cath/common/boost_addenda/range/range_concept_type_aliases.hpp"
#include "cath/common/exception/invalid_argument_exception.hpp"

// clang-format off
namespace cath::common { template <class RNG> class equal_group_itr; }
// clang-format on

namespace cath::common {
	namespace detail {

		/// \brief TODOCUMENT
		///
		/// For a const_iterator, use a const RNG type.
		template <typename RNG>
		using equal_group_itr_impl = boost::iterator_adaptor<equal_group_itr<RNG>,         // Derived
		                                                     range_iterator_t<RNG>,        // Base
		                                                     boost::sub_range<RNG>,        // Value
		                                                     boost::forward_traversal_tag, // CategoryOrTraversal. ///< \todo There's no fundamental reason this couldn't be bidirectional or random-access; only that implementing the extras would be more work that hasn't yet been warranted.
		                                                     boost::sub_range<RNG>         // Reference
		                                                     >;

	} // namespace detail

	/// \brief TODOCUMENT
	template <class RNG>
	class equal_group_itr final : public detail::equal_group_itr_impl<RNG> {
	private:
		friend class boost::iterator_core_access;

		/// \brief TODOCUMENT
		using super              = detail::equal_group_itr_impl<RNG>;

		/// \brief TODOCUMENT
		using sub_range          = boost::sub_range<RNG>;

		/// \brief TODOCUMENT
		using base_iterator_type = range_iterator_t<RNG>;

		/// \brief TODOCUMENT
		using orig_value_type    = range_value_t<RNG>;

		/// \brief TODOCUMENT
		using ineq_fn_type       = std::function<bool(const orig_value_type &, const orig_value_type &)>;

	public:
		/// \brief TODOCUMENT
		using difference_type    = typename super::difference_type;

	private:
		/// \brief TODOCUMENT
		base_iterator_type equal_group_end_itr;

		/// \brief TODOCUMENT
		base_iterator_type end_itr;

		/// \brief TODOCUMENT
		ineq_fn_type unequal_function;

		base_iterator_type checked_next(const base_iterator_type &,
		                                const base_iterator_type &,
		                                const ineq_fn_type &,
		                                const bool & = false);

		sub_range dereference() const;
		bool equal(const equal_group_itr &) const;
		void increment();

	public:
		equal_group_itr() = default;

		template <typename FN>
		equal_group_itr(const base_iterator_type &,
		                const base_iterator_type &,
		                FN = std::not_equal_to<orig_value_type>() );

		template <typename FN>
		explicit equal_group_itr(RNG &,
		                         FN = std::not_equal_to<orig_value_type>() );

		base_iterator_type get_begin_itr() const;
	};

	/// \brief TODOCUMENT
	///
	/// Note that this isn't very efficient for ranges that aren't random access.
	/// If that turns out to be an issue, then this could be changed.
	template <class RNG>
	typename equal_group_itr<RNG>::base_iterator_type equal_group_itr<RNG>::checked_next(const base_iterator_type &prm_itr,              ///< TODOCUMENT
	                                                                                     const base_iterator_type &prm_end_itr,          ///< TODOCUMENT
	                                                                                     const ineq_fn_type       &prm_unequal_function, ///< TODOCUMENT
	                                                                                     const bool               &prm_permit_end_itr    ///< Whether to permit ( prm_itr == prm_end_itr ) and just return prm_itr unincremented, or just throw.
	                                                                                                                                     ///< The former is needed for when constructing an end equal_group_itr
	                                                                                     ) {
		// Check whether this is an end_itr and if so then either return it unincremented
		// or throw, depending on the value of prm_permit_end_itr
		if ( prm_itr == prm_end_itr ) {
			if ( prm_permit_end_itr ) {
				return prm_itr;
			}
			BOOST_THROW_EXCEPTION(invalid_argument_exception("Attempt to increment end equal_group_itr"));
		}

		const auto next_itr = std::adjacent_find(
			prm_itr,
			prm_end_itr,
			prm_unequal_function
		);
		return ( next_itr == prm_end_itr ) ? prm_end_itr
		                                   : std::next( next_itr );
	}

	/// \brief TODOCUMENT
	template <class RNG>
	typename equal_group_itr<RNG>::sub_range equal_group_itr<RNG>::dereference() const {
		return { this->base(), equal_group_end_itr };
	}

	/// \brief TODOCUMENT
	///
	/// \todo This doesn't check the
	template <class RNG>
	bool equal_group_itr<RNG>::equal(const equal_group_itr &prm_equal_group_itr ///< TODOCUMENT
	                                 ) const {
		return ( get_begin_itr() == prm_equal_group_itr.get_begin_itr() );
	}

	/// \brief TODOCUMENT
	template <class RNG>
	void equal_group_itr<RNG>::increment() {
		this->base_reference() = equal_group_end_itr;
		if ( equal_group_end_itr != end_itr ) {
			equal_group_end_itr = checked_next( equal_group_end_itr, end_itr, unequal_function );
		}
	}

	/// \brief TODOCUMENT
	template <class RNG>
	template <class FN>
	equal_group_itr<RNG>::equal_group_itr(const base_iterator_type &prm_begin_itr,       ///< TODOCUMENT
	                                      const base_iterator_type &prm_end_itr,         ///< TODOCUMENT
	                                      FN                        prm_unequal_function ///< TODOCUMENT
	                                      ) : super              ( prm_begin_itr          ),
	                                          equal_group_end_itr( checked_next(
	                                          	prm_begin_itr,
	                                          	prm_end_itr,
	                                          	prm_unequal_function,
	                                          	true
	                                          ) ),
	                                          end_itr            ( prm_end_itr            ),
	                                          unequal_function   ( prm_unequal_function   ) {
		if ( ! unequal_function ) {
			BOOST_THROW_EXCEPTION(invalid_argument_exception("Unable to construct equal_group_itr with null unequal_function"));
		}
	}

	/// \brief TODOCUMENT
	template <class RNG>
	template <class FN>
	equal_group_itr<RNG>::equal_group_itr(RNG &prm_range,             ///< TODOCUMENT
	                                      FN   prm_less_than_function ///< TODOCUMENT
	                                      ) : equal_group_itr(
	                                          	std::begin( prm_range ),
	                                          	std::end  ( prm_range ),
	                                          	prm_less_than_function
	                                          ){
	}

	/// \brief TODOCUMENT
	template <class RNG>
	typename equal_group_itr<RNG>::base_iterator_type equal_group_itr<RNG>::get_begin_itr() const {
		return this->base();
	}

} // namespace cath::common

#endif // CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_BOOST_ADDENDA_RANGE_ADAPTOR_ITERATOR_EQUAL_GROUP_ITR_HPP
