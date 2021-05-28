/// \file
/// \brief The cross_itr class header

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

#ifndef _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_BOOST_ADDENDA_RANGE_UTILITY_ITERATOR_CROSS_ITR_HPP
#define _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_BOOST_ADDENDA_RANGE_UTILITY_ITERATOR_CROSS_ITR_HPP

#include <boost/iterator/iterator_facade.hpp>

#include "cath/common/boost_addenda/iterator/iterator_traits_type_aliases.hpp"
#include "cath/common/boost_addenda/range/range_concept_type_aliases.hpp"
#include "cath/common/config.hpp"
#include "cath/common/cpp17/apply.hpp"
#include "cath/common/exception/invalid_argument_exception.hpp"

// clang-format off
namespace cath::common { template <typename... RNGs> class cross_itr; }
// clang-format on

namespace cath::common {
	namespace detail {

		/// \brief Convenience type alias for the iterator_facade through which cross_itr is implemented
		///
		/// \todo This currently requires all inputs to meet forward t all inputs to provide forward traversal.
		///       Conversely, it doesn't provide bidirectional/random-access where the inputs would allow.
		///       Consider finding a way to set the traversal to the weakest traversal from the inputs
		///       and add bidirectional/random-access functionality where it'd be useful
		template <typename... RNGs>
		using cross_itr_impl = boost::iterator_facade<cross_itr<RNGs...>,                    // Derived
		                                              std::tuple<range_value_t<RNGs>...>,    // Value
		                                              boost::forward_traversal_tag,          // CategoryOrTraversal
		                                              std::tuple<range_reference_t<RNGs>...> // Reference
		                                              >;

		/// \brief TODOCUMENT
		struct dereferencer final {
			template <typename... Is>
			constexpr auto operator()(const Is &... prm_iters
			                          ) {
				return std::tuple<iterator_reference_t<Is>... >( ( *prm_iters )... );
			}
		};

		/// \brief TODOCUMENT
		///
		/// \returns Whether this iterator was incremented rather than wrapped
		template <typename I>
		inline bool increment_or_wrap(I       &prm_itr,       ///< TODOCUMENT
		                              const I &prm_begin_itr, ///< TODOCUMENT
		                              const I &prm_end_itr    ///< TODOCUMENT
		                              ) {
			if constexpr ( IS_IN_DEBUG_MODE ) {
				if ( prm_begin_itr == prm_end_itr ) {
					BOOST_THROW_EXCEPTION( invalid_argument_exception(
					  "Unable to do performed cross iteration over ranges including an empty range" ) );
				}
			}
			++prm_itr;
			const bool wrapped = ( prm_itr == prm_end_itr );
			if ( wrapped ) {
				prm_itr = prm_begin_itr;
			}
			return ! wrapped;
		}

		/// \brief TODOCUMENT
		template <size_t N, typename... Is>
		struct increment_impl final {
			static inline void increment(std::tuple<Is...> &prm_iters,       ///< TODOCUMENT
			                             std::tuple<Is...> &prm_begin_iters, ///< TODOCUMENT
			                             std::tuple<Is...> &prm_end_iters,   ///< TODOCUMENT
			                             const bool        &prm_complete
			                             ) {
				const bool is_complete = prm_complete ? true
				                                      : increment_or_wrap(
				                                        	std::get<N>( prm_iters       ),
				                                        	std::get<N>( prm_begin_iters ),
				                                        	std::get<N>( prm_end_iters   )
				                                        );
				increment_impl< N - 1, Is...>::increment( prm_iters, prm_begin_iters, prm_end_iters, is_complete );
			}
		};

		/// \brief TODOCUMENT
		template <typename... Is>
		struct increment_impl<0, Is...> final {
			static inline void increment(std::tuple<Is...> &prm_iters,       ///< TODOCUMENT
			                             std::tuple<Is...> &prm_begin_iters, ///< TODOCUMENT
			                             std::tuple<Is...> &prm_end_iters,   ///< TODOCUMENT
			                             const bool        &prm_complete     ///< TODOCUMENT
			                             ) {
				if ( ! prm_complete ) {
					const bool complete = increment_or_wrap(
						std::get<0>( prm_iters       ),
						std::get<0>( prm_begin_iters ),
						std::get<0>( prm_end_iters   )
					);
					if ( ! complete ) {
						prm_iters = prm_end_iters;
					}
				}
			}
		};

		/// \brief TODOCUMENT
		template <typename... Is>
		void increment(std::tuple<Is...> &prm_iters,       ///< TODOCUMENT
		               std::tuple<Is...> &prm_begin_iters, ///< TODOCUMENT
		               std::tuple<Is...> &prm_end_iters    ///< TODOCUMENT
		               ) {
			if ( sizeof...( Is ) > 0 ) {
				increment_impl< sizeof...( Is ) - 1, Is... >::increment( prm_iters, prm_begin_iters, prm_end_iters, false );
			}
		}

	} // namespace detail

	/// \brief TODOCUMENT
	///
	/// Invariants:
	///  * The client must preserve the validity of the original ranges and their one-past-end iterators throughout the cross_itr's lifetime
	template <typename... RNGs>
	class cross_itr final : public detail::cross_itr_impl<RNGs...> {
	private:
		friend class boost::iterator_core_access;

		/// \brief TODOCUMENT
		using super               = detail::cross_itr_impl<RNGs...>;

		/// \brief TODOCUMENT
		using reference_type      = std::tuple<range_reference_t<RNGs>...>;

		/// \brief TODOCUMENT
		using iterator_type       = std::tuple<range_iterator_t<RNGs>...>;

	private:
		/// \brief TODOCUMENT
		iterator_type the_iterators;

		/// \brief TODOCUMENT
		iterator_type the_begin_iterators;

		/// \brief TODOCUMENT
		iterator_type the_end_iterators;

		reference_type dereference() const;
		void increment();

		template <typename... OTHER_RNGs>
		bool equal(const cross_itr<OTHER_RNGs...> &) const;

	public:
		cross_itr();
		cross_itr(const iterator_type &,
		          const iterator_type &);
		explicit cross_itr(RNGs &...);
	};

	/// \brief TODOCUMENT
	template <typename... RNGs>
	auto cross_itr<RNGs...>::dereference() const -> reference_type {
		/// \TODO Come C++17, use ::std::apply
		return ::cath::common::apply( detail::dereferencer(), the_iterators );
	}

	/// \brief TODOCUMENT
	template <typename... RNGs>
	void cross_itr<RNGs...>::increment() {
		detail::increment( the_iterators, the_begin_iterators, the_end_iterators );
	}

	/// \brief TODOCUMENT
	template <typename... RNGs>
	template <typename... OTHER_RNGs>
	bool cross_itr<RNGs...>::equal(const cross_itr<OTHER_RNGs...> &prm_cross_itr ///< TODOCUMENT
	                               ) const {
		return ( the_iterators == prm_cross_itr.the_iterators );
	}

	/// \brief Ctor of singular iterator
	/// \brief TODOCUMENT
	template <typename... RNGs>
	cross_itr<RNGs...>::cross_itr() = default;

	/// \brief Ctor from iterators
	template <typename... RNGs>
	cross_itr<RNGs...>::cross_itr(const iterator_type &prm_begin, ///< Begin iterator for the original range
	                              const iterator_type &prm_end    ///< End   iterator for the original range
	                              ) : the_iterators      ( prm_begin ),
	                                  the_begin_iterators( prm_begin ),
	                                  the_end_iterators  ( prm_end   ) {
	}

	/// \brief Ctor from ranges
	template <typename... RNGs>
	cross_itr<RNGs...>::cross_itr(RNGs &... prm_ranges ///< The ranges over which this cross_itr should act
	                              ) : cross_itr(
	                                  	make_tuple( std::begin( prm_ranges )... ),
	                                  	make_tuple( std::end  ( prm_ranges )... )
	                                  ) {
	}

	/// \brief TODOCUMENT
	template <typename... RNGs>
	cross_itr<RNGs...> make_cross_itr(RNGs &... prm_ranges ///< The ranges over which this cross_itr should act
	                                  ) {
		return cross_itr<RNGs...>(
			std::make_tuple( std::begin( prm_ranges )... ),
			std::make_tuple( std::end  ( prm_ranges )... )
		);
	}

	/// \brief TODOCUMENT
	template <typename... RNGs>
	cross_itr<RNGs...> make_end_cross_itr(RNGs &... prm_ranges ///< The ranges over which this cross_itr should act
	                                      ) {
		return cross_itr<RNGs...>(
			std::make_tuple( std::end( prm_ranges )... ),
			std::make_tuple( std::end( prm_ranges )... )
		);
	}

	namespace detail {
		template <typename... T> class TD;
		/// \brief TODOCUMENT
		struct cross_itr_maker final {
			template <typename... RNGs>
			auto operator()(RNGs &... prm_ranges
			                ) {
				return make_cross_itr<std::remove_reference_t<RNGs>...>( prm_ranges... );
			}
		};

		/// \brief TODOCUMENT
		struct end_cross_itr_maker final {
			template <typename... RNGs>
			auto operator()(RNGs &... prm_ranges
			                ) {
				return make_end_cross_itr<std::remove_reference_t<RNGs>...>( prm_ranges... );
			}
		};

		template <typename... T> class TD;
		/// \brief TODOCUMENT
		struct const_cross_itr_maker final {
			template <typename... RNGs>
			auto operator()(RNGs &... prm_ranges
			                ) {
				return make_cross_itr<const std::remove_reference_t<RNGs>...>( prm_ranges... );
			}
		};

		/// \brief TODOCUMENT
		struct end_const_cross_itr_maker final {
			template <typename... RNGs>
			auto operator()(RNGs &... prm_ranges
			                ) {
				return make_end_cross_itr<const std::remove_reference_t<RNGs>...>( prm_ranges... );
			}
		};

	} // namespace detail

	/// \brief Ctor from a tuple of ranges
	template <typename... RNGs>
	cross_itr<std::remove_reference_t<RNGs>...> make_cross_itr(std::tuple<RNGs ...> &prm_tuple ///< The ranges over which this cross_itr should act
	                                                           ) {

		/// \TODO Come C++17, use ::std::apply
		return ::cath::common::apply( detail::cross_itr_maker(), prm_tuple );
	}

	/// \brief Ctor from a tuple of ranges
	template <typename... RNGs>
	cross_itr<std::remove_reference_t<RNGs>...> make_end_cross_itr(std::tuple<RNGs ...> &prm_tuple ///< The ranges over which this cross_itr should act
	                                                               ) {
		/// \TODO Come C++17, use ::std::apply
		return ::cath::common::apply( detail::end_cross_itr_maker(), prm_tuple );
	}

	/// \brief Ctor from a tuple of ranges
	///
	/// \todo See notes for cross(const tuple &)
	template <typename... RNGs>
	cross_itr<const std::remove_reference_t<RNGs>...> make_cross_itr(const std::tuple<RNGs ...> &prm_tuple ///< The ranges over which this cross_itr should act
	                                                                 ) {
//		// If RNGs contains non-const references, they shouldn't all be made into const references
//		// just because calling get<>( ) on the const tuple returns a const reference.
//		//
//		// So make a local, non-const tuple of references and build from that instead
//		std::tuple<std::add_lvalue_reference_t<RNGs>...> non_const_ref_copy( prm_tuple );

		/// \TODO Come C++17, use ::std::apply
		return ::cath::common::apply( detail::const_cross_itr_maker(), prm_tuple );
	}

	/// \brief Ctor from a tuple of ranges
	///
	/// \todo See notes for cross(const tuple &)
	template <typename... RNGs>
	cross_itr<const std::remove_reference_t<RNGs>...> make_end_cross_itr(const std::tuple<RNGs ...> &prm_tuple ///< The ranges over which this cross_itr should act
	                                                                     ) {
//		// If RNGs contains non-const references, they shouldn't all be made into const references
//		// just because calling get<>( ) on the const tuple returns a const reference.
//		//
//		// So make a local, non-const tuple of references and build from that instead
//		std::tuple<std::add_lvalue_reference_t<RNGs>...> non_const_ref_copy( prm_tuple );

		/// \TODO Come C++17, use ::std::apply
		return ::cath::common::apply( detail::end_const_cross_itr_maker(), prm_tuple );
	}

	template <typename... Ts>
	using cross_range = boost::iterator_range<cross_itr<Ts...>>;

	/// \brief TODOCUMENT
	template <typename... RNGs>
	cross_range<RNGs...> cross(RNGs &... prm_ranges ///< The ranges over which this cross_itr should act
	                           ) {
		return {
			make_cross_itr    ( prm_ranges... ),
			make_end_cross_itr( prm_ranges... ),
		};
	}

	/// \brief Ctor from a tuple of ranges
	template <typename... RNGs>
	cross_range<std::remove_reference_t<RNGs>...> cross(std::tuple<RNGs ...> &prm_tuple ///< The ranges over which this cross_itr should act
	                                                    ) {
		return {
			make_cross_itr    ( prm_tuple ),
			make_end_cross_itr( prm_tuple ),
		};
	}

	/// \brief Ctor from a tuple of ranges
	///
	/// \todo This should be able to take non-const reference types for some (or all) of the RNGs and create a cross_range
	///       that allows non-const iteration over those types' ranges. That won't happen if prm_tuple is used directly
	///       because the get<>( ) call on the const tuple will return a *const* reference. This can be solved by taking a
	///       non-const copy of prm_tuple and using that instead. However it's important that this copy not involve
	///       taking local copies of any non-reference types in RNGs (because that would create a cross_range over
	///       a soon-to-die temporary). It would be a really nice solution to copy prm_tuple to a non-const tuple
	///       of references to ranges:
	///       ~~~~~.cpp
	///       std::tuple<std::add_lvalue_reference_t<RNGs>...> non_const_ref_copy( prm_tuple );
	///       ~~~~~
	///       which would avoid adding spurious consts to non-const RNGs and would avoid copying non-reference RNGs.
	///       Unfortunately, this doesn't work because the tuple constructor won't allow it. It appears that there
	///       may be no good reason for this and it may be a flaw in the C++11/C++17 tuple ctor.
	///       (see http://stackoverflow.com/a/24919135)
	///       Come C++17, try again.
	///
	/// For now, if you want to use this for anything other than non-const RNGs, use the version above that
	/// takes a non-const reference to a tuple
	template <typename... RNGs>
	cross_range<const std::remove_reference_t<RNGs>...> cross(const std::tuple<RNGs ...> &prm_tuple ///< The ranges over which this cross_itr should act
	                                                          ) {
//		// If RNGs contains non-const references, they shouldn't all be made into const references
//		// just because calling get<>( ) on the const tuple returns a const reference.
//		//
//		// So make a local, non-const tuple of references and build from that instead
//		std::tuple<std::add_lvalue_reference_t<RNGs>...> non_const_ref_copy( prm_tuple );
		return {
			make_cross_itr    ( prm_tuple ),
			make_end_cross_itr( prm_tuple ),
		};
	}

} // namespace cath::common

#endif // _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_BOOST_ADDENDA_RANGE_UTILITY_ITERATOR_CROSS_ITR_HPP
