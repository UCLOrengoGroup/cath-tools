/// \file
/// \brief The hit_extras class header

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

#ifndef _CATH_TOOLS_SOURCE_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_HIT_EXTRAS_HPP
#define _CATH_TOOLS_SOURCE_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_HIT_EXTRAS_HPP

#include <optional>
#include <type_traits>
#include <variant>
#include <vector>

#include <boost/algorithm/string/join.hpp>
#include <boost/range/adaptor/transformed.hpp>

#include "cath/common/cpp20/constexpr_invoke.hpp"
#include "cath/common/exception/invalid_argument_exception.hpp"

namespace cath::rslv {

	/// \brief The category of the item of extra information for a hit
	enum class hit_extra_cat : char {
		ALND_RGNS, ///< Aligned regions as parsed from, eg, hmmsearch output
		COND_EVAL, ///< Conditional evalue from HMMER
		INDP_EVAL  ///< Independent evalue from HMMER
	};

	std::string to_string(const hit_extra_cat &);

	namespace detail {

		/// \brief Implement a lookup from hit_extra_cat value to the type of the associated information
		template <hit_extra_cat> struct type_of_hit_extra_cat;
		
		/// \brief Implement a lookup from hit_extra_cat value to the type of the associated information
		template <> struct type_of_hit_extra_cat<hit_extra_cat::ALND_RGNS> { using type = std::string; };
		
		/// \brief Implement a lookup from hit_extra_cat value to the type of the associated information
		template <> struct type_of_hit_extra_cat<hit_extra_cat::COND_EVAL> { using type = double;      };
		
		/// \brief Implement a lookup from hit_extra_cat value to the type of the associated information
		template <> struct type_of_hit_extra_cat<hit_extra_cat::INDP_EVAL> { using type = double;      };

	} // namespace detail

	/// \brief A lookup from hit_extra_cat value to the type of the associated information
	template <hit_extra_cat Cat>
	using type_of_hit_extra_cat_t    = typename detail::type_of_hit_extra_cat<Cat>::type;

	/// \brief Variant for storing any of the types of information associated with hit_extra_cat
	///
	/// Exercise: this could be auto generated from an ALL_HIT_EXTRA_CATS std::array using type_of_hit_extra_cat_t and some unique-ing
	using hit_extra_variant          = ::std::variant<std::string, double>;

	/// \brief Type alias for a pair of hit_extra_cat and a hit_extra_variant
	using hit_extra_cat_var_pair     = std::pair<hit_extra_cat, hit_extra_variant>;

	/// \brief Type alias for a vector of hit_extra_cat_var_pair values
	using hit_extra_cat_var_pair_vec = std::vector<hit_extra_cat_var_pair>;


	/// \brief Invoke the specified callable on the type 
	///
	/// \tparam Fn must support being called with a single argument of any of the types if pieces of hit_extra information
	///
	/// \pre The hit_extra_cat and associated hit_extra_variant must match
	//       (ie type_of_hit_extra_cat_t of the hit_extra_cat must match the type stored in the hit_extra_variant)
	template <typename Fn>
	constexpr void invoke_for_hit_extra_info(Fn                           &&prm_fn,  ///< The function to invoke
	                                         const hit_extra_cat_var_pair  &prm_pair ///< The pair of hit_extra_cat and associated hit_extra_variant
	                                         ) {
		switch ( prm_pair.first ) {
			case( hit_extra_cat::ALND_RGNS ) : { common::constexpr_invoke( std::forward<Fn>( prm_fn ), ::std::get< type_of_hit_extra_cat_t< hit_extra_cat::ALND_RGNS > >( prm_pair.second ) ); return; }
			case( hit_extra_cat::COND_EVAL ) : { common::constexpr_invoke( std::forward<Fn>( prm_fn ), ::std::get< type_of_hit_extra_cat_t< hit_extra_cat::COND_EVAL > >( prm_pair.second ) ); return; }
			case( hit_extra_cat::INDP_EVAL ) : { common::constexpr_invoke( std::forward<Fn>( prm_fn ), ::std::get< type_of_hit_extra_cat_t< hit_extra_cat::INDP_EVAL > >( prm_pair.second ) ); return; }
		}
		BOOST_THROW_EXCEPTION(common::invalid_argument_exception("Value of hit_extra_cat not recognised whilst trying to do something with it"));
	}

	namespace detail {

		/// \brief Helper class to stringify the value passed to the function operator into the string
		///        passed to the ctor
		///
		/// This is for use with invoke_for_hit_extra_info
		struct info_stringifier final {

			/// \brief A reference to the string to be populated
			std::reference_wrapper<std::string> result_string;

			/// \brief Ctor to pass the string to be populated into result_string 
			explicit info_stringifier(std::string &prm_string ///< The string to be populated
			                          ) : result_string{ prm_string } {
			}

			/// \brief Stringify the specified value to result_string
			template <typename T>
			void operator()(const T &prm_value ///< The value to stringify
			              ) {

				result_string.get() = std::to_string( prm_value );
			}

		};

		/// \brief Specialisation for stringifying a string
		template <>
		inline void info_stringifier::operator()<std::string>(const std::string &prm_value ///< The string to stringify
		                                                      ) {
			result_string.get() = prm_value;
		}

	} // namespace detail

	/// \brief Return a string of the value information in the specified hit_extra_cat_var_pair
	inline std::string string_of_info(const hit_extra_cat_var_pair &prm_pair ///< The hit_extra_cat_var_pair to describe
	                                  ) {
		std::string result;
		detail::info_stringifier the_stringifier{ result };
		invoke_for_hit_extra_info( the_stringifier, prm_pair );
		return result;
	}





	/// \brief A store of hit extra information, where each piece is stored with the hit_extra_cat and the associated type in a hit_extra_variant
	struct hit_extras_store final {
	private:
		/// \brief The actual store of hit_extra_cat_var_pairs
		hit_extra_cat_var_pair_vec extras;

	public:
		/// \brief Type alias as part of making this a range
		using const_iterator = hit_extra_cat_var_pair_vec::const_iterator;

		/// \brief Add the specified piece of information
		template <hit_extra_cat Cat>
		hit_extras_store & push_back(type_of_hit_extra_cat_t<Cat> prm_extra ///< The piece of information to store
		                             ) {
			extras.emplace_back( Cat, std::move( prm_extra ) );
			return *this;
		}

		/// \brief Standard const begin() operator to provide range access
		[[nodiscard]] const_iterator begin() const {
			return ::std::cbegin( extras );
		}

		/// \brief Standard const end() operator to provide range access
		[[nodiscard]] const_iterator end() const {
			return ::std::cend( extras );
		}

		/// \brief Return whether this is empty
		[[nodiscard]] bool empty() const {
			return extras.empty();
		}

		/// \brief Return the number of pieces of information currently being stored
		[[nodiscard]] size_t size() const {
			return extras.size();
		}
	};

	/// \brief Get the first value in the specified hit_extras_store that matches the specified hit_extra_cat
	///        or return nullopt if there are no such matches
	///
	/// \relates hit_extras_store
	template <hit_extra_cat Cat>
	::std::optional<type_of_hit_extra_cat_t<Cat>> get_first(const hit_extras_store &prm_store ///< The hit_extras_store to query
	                                                        ) {
		for (const auto &extra_pair : prm_store) {
			if ( extra_pair.first == Cat ) {
				return ::std::make_optional( ::std::get<type_of_hit_extra_cat_t<Cat>>( extra_pair.second ) );
			}
		}
		return ::std::nullopt;
	}

	/// \brief Generate a string describing the specified hit_extras_store
	///
	/// \relates hit_extras_store
	inline std::string to_string(const hit_extras_store &prm_store ///< The hit_extras_store to describe
	                             ) {
		return "hit_extras_store["
			+ boost::algorithm::join(
				prm_store
					| boost::adaptors::transformed( [] (const hit_extra_cat_var_pair &x) {
						return to_string( x.first ) + ":" + string_of_info( x );
					} ),
				","
			)
			+ "]";
	}

	/// \brief Insert a description of the specified hit_extras_store into the specified ostream
	///
	/// \relates hit_extras_store
	inline std::ostream & operator<<(std::ostream           &prm_os,   ///< The ostream into which the description should be inserted
	                                 const hit_extras_store &prm_store ///< The hit_extras_store to describe
	                                 ) {
		prm_os << to_string( prm_store );
		return prm_os;
	}

} // namespace cath::rslv

#endif // _CATH_TOOLS_SOURCE_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_HIT_EXTRAS_HPP
