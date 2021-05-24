/// \file
/// \brief The prog_opt_num_range header

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

#ifndef _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_PROGRAM_OPTIONS_PROG_OPT_NUM_RANGE_HPP
#define _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_PROGRAM_OPTIONS_PROG_OPT_NUM_RANGE_HPP

#include <string>
#include <vector>

#include <boost/any.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/program_options.hpp>

#include "cath/common/type_aliases.hpp"
#include "cath/common/type_traits.hpp"

namespace cath::common {

	/// \brief A type to invoke a program Options validation that rejects values out the range
	///        specified by the template arguments
	template <typename Num,          ///< The type in which the validation should occur
	          int64_t MinVal,        ///< The minimum permissible value before aborting validation
	          int64_t MaxVal,        ///< The maximum permissible value before aborting validation
	          typename ConvNum = Num ///< The type to which the value string should first be cast
	          >
	class prog_opt_num_range final {
	private:
		/// \brief The wrapped Num value
		Num value;

		static_assert( std::is_arithmetic_v<Num>,
			"prog_opt_num_range can only be instantiated for arithmetic types (eg int, double)");

	  public:
		/// \brief Ctor from Num
		///
		/// \param prm_value The value to set
		// clang-format off
		constexpr prog_opt_num_range( const Num &prm_value ) :value( prm_value ) { // NOLINT(google-explicit-constructor)
		}

		/// \brief Num conversion operator
		constexpr operator Num() const { // NOLINT(google-explicit-constructor)
			return value;
		}
		// clang-format on
	};

	/// \brief Provide Boost program_options validation for prog_opt_num_range
	template <typename Num,          ///< The type in which the validation should occur
	          int64_t MinVal,        ///< The minimum permissible value before aborting validation
	          int64_t MaxVal,        ///< The maximum permissible value before aborting validation
	          typename ConvNum = Num ///< The type to which the value string should first be cast
	          >
	void validate(boost::any    &prm_value,         ///< The value to populate
	              const str_vec &prm_value_strings, ///< The string values to validate
	              prog_opt_num_range<Num, MinVal, MaxVal, ConvNum> * dummy_var1, int /*dummy_var2*/
	              ) {
		using ponr_t = common::remove_cvref_t< std::remove_pointer_t< decltype( dummy_var1 ) > >;

		// Standard validate boilerplate:
		//  * Make sure no previous assignment to 'a' was made.
		//  * Extract the first string from 'prm_value_strings'.
		//    (If there is more than one string, it's an error, and exception will be thrown.)
		boost::program_options::check_first_occurrence( prm_value );
		const std::string &value_string = boost::program_options::get_single_string( prm_value_strings );

		// Attempt to lexical_cast value_string and if it fails, throw an invalid_option_value exception
		try {
			const Num result = boost::numeric_cast<Num>( boost::lexical_cast<ConvNum>( value_string ) );
			if ( boost::numeric_cast<int64_t>( result ) < MinVal || boost::numeric_cast<int64_t>( result ) > MaxVal ) {
				throw 0;
			}
			prm_value = ponr_t{ result };
		}
		catch (...) {
			BOOST_THROW_EXCEPTION( boost::program_options::invalid_option_value( value_string ) );
		}

		
	}

} // namespace cath::common

#endif // _CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_PROGRAM_OPTIONS_PROG_OPT_NUM_RANGE_HPP
