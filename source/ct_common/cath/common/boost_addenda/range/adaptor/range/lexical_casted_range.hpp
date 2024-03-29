/// \file
/// \brief The lexical_casted_range class header

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

#ifndef CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_BOOST_ADDENDA_RANGE_ADAPTOR_RANGE_LEXICAL_CASTED_RANGE_HPP
#define CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_BOOST_ADDENDA_RANGE_ADAPTOR_RANGE_LEXICAL_CASTED_RANGE_HPP

#include "cath/common/boost_addenda/range/adaptor/iterator/lexical_cast_itr.hpp"

namespace cath::common::detail {

	/// \brief Range wrapper over RNG whose dereference operator returns the result of performing lexical_cast<T>() on the
	///        dereferenced value from the original value
	///
	/// For a const lexical_casted_range, use a const RNG type.
	///
	/// This is implemented as an iterator_range of transform_iterator using lexical_cast_value
	template <typename T, typename RNG>
	class lexical_casted_range final : public boost::iterator_range<lexical_cast_itr<T, RNG> > {
	private:
		/// \brief Convenience type alias for lexical_cast_itr for T and RNG
		using lex_casted_iterator = lexical_cast_itr<T, RNG>;

		/// \brief Convenience type alias for the parent iterator_range class through which this is implemented
		using super               = boost::iterator_range<lex_casted_iterator>;

	public:
		explicit lexical_casted_range(const RNG &);
	 };

	/// \brief Ctor from a range
	template <typename T, typename RNG>
	lexical_casted_range<T, RNG>::lexical_casted_range(const RNG &prm_range ///< The range over which to apply this lexical_casted_range
	                                                   ) : super(
	                                                       	lex_casted_iterator( std::begin( prm_range ) ),
	                                                       	lex_casted_iterator( std::end  ( prm_range ) )
	                                                       ) {
	}

} // namespace cath::common::detail


#endif // CATH_TOOLS_SOURCE_CT_COMMON_CATH_COMMON_BOOST_ADDENDA_RANGE_ADAPTOR_RANGE_LEXICAL_CASTED_RANGE_HPP
