/// \file
/// \brief The length_getter_enum header

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

#ifndef CATH_TOOLS_SOURCE_CT_UNI_CATH_SCORE_LENGTH_GETTER_LENGTH_GETTER_ENUM_HPP
#define CATH_TOOLS_SOURCE_CT_UNI_CATH_SCORE_LENGTH_GETTER_LENGTH_GETTER_ENUM_HPP

#include <type_traits>

#include <boost/mpl/find_if.hpp>
#include <boost/mpl/vector.hpp>

// clang-format off
namespace cath::score { class geometric_mean_length_getter; }
namespace cath::score { class length_of_first_getter; }
namespace cath::score { class length_of_longer_getter; }
namespace cath::score { class length_of_second_getter; }
namespace cath::score { class length_of_shorter_getter; }
namespace cath::score { class mean_length_getter; }
namespace cath::score { class num_aligned_length_getter; }
// clang-format on

namespace cath::score::detail {

	/// \brief TODOCUMENT
	enum class length_getter_enum : char {
		FIRST,          ///< TODOCUMENT
		GEOMETRIC_MEAN, ///< TODOCUMENT
		LONGER,         ///< TODOCUMENT
		MEAN,           ///< TODOCUMENT
		NUM_ALIGNED,    ///< TODOCUMENT
		SECOND,         ///< TODOCUMENT
		SHORTER         ///< TODOCUMENT
	};

	/// \brief Compile type lookup table between length_getter types and the corresponding length_getter_enum values
	using enums_of_length_getter_types = boost::mpl::vector<
		boost::mpl::pair< length_of_first_getter,       std::integral_constant< length_getter_enum, length_getter_enum::FIRST          > >,
		boost::mpl::pair< geometric_mean_length_getter, std::integral_constant< length_getter_enum, length_getter_enum::GEOMETRIC_MEAN > >,
		boost::mpl::pair< length_of_longer_getter,      std::integral_constant< length_getter_enum, length_getter_enum::LONGER         > >,
		boost::mpl::pair< mean_length_getter,           std::integral_constant< length_getter_enum, length_getter_enum::MEAN           > >,
		boost::mpl::pair< num_aligned_length_getter,    std::integral_constant< length_getter_enum, length_getter_enum::NUM_ALIGNED    > >,
		boost::mpl::pair< length_of_second_getter,      std::integral_constant< length_getter_enum, length_getter_enum::SECOND         > >,
		boost::mpl::pair< length_of_shorter_getter,     std::integral_constant< length_getter_enum, length_getter_enum::SHORTER        > >
	>;

	/// \brief Compile-type metafunction to find the type associated with a length_getter_enum value
	///
	/// Use via length_getter_of_enum_t
	template <length_getter_enum E>
	class length_getter_of_enum_impl final {
	private:
		using enum_integral = std::integral_constant<length_getter_enum, E>;
		using second_mf     = typename boost::mpl::second< ::boost::mpl::placeholders::_1 >;
		using find_iter     = typename boost::mpl::find_if<enums_of_length_getter_types, ::std::is_same< second_mf, enum_integral > >::type ;

	public:
		using type          = typename boost::mpl::first<typename boost::mpl::deref<find_iter>::type>::type;
	};

	/// \brief Interface to compile-type metafunction to find the type associated with a length_getter_enum value
	///
	/// Use like this: `using my_length_getter_type = length_getter_of_enum_t<length_getter_enum::MEAN>;`
	template <length_getter_enum E>
	using length_getter_of_enum_t = typename length_getter_of_enum_impl<E>::type;

	/// \brief Compile-type metafunction to find the enum associated with a length_getter type
	///
	/// Use via ?????
	template <typename T>
	class enum_of_length_getter_impl final {
	private:
		using first_mf      = typename boost::mpl::first< ::boost::mpl::placeholders::_1 >;
		using find_iter     = typename boost::mpl::find_if<enums_of_length_getter_types, ::std::is_same< first_mf, T> >::type ;
		using enum_integral = typename boost::mpl::second<typename boost::mpl::deref<find_iter>::type>::type;

	public:
		static constexpr length_getter_enum value = enum_integral();
	};

	template <typename T>
	constexpr length_getter_enum enum_of_length_getter_v = enum_of_length_getter_impl<T>::value;

} // namespace cath::score::detail

#endif // CATH_TOOLS_SOURCE_CT_UNI_CATH_SCORE_LENGTH_GETTER_LENGTH_GETTER_ENUM_HPP
