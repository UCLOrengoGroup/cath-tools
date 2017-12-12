/// \file
/// \brief The common type_aliases header

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

#ifndef _CATH_TOOLS_SOURCE_SRC_COMMON_COMMON_TYPE_ALIASES_H
#define _CATH_TOOLS_SOURCE_SRC_COMMON_COMMON_TYPE_ALIASES_H

#include <boost/optional/optional_fwd.hpp>
#include <boost/tuple/tuple.hpp>

#include <deque>
#include <iosfwd>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace cath { namespace common { template <typename> class clone_ptr; } }
namespace cath { namespace common { template <typename T> class vector_of_vector; } }

namespace cath {

	using ostream_ref                   = std::reference_wrapper<std::ostream>;
	using ostream_ref_opt               = boost::optional<ostream_ref>;

	using bool_size_str_tpl             = std::tuple<bool, size_t, std::string>;
	using bool_size_str_tpl_vec         = std::vector<bool_size_str_tpl>;

	using diff_diff_pair                = std::pair<ptrdiff_t, ptrdiff_t>;
	using size_size_doub_tpl            = std::tuple<size_t, size_t, double>;
	using size_size_doub_tpl_vec        = std::vector<size_size_doub_tpl>;
	using size_size_tpl                 = std::tuple<size_t, size_t>;
	using str_str_str_str_tpl           = std::tuple<std::string, std::string, std::string, std::string>;
	using str_str_str_str_tpl_vec       = std::vector<str_str_str_str_tpl>;

	using size_bool_pair                = std::pair<size_t, bool>;
	using size_bool_pair_vec            = std::vector<size_bool_pair>;

	using str_set                       = std::set<std::string>;

	using str_doub_pair                 = std::pair<std::string, double>;
	using str_doub_pair_vec             = std::vector<str_doub_pair>;

	using str_bool_pair                 = std::pair<std::string, bool>;
	using str_bool_pair_vec             = std::vector<str_bool_pair>;

	/// \brief A type alias for a reference_wrapper of const string
	using string_cref                   = std::reference_wrapper<const std::string>;

	/// \brief Type alias for a tuple of string, size_t and double
	using str_size_doub_tpl             = std::tuple<std::string, size_t, double>;

	/// \brief Type alias for a vector of str_size_doub_tpl
	using str_size_doub_tpl_vec         = std::vector<str_size_doub_tpl>;

	using bool_deq                      = std::deque<bool>;
	using bool_deq_itr                  = bool_deq::iterator;
	using bool_deq_citr                 = bool_deq::const_iterator;
	using bool_deq_vec                  = std::vector<bool_deq>;

	using char_vec                      = std::vector<char>;

	using doub_vec                      = std::vector<double>;
	using doub_vec_vec                  = std::vector<doub_vec>;

	/// \brief Type alias for an optional double
	using doub_opt                      = boost::optional<double>;

	using diff_vec                      = std::vector<ptrdiff_t>;
	using diff_vec_vec                  = std::vector<diff_vec>;

	/// \brief Alias template for a vector of reference_wrappers of T
	template <typename T>
	using ref_vec                       = std::vector<std::reference_wrapper<T>>;

	template <typename T>
	using uptr_vec                      = std::vector<std::unique_ptr<T>>;

	using size_opt                      = boost::optional<size_t>;
	using size_opt_vec                  = std::vector<boost::optional<size_t> >;

	//using uint_uint_pair                = std::pair<unsigned int, unsigned int>;
	//using uint_uint_pair_vec            = std::vector<uint_uint_pair>;

	using str_str_pair                  = std::pair<std::string, std::string>;
	using str_str_pair_vec              = std::vector<str_str_pair>;

	using str_str_pair_doub_map         = std::map  <str_str_pair, double>;
	using str_str_pair_doub_pair        = std::pair <str_str_pair, double>;

	using str_str_pair_bool_map         = std::map  <str_str_pair, bool>;
	using str_str_pair_bool_pair        = std::pair <str_str_pair, bool>;

	using str_str_str_pair_map          = std::map<std::string, str_str_pair>;

	using str_str_pair_size_map         = std::map<str_str_pair, size_t>;

	using size_size_pair                = std::pair<size_t, size_t>;
	using size_size_pair_vec            = std::vector<size_size_pair>;
	using size_size_pair_doub_map       = std::map<size_size_pair, double>;
	using size_size_pair_doub_map_value = size_size_pair_doub_map::value_type;

	using str_str_map                   = std::map<std::string, std::string>;

	using size_size_pair_opt            = boost::optional<size_size_pair>;

	using str_doub_map                  = std::map <std::string, double>;

	using str_size_pair                 = std::pair<std::string, size_t>;
	using str_size_map                  = std::map <std::string, size_t>;
	using str_size_pair_vec             = std::vector<str_size_pair>;

	using str_opt                       = boost::optional<std::string>;

	/// \brief A type alias for a vector of name_set objects
	using str_opt_vec                   = std::vector<str_opt>;

	using doub_doub_pair                = std::pair<double, double>;
	using doub_doub_pair_vec            = std::vector<doub_doub_pair>;
	using doub_doub_pair_vec_itr        = doub_doub_pair_vec::iterator;
	using doub_doub_pair_vec_citr       = doub_doub_pair_vec::const_iterator;

	using doub_size_pair                = std::pair<double, size_t>;
	using size_doub_pair                = std::pair<size_t, double>;

	using str_citr                      = std::string::const_iterator;

	/// \brief A type alias for a pair of str_ctirs
	using str_citr_str_citr_pair        = std::pair<str_citr, str_citr>;

	using str_vec                       = std::vector<std::string>;
	using str_vec_citr                  = str_vec::const_iterator;
	using str_vec_vec                   = std::vector<str_vec>;

	using str_str_vec_pair              = std::pair<std::string, str_vec>;
	using str_str_vec_map               = std::map <std::string, str_vec>;

	using str_deq                       = std::deque<std::string>;
	using str_deq_citr                  = str_vec::const_iterator;

	using str_size_type                 = std::string::size_type;

	using size_vec                      = std::vector<size_t>;
	using size_vec_vec                  = std::vector<size_vec>;

	using size_vec_size_vec_pair        = std::pair<size_vec, size_vec>;

	using size_deq                      = std::deque<size_t>;

	using size_size_map                 = std::map<size_t, size_t>;

	using char_opt                      = boost::optional<char>;

	using int_vec                       = std::vector<int>;
	using int_vec_vec                   = std::vector<int_vec>;

	using size_set                      = std::set<size_t>;

	using size_size_vec_pair            = std::pair<size_t, size_vec>;
	using size_size_vec_map             = std::map<size_t, size_vec>;

	using size_size_size_size_tpl       = std::tuple<size_t, size_t, size_t, size_t>;
} // namespace cath

namespace cath {
	namespace common {

		using bool_vec_of_vec = vector_of_vector<bool>; // WARNING: Remember std has a stupid template specialisation for vector<bool>

		/// \brief Alias template for a vector of clone_ptr of T
		template <typename T>
		using clptr_vec = std::vector< common::clone_ptr< T > >;

	} // namespace common
} // namespace cath

namespace cath {
	class protein;

	using float_score_type             = double;
	using float_score_vec              = std::vector<float_score_type>;
	using float_score_vec_vec          = std::vector<float_score_vec>;
	using float_score_float_score_pair = std::pair<float_score_type, float_score_type>;

	using score_opt                    = boost::optional<double>;
	using score_opt_vec                = std::vector<score_opt>;
	using score_opt_vec_vec            = std::vector<score_opt_vec>;

	using prot_prot_pair               = std::pair<protein, protein>;

	/// \brief The type of the scores used in aligning with dynamic-programming etc
	///
	/// \todo Switch this to float/double; if that's too broad-brush then templatise on score_type where appropriate
	using score_type                   = int;
	using score_vec                    = std::vector<score_type>;
	using score_vec_vec                = std::vector<score_vec>;
	using score_vec_of_vec             = common::vector_of_vector<score_type>;

	using str_str_score_tpl            = std::tuple<std::string, std::string, score_type>;
} // namespace cath

#endif
