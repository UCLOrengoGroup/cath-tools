/// \file
/// \brief The res_pair_keyer i/o header

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

#ifndef _CATH_TOOLS_SOURCE_UNI_SCAN_RES_PAIR_KEYER_DETAIL_RES_PAIR_KEYER_IO_H
#define _CATH_TOOLS_SOURCE_UNI_SCAN_RES_PAIR_KEYER_DETAIL_RES_PAIR_KEYER_IO_H

#include <boost/algorithm/string/join.hpp>
#include <boost/lexical_cast.hpp>

#include "common/cpp17/apply.hpp"
#include "common/type_aliases.hpp"

#include <iomanip>
#include <sstream>

namespace cath {
	namespace scan {
		namespace detail {

			/// \brief TODOCUMENT
			template <typename T>
			inline std::string output_key_part_impl(const T &arg_value ///< TODOCUMENT
			                                        ) {
				std::ostringstream out_ss;
				out_ss << std::right << std::setw( 3 ) << arg_value;
				return out_ss.str();
			}

			/// \brief TODOCUMENT
			template <typename T>
			inline std::string output_key_part(const T &arg_value ///< TODOCUMENT
			                                   ) {
				return output_key_part_impl( arg_value );
			}

			/// \brief TODOCUMENT
			template <>
			inline std::string output_key_part<uint8_t>(const uint8_t &arg_value ///< TODOCUMENT
			                                            ) {
				return output_key_part_impl( static_cast<size_t>( arg_value ) );
			}

			struct key_parts_outputter final {
				/// \brief TODOCUMENT
				template <typename... Ks>
				std::string operator()(const Ks &... arg_key_parts ///< TODOCUMENT
				                       ) {
					const str_vec key_part_strings = { output_key_part( arg_key_parts )... };
					return "res_pair_key[ " + boost::algorithm::join( key_part_strings, ", " ) + " ]";
				}
			};

			/// \brief TODOCUMENT
			template <typename... Ks>
			std::string output_key(const std::tuple<Ks...> &arg_key ///< TODOCUMENT
			                       ) {
				return common::apply( key_parts_outputter(), arg_key );
			}



			/// \brief TODOCUMENT
			struct keyer_parts_outputter final {
				/// \brief TODOCUMENT
				template <typename... KPs>
				std::string operator()(const KPs &... arg_keyer_parts ///< TODOCUMENT
				                       ) {
					const str_vec keyer_part_strings = { arg_keyer_parts.get_name()... };
					return boost::algorithm::join( keyer_part_strings, ", " );
				}
			};

			/// \brief TODOCUMENT
			template <typename... KPs>
			std::string output_keyer_parts(const std::tuple<KPs...> &arg_keyer_parts ///< TODOCUMENT
			                               ) {
				return common::apply( keyer_parts_outputter(), arg_keyer_parts );
			}

		} // namespace detail
	} // namespace scan
} // namespace cath

#endif

