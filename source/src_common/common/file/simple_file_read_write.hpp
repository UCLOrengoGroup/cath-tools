/// \file
/// \brief The simple_file_read_write header

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

#ifndef _CATH_TOOLS_SOURCE_SRC_COMMON_COMMON_FILE_SIMPLE_FILE_READ_WRITE_HPP
#define _CATH_TOOLS_SOURCE_SRC_COMMON_COMMON_FILE_SIMPLE_FILE_READ_WRITE_HPP

#include <boost/concept/assert.hpp>
#include <boost/concept_archetype.hpp>
#include <boost/core/ignore_unused.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/range.hpp>
#include <boost/range/iterator_range.hpp>

#include "common/cpp14/cbegin_cend.hpp"
#include "common/cpp17/apply.hpp"
#include "common/exception/out_of_range_exception.hpp"
#include "common/file/open_fstream.hpp"
#include "common/type_aliases.hpp"

#include <fstream>

namespace cath {
	namespace common {
		namespace detail {

			/// \brief Stream parts from the ctor-specified istream to the specified argument references
			///        in the correct order
			///
			/// Motivation: to be used in apply() for populating a tuple from an istream.
			class tuple_parts_istreamer final {
			private:
				/// \brief Reference to the istream from which the parts should be read
				std::istream &the_istream;

				/// \brief Stream an individual part from the_istream
				template <typename T>
				void istream_part(T &prm_value ///< The part to be populated from the istream
				                  ) {
					the_istream >> prm_value;
				}

			public:
				/// \brief Ctor from istream from which the parts should be read
				explicit tuple_parts_istreamer(std::istream &prm_istream ///< The istream from which the parts should be read
				                               ) : the_istream( prm_istream ) {
				}

				/// \brief Populate the arguments from the_istream
				///
				/// This is the function that's called by transform tuple
				template <typename... Ks>
				void operator()(Ks &... prm_parts ///< Non-const references to the parts to be populated
				                ) {
					// Use array initialisation to ensure the correct order of evaluation
					// (function calls don't guarantee evaluation order of parts and
					//  though initializer_lists should, GCC had a bug undermining around
					//  versions v4.7 and v4.8)
					bool b[] = { true, ( istream_part( prm_parts ), true )... };
					boost::ignore_unused( b );
				}
			};

			/// \brief Explicit specialisation of istream_part for bool that uses boolalpha for
			///        the single element
			template <>
			inline void tuple_parts_istreamer::istream_part<bool>(bool &prm_value ///< The bool value to be populated from the_istream
			                                                      ) {
				the_istream >> std::boolalpha >> prm_value >> std::noboolalpha;
			}

			/// \brief Template for reading a line from an istream into an arbitrary type T
			///
			/// \tparam T must be default constructible
			template <typename T>
			struct type_line_reader final {

				/// \brief Read a line from the specified istream into a T and return the results
				static T read_line(std::istream &prm_is ///< The istream containing the line to be read
				                   ) {
					T value;
					if ( prm_is >> value ) {
						return value;
					}
					BOOST_THROW_EXCEPTION(out_of_range_exception("Unable to parse value from input stream"));
				}
			};

			/// \brief Partial specialisation for reading a line from a tuple
			///
			/// \tparam std::tuple<Ts...> must be default constructible
			///         (which presumably requires all Ts... are default constructible)
			///
			/// Unlike pair<> (below) this is not currently able to handle any of the Ts
			/// being complex (pair, tuple) types themselves.
			///
			/// \todo If there's need, consider attempting to achieve that recursion
			///       allowing, eg, tuple<<tuple...>, ...>
			template <typename... Ts>
			struct type_line_reader<std::tuple<Ts...>> final {

				/// \brief Read a line from the specified istream into a T
				static std::tuple<Ts...> read_line(std::istream &prm_is ///< The istream containing the line to be read
				                                   ) {
					std::tuple<Ts...> new_tuple;
					/// \TODO Come C++17, use ::std::apply
					::cath::common::apply( tuple_parts_istreamer( prm_is ), new_tuple );
					return new_tuple;
				}
			};

			/// \brief TODOCUMENT
			///
			/// \tparam T must be default constructible
			/// \tparam U must be default constructible
			///
			/// This recursively uses 
			template <typename T, typename U>
			struct type_line_reader<std::pair<T, U>> final {

				/// \brief TODOCUMENT
				static std::pair<T, U> read_line(std::istream &prm_is ///< The istream containing the line to be read
				                                 ) {
					std::pair<T, U> value;
					value.first  = type_line_reader<T>::read_line( prm_is );
					value.second = type_line_reader<U>::read_line( prm_is );
					return value;
				}
			};

			/// \brief TODOCUMENT
			template <>
			struct type_line_reader<bool> final {

				/// \brief TODOCUMENT
				static bool read_line(std::istream &prm_is ///< The istream containing the line to be read
				                      ) {
					bool value = false;
					if ( prm_is >> std::boolalpha >> value >> std::noboolalpha ) {
						return value;
					}
					BOOST_THROW_EXCEPTION(out_of_range_exception("Unable to parse bool value from input stream"));
				}
			};



			/// \brief TODOCUMENT
			class tuple_parts_ostreamer final {
			private:
				/// \brief TODOCUMENT
				std::ostream &the_ostream;

				/// \brief TODOCUMENT
				template <typename T>
				void ostream_part(const T &prm_value) {
					the_ostream <<  prm_value;
				}

				/// \brief TODOCUMENT
				void ostream_part(const bool &prm_value) {
					the_ostream << std::boolalpha << prm_value << std::noboolalpha;
				}

			public:
				/// \brief TODOCUMENT
				explicit tuple_parts_ostreamer(std::ostream &prm_ostream
				                               ) : the_ostream ( prm_ostream ) {
				}

				/// \brief TODOCUMENT
				template <typename T1, typename T2, typename... Ts>
				void operator()(const T1 &     prm_value_1,         ///< TODOCUMENT
				                const T2 &     prm_value_2,         ///< TODOCUMENT
				                const Ts & ... prm_values_3_onwards ///< TODOCUMENT
				                ) {
					ostream_part( prm_value_1 );
					the_ostream << " ";
					(*this)( prm_value_2, prm_values_3_onwards... );
				}

				/// \brief TODOCUMENT
				template <typename T>
				void operator()(const T &prm_value ///< TODOCUMENT
				                ) {
					ostream_part( prm_value );
				}
			};



			/// \brief TODOCUMENT
			class type_line_writer final {
			private:
				/// \brief TODOCUMENT
				template <typename T>
				void write_part(std::ostream &prm_os,   ///< TODOCUMENT
				                const T      &prm_value ///< TODOCUMENT
				                ) const {
					write_part( prm_os, std::make_tuple( prm_value) );
				}

				/// \brief TODOCUMENT
				template <typename T, typename U>
				void write_part(std::ostream          &prm_os,   ///< TODOCUMENT
				                const std::pair<T, U> &prm_value ///< TODOCUMENT
				                ) const {
					write_part( prm_os, prm_value.first  );
					prm_os << " ";
					write_part( prm_os, prm_value.second );
				}

				/// \brief TODOCUMENT
				template <typename... Ts>
				void write_part(std::ostream            &prm_os,   ///< TODOCUMENT
				                const std::tuple<Ts...> &prm_value ///< TODOCUMENT
				                ) const {
					/// \TODO Come C++17, use ::std::apply
					::cath::common::apply( tuple_parts_ostreamer( prm_os ), prm_value );
				}

			public:
				/// \brief TODOCUMENT
				template <typename T>
				void operator()(std::ostream &prm_os,   ///< TODOCUMENT
				                const T      &prm_value ///< TODOCUMENT
				                ) const {
					write_part( prm_os, prm_value );
					prm_os << "\n";
				}
			};


		} // namespace detail

		/// \brief TODOCUMENT
		template <typename T>
		std::vector<T> read_file(const boost::filesystem::path &prm_file ///< TODOCUMENT
		                         ) {
			std::ifstream input_stream;
			open_ifstream( input_stream, prm_file );
			std::vector<T> line_entries;
			std::string line_string;
			while ( getline( input_stream, line_string ) ) {
				std::istringstream line_ss( line_string );
				line_entries.push_back( detail::type_line_reader<T>::read_line( line_ss ) );
			}
			input_stream.close();
			return line_entries;
		}

		/// \brief TODOCUMENT
		template <typename ITER>
		void write_file(const boost::filesystem::path &prm_file,  ///< TODOCUMENT
		                const ITER                    &prm_begin, ///< TODOCUMENT
		                const ITER                    &prm_end    ///< TODOCUMENT
		                ) {
			BOOST_CONCEPT_ASSERT(( boost::forward_iterator_archetype<ITER> ));
			std::ofstream output_stream;
			open_ofstream( output_stream, prm_file );
			for (const auto &value : boost::iterator_range<ITER>( prm_begin, prm_end ) ) {
				detail::type_line_writer()( output_stream, value );
			}
			output_stream << std::flush;
			output_stream.close();
		}

		/// \brief TODOCUMENT
		template <typename SinglePassRange>
		void write_file(const boost::filesystem::path &prm_file, ///< TODOCUMENT
		                const SinglePassRange         &prm_range ///< TODOCUMENT
		                ) {
			BOOST_RANGE_CONCEPT_ASSERT(( boost::SinglePassRangeConcept<const SinglePassRange> ));
			cath::common::write_file(
				prm_file,
				common::cbegin( prm_range ),
				common::cend  ( prm_range )
			);
		}

		/// \brief Write a single string to a file
		inline void write_file(const boost::filesystem::path &prm_file,  ///< The file to which the string should be written
		                       const std::string             &prm_string ///< The string to write
		                       ) {
			cath::common::write_file(
				prm_file,
				str_vec{ { prm_string } }
			);
		}
	} // namespace common
} // namespace cath

#endif
