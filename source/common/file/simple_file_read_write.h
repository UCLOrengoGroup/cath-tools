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

#ifndef SIMPLE_FILE_READ_WRITE_H_INCLUDED
#define SIMPLE_FILE_READ_WRITE_H_INCLUDED

#include <boost/concept/assert.hpp>
#include <boost/concept_archetype.hpp>
#include <boost/core/ignore_unused.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/range.hpp>
#include <boost/range/iterator_range.hpp>

#include "common/algorithm/transform_tuple.h"
#include "common/cpp14/cbegin_cend.h"
#include "common/file/open_fstream.h"
#include "common/type_aliases.h"
#include "exception/out_of_range_exception.h"

#include <fstream>

namespace cath {
	namespace common {
		namespace detail {

			/// \brief Stream parts from the ctor-specified istream to the specified argument references
			///        in the correct order
			///
			/// Motivation: to be used in transform_tuple() for populating a tuple from an istream.
			class tuple_parts_istreamer final {
			private:
				/// \brief Reference to the istream from which the parts should be read
				std::istream &the_istream;

				/// \brief Stream an individual part from the_istream
				template <typename T>
				void istream_part(T &arg_value ///< The part to be populated from the istream
				                  ) {
					the_istream >> arg_value;
				}

			public:
				/// \brief Ctor from istream from which the parts should be read
				tuple_parts_istreamer(std::istream &arg_istream ///< The istream from which the parts should be read
				                      ) : the_istream( arg_istream ) {
				}

				/// \brief Populate the arguments from the_istream
				///
				/// This is the function that's called by transform tuple
				template <typename... Ks>
				void operator()(Ks &... arg_parts ///< Non-const references to the parts to be populated
				                ) {
					// Use array initialisation to ensure the correct order of evaluation
					// (function calls don't guarantee evaluation order of parts and
					//  though initializer_lists should, GCC had a bug undermining around
					//  versions v4.7 and v4.8)
					bool b[] = { true, ( istream_part( arg_parts ), true )... };
					boost::ignore_unused( b );
				}
			};

			/// \brief Explicit specialisation of istream_part for bool that uses boolalpha for
			///        the single element
			template <>
			inline void tuple_parts_istreamer::istream_part<bool>(bool &arg_value ///< The bool value to be populated from the_istream
			                                                      ) {
				the_istream >> std::boolalpha >> arg_value >> std::noboolalpha;
			}

			/// \brief Template for reading a line from an istream into an arbitrary type T
			///
			/// \tparam T must be default constructible
			template <typename T>
			struct type_line_reader final {

				/// \brief Read a line from the specified istream into a T and return the results
				static T read_line(std::istream &arg_is ///< The istream containing the line to be read
				                   ) {
					T value;
					if ( arg_is >> value ) {
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
				static std::tuple<Ts...> read_line(std::istream &arg_is ///< The istream containing the line to be read
				                                   ) {
					std::tuple<Ts...> new_tuple;
					common::transform_tuple( new_tuple, tuple_parts_istreamer( arg_is ) );
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
				static std::pair<T, U> read_line(std::istream &arg_is ///< The istream containing the line to be read
				                                 ) {
					std::pair<T, U> value;
					value.first  = type_line_reader<T>::read_line( arg_is );
					value.second = type_line_reader<U>::read_line( arg_is );
					return value;
				}
			};

			/// \brief TODOCUMENT
			template <>
			struct type_line_reader<bool> final {

				/// \brief TODOCUMENT
				static bool read_line(std::istream &arg_is ///< The istream containing the line to be read
				                      ) {
					bool value = false;
					if ( arg_is >> std::boolalpha >> value >> std::noboolalpha ) {
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
				void ostream_part(const T &arg_value) {
					the_ostream <<  arg_value;
				}

				/// \brief TODOCUMENT
				void ostream_part(const bool &arg_value) {
					the_ostream << std::boolalpha << arg_value << std::noboolalpha;
				}

			public:
				/// \brief TODOCUMENT
				tuple_parts_ostreamer(std::ostream &arg_ostream
				                      ) : the_ostream ( arg_ostream ) {
				}

				/// \brief TODOCUMENT
				template <typename T1, typename T2, typename... Ts>
				void operator()(const T1 &     arg_value_1,         ///< TODOCUMENT
				                const T2 &     arg_value_2,         ///< TODOCUMENT
				                const Ts & ... arg_values_3_onwards ///< TODOCUMENT
				                ) {
					ostream_part( arg_value_1 );
					the_ostream << " ";
					(*this)( arg_value_2, arg_values_3_onwards... );
				}

				/// \brief TODOCUMENT
				template <typename T>
				void operator()(const T &arg_value ///< TODOCUMENT
				                ) {
					ostream_part( arg_value );
				}
			};



			/// \brief TODOCUMENT
			class type_line_writer final {
			private:
				/// \brief TODOCUMENT
				template <typename T>
				void write_part(std::ostream &arg_os,   ///< TODOCUMENT
				                const T      &arg_value ///< TODOCUMENT
				                ) const {
					write_part( arg_os, std::make_tuple( arg_value) );
				}

				/// \brief TODOCUMENT
				template <typename T, typename U>
				void write_part(std::ostream          &arg_os,   ///< TODOCUMENT
				                const std::pair<T, U> &arg_value ///< TODOCUMENT
				                ) const {
					write_part( arg_os, arg_value.first  );
					arg_os << " ";
					write_part( arg_os, arg_value.second );
				}

				/// \brief TODOCUMENT
				template <typename... Ts>
				void write_part(std::ostream            &arg_os,   ///< TODOCUMENT
				                const std::tuple<Ts...> &arg_value ///< TODOCUMENT
				                ) const {
					common::transform_tuple( arg_value, tuple_parts_ostreamer( arg_os ) );
				}

			public:
				/// \brief TODOCUMENT
				template <typename T>
				void operator()(std::ostream &arg_os,   ///< TODOCUMENT
				                const T      &arg_value ///< TODOCUMENT
				                ) const {
					write_part( arg_os, arg_value );
					arg_os << "\n";
				}
			};


		}

		/// \brief TODOCUMENT
		template <typename T>
		std::vector<T> read_file(const boost::filesystem::path &arg_file ///< TODOCUMENT
		                         ) {
			std::ifstream input_stream;
			open_ifstream( input_stream, arg_file );
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
		void write_file(const boost::filesystem::path &arg_file,  ///< TODOCUMENT
		                const ITER                    &arg_begin, ///< TODOCUMENT
		                const ITER                    &arg_end    ///< TODOCUMENT
		                ) {
			BOOST_CONCEPT_ASSERT(( boost::forward_iterator_archetype<ITER> ));
			std::ofstream output_stream;
			open_ofstream( output_stream, arg_file );
			for (const auto &value : boost::iterator_range<ITER>( arg_begin, arg_end ) ) {
				detail::type_line_writer()( output_stream, value );
			}
			output_stream << std::flush;
			output_stream.close();
		}

		/// \brief TODOCUMENT
		template <typename SinglePassRange>
		void write_file(const boost::filesystem::path &arg_file, ///< TODOCUMENT
		                const SinglePassRange         &arg_range ///< TODOCUMENT
		                ) {
			BOOST_RANGE_CONCEPT_ASSERT(( boost::SinglePassRangeConcept<const SinglePassRange> ));
			cath::common::write_file(
				arg_file,
				common::cbegin( arg_range ),
				common::cend  ( arg_range )
			);
		}

		/// \brief Write a single string to a file
		inline void write_file(const boost::filesystem::path &arg_file,  ///< The file to which the string should be written
		                       const std::string             &arg_string ///< The string to write
		                       ) {
			cath::common::write_file(
				arg_file,
				str_vec{ { arg_string } }
			);
		}
	}
}

#endif
