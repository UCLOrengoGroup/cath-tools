/// \file
/// \brief The length_getter class header

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

#ifndef _CATH_TOOLS_SOURCE_UNI_SCORE_LENGTH_GETTER_LENGTH_GETTER_HPP
#define _CATH_TOOLS_SOURCE_UNI_SCORE_LENGTH_GETTER_LENGTH_GETTER_HPP

#include <boost/logic/tribool_fwd.hpp>
#include <boost/ptr_container/ptr_vector.hpp>

#include "score/aligned_pair_score/aligned_pair_score.hpp"

namespace cath {
	namespace align {
		class alignment;
	} // namespace align
	class protein;

	namespace score {

		/// \brief TODOCUMENT
		enum class length_getter_category : char {
			LONGER,  ///< TODOCUMENT
			SHORTER, ///< TODOCUMENT
			OTHER    ///< TODOCUMENT
		};

		/// \brief TODOCUMENT
		using length_getter_category_vec = std::vector<length_getter_category>;



		/// \brief TODOCUMENT
		class length_getter : private cath::common::polymorphic_less_than_comparable<length_getter>,
		                      private boost::equivalent<length_getter,
		                              boost::totally_ordered<length_getter> > {
		private:
			/// \brief Pure virtual method with which each concrete length_getter must define how to create a clone of itself
			virtual std::unique_ptr<length_getter> do_clone() const = 0;

			virtual boost::logic::tribool do_higher_is_better() const = 0;

			/// \brief TODOCUMENT
			virtual size_t do_get_length(const align::alignment &,
			                             const protein &,
			                             const protein &) const = 0;

			/// \brief TODOCUMENT
			virtual length_getter_category do_get_length_getter_category() const = 0;

			/// \brief TODOCUMENT
			virtual std::string do_id_name() const = 0;

			/// \brief TODOCUMENT
			virtual str_bool_pair_vec do_short_name_suffixes() const = 0;

			/// \brief TODOCUMENT
			virtual std::string do_long_name() const = 0;

			/// \brief TODOCUMENT
			virtual std::string do_description() const = 0;

			/// \brief TODOCUMENT
			virtual const std::string do_description_brackets_string() const = 0;

			/// \brief TODOCUMENT
			virtual bool do_less_than_with_same_dynamic_type(const length_getter &) const = 0;

		public:
			length_getter() = default;
			virtual ~length_getter() noexcept = default;

			length_getter(const length_getter &) = default;
			length_getter(length_getter &&) noexcept = default;
			length_getter & operator=(const length_getter &) = default;
			length_getter & operator=(length_getter &&) noexcept = default;

			std::unique_ptr<length_getter> clone() const;

			boost::logic::tribool higher_is_better() const;

			size_t get_length(const align::alignment &,
			                  const protein &,
			                  const protein &) const;

			length_getter_category get_length_getter_category() const;

			std::string human_friendly_short_name() const;
			std::string full_short_name() const;

			std::string id_name() const;
			str_bool_pair_vec short_name_suffixes() const;
			std::string long_name() const;
			std::string description() const;
			const std::string description_brackets_string() const;

			bool less_than_with_same_dynamic_type(const length_getter &) const;
		};

		str_bool_pair_vec length_getter_as_short_name_suffixes(const length_getter &,
		                                                       const bool & = true);
		str_bool_pair_vec length_getter_as_short_name_suffixes(const length_getter &,
		                                                       const length_getter_category_vec &);

		boost::ptr_vector<length_getter> get_all_length_getters();

		score_value get_length_score(const length_getter &,
		                             const align::alignment &,
		                             const protein &,
		                             const protein &);

		/// \brief Function to make length_getter meet the Clonable concept (used in ptr_container)
        ///
        /// NOTE: Don't call this yourself. Call the object's clone() method instead because that returns a
        ///       smart pointer rather than the raw pointer this has to return to meet the Clonable concept.
        ///
        /// This gets the smart pointer from the clone() method and then calls release on it.
        ///
        /// \returns A raw pointer to a new copy of the length_getter argument, with the same dynamic type.
        ///          The caller is responsible for deleting this new object.
        inline length_getter * new_clone(const length_getter &prm_length_getter ///< The length_getter to clone
                                         ) {
                return prm_length_getter.clone().release();
        }
	} // namespace score
} // namespace cath

#endif
