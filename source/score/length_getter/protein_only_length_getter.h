/// \file
/// \brief The protein_only_length_getter class header

/// \copyright
/// CATH Binaries - Protein structure comparison tools such as SSAP and SNAP
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

#ifndef PROTEIN_ONLY_LENGTH_GETTER_H_INCLUDED
#define PROTEIN_ONLY_LENGTH_GETTER_H_INCLUDED

#include <boost/ptr_container/ptr_vector.hpp>

#include "score/length_getter/length_getter.h"

namespace cath {
	namespace score {

		/// \brief TODOCUMENT
		class protein_only_length_getter : public length_getter {
		private:
			/// \brief TODOCUMENT
			virtual std::unique_ptr<protein_only_length_getter> do_protein_only_clone() const = 0;

			/// \brief TODOCUMENT
			virtual size_t do_get_length(const protein &,
			                             const protein &) const = 0;

			/// \brief TODOCUMENT
			virtual std::string do_get_choice_adjective() const = 0;



			virtual std::unique_ptr<length_getter> do_clone() const override final;

			virtual boost::logic::tribool do_higher_is_better() const override final;

			virtual size_t do_get_length(const align::alignment &,
			                             const protein &,
			                             const protein &) const override final;

			virtual std::string do_id_name() const override final;

			virtual str_bool_pair_vec do_short_name_suffixes() const override final;

			virtual std::string do_long_name() const override final;

			virtual std::string do_description() const override final;

//			virtual std::string do_short_suffix_string() const override final;

//			virtual std::string do_long_suffix_string() const override final;

			virtual const std::string do_description_brackets_string() const override final;

		public:
			virtual ~protein_only_length_getter() noexcept = default;

			std::unique_ptr<protein_only_length_getter> protein_only_clone() const;

			size_t get_prot_only_length(const protein &,
			                            const protein &) const;

			std::string get_choice_adjective() const;
		};

		boost::ptr_vector<protein_only_length_getter> get_all_protein_only_length_getters();

		score_value get_length_score(const protein_only_length_getter &,
		                             const protein &,
		                             const protein &);

		/// \brief Function to make protein_only_length_getter meet the Clonable concept (used in ptr_container)
        ///
        /// NOTE: Don't call this yourself. Call the object's clone() method instead because that returns a
        ///       smart pointer rather than the raw pointer this has to return to meet the Clonable concept.
        ///
        /// This gets the smart pointer from the clone() method and then calls release on it.
        ///
        /// \returns A raw pointer to a new copy of the protein_only_length_getter argument, with the same dynamic type.
        ///          The caller is responsible for deleting this new object.
        inline protein_only_length_getter * new_clone(const protein_only_length_getter &arg_protein_only_length_getter ///< The protein_only_length_getter to clone
                                                      ) {
                return arg_protein_only_length_getter.protein_only_clone().release();
        }
	}
}

#endif