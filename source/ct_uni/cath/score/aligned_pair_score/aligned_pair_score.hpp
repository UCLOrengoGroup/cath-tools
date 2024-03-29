/// \file
/// \brief The aligned_pair_score class header

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

#ifndef CATH_TOOLS_SOURCE_CT_UNI_CATH_SCORE_ALIGNED_PAIR_SCORE_ALIGNED_PAIR_SCORE_HPP
#define CATH_TOOLS_SOURCE_CT_UNI_CATH_SCORE_ALIGNED_PAIR_SCORE_ALIGNED_PAIR_SCORE_HPP

#include <boost/logic/tribool_fwd.hpp>
#include <boost/operators.hpp>
#include <boost/serialization/nvp.hpp>

#include "cath/common/polymorphic_comparison/polymorphic_less_than_comparable.hpp"
#include "cath/common/type_aliases.hpp"
#include "cath/score/score_type_aliases.hpp"

#include <memory>

// clang-format off
namespace cath { class protein; }
namespace cath::align { class alignment; }
namespace cath::score::detail { str_aligned_pair_score_pmap get_aligned_pair_score_of_id_name(); }
// clang-format on

namespace cath::score {

	/// \brief Provide ABC interface for classes that calculate scores for pair alignments and their associated proteins
	class aligned_pair_score : private cath::common::polymorphic_less_than_comparable<aligned_pair_score>,
	                           private boost::equivalent<aligned_pair_score,
	                                   boost::totally_ordered<aligned_pair_score> > {
	private:
		friend class boost::serialization::access;

		template<class archive> void serialize(archive      &/*ar*/,     ///< TODOCUMENT
		                                       const size_t  /*version*/ ///< TODOCUMENT
		                                       ) {
		}

		/// \brief Pure virtual method with which each concrete policy must define how to create a clone of itself
		[[nodiscard]] virtual std::unique_ptr<aligned_pair_score> do_clone() const = 0;

		/// \brief Pure virtual method with which each concrete aligned_pair_score must define whether a higher score
		///        (rather than a lower score) generally reflects a better alignment and/or more similar structures.
		[[nodiscard]] virtual boost::logic::tribool do_higher_is_better() const = 0;

		/// \brief Pure virtual method with which each concrete aligned_pair_score must define the method of calculating a score for
		///        a pair alignment and two associated proteins.
		[[nodiscard]] virtual score_value do_calculate( const align::alignment &, const protein &, const protein & ) const = 0;

		/// \brief Pure virtual method with which each concrete aligned_pair_score must define a free text description, describing the score.
		[[nodiscard]] virtual std::string do_description() const = 0;

		/// \brief Pure virtual method with which each concrete aligned_pair_score must define a unique id_name (which cannot contain a full-stop)
		[[nodiscard]] virtual std::string do_id_name() const = 0;

		/// \brief Pure virtual method with which each concrete aligned_pair_score must define the strings to be appended to id_name()
		///        with full-stops to the get the full short name. Each string also comes with a bool that specifies
		///        whether the string should be included in human-friendly short_name (ie it indicates a non-default value)
		///        as well as the full short_name.
		[[nodiscard]] virtual str_bool_pair_vec do_short_name_suffixes() const = 0;

		/// \brief TODOCUMENT
		[[nodiscard]] virtual std::string do_long_name() const;

		/// \brief TODOCUMENT
		[[nodiscard]] virtual std::string do_reference() const;

		//			/// \brief Pure virtual method with which each concrete aligned_pair_score must define how to build an aligned_pair_score
		//			///        of that concrete type from a short_name_spec string
		//			virtual std::unique_ptr<aligned_pair_score> do_build_from_short_name_spec(const std::string &) const = 0;

		/// \brief TODOCUMENT
		[[nodiscard]] virtual bool do_less_than_with_same_dynamic_type( const aligned_pair_score & ) const = 0;

	  public:
		aligned_pair_score() = default;
		virtual ~aligned_pair_score() noexcept = default;
		[[nodiscard]] std::unique_ptr<aligned_pair_score> clone() const;

		aligned_pair_score(const aligned_pair_score &) = default;
		aligned_pair_score(aligned_pair_score &&) noexcept = default;
		aligned_pair_score & operator=(const aligned_pair_score &) = default;
		aligned_pair_score & operator=(aligned_pair_score &&) noexcept = default;

		[[nodiscard]] boost::logic::tribool higher_is_better() const;
		[[nodiscard]] score_value calculate( const align::alignment &, const protein &, const protein & ) const;

		[[nodiscard]] std::string       id_name() const;
		[[nodiscard]] str_bool_pair_vec short_name_suffixes() const;
		[[nodiscard]] std::string       human_friendly_short_name() const;
		[[nodiscard]] std::string       full_short_name() const;
		[[nodiscard]] std::string       long_name() const;
		[[nodiscard]] std::string       description() const;
		[[nodiscard]] std::string       reference() const;

		// std::unique_ptr<aligned_pair_score> build_from_short_name_spec(const std::string &) const;

		[[nodiscard]] bool less_than_with_same_dynamic_type( const aligned_pair_score & ) const;
	};

//	std::unique_ptr<aligned_pair_score> make_aligned_pair_score_from_full_short_name(const std::string &);

	std::ostream & operator<<(std::ostream &,
	                          const aligned_pair_score &);

	/// \brief Function to make aligned_pair_score meet the Clonable concept (used in ptr_container)
	///
	/// NOTE: Don't call this yourself. Call the object's clone() method instead because that returns a
	///       smart pointer rather than the raw pointer this has to return to meet the Clonable concept.
	///
	/// This gets the smart pointer from the clone() method and then calls release on it.
	///
	/// \returns A raw pointer to a new copy of the aligned_pair_score argument, with the same dynamic type.
	///          The caller is responsible for deleting this new object.
	inline aligned_pair_score * new_clone(const aligned_pair_score &prm_aligned_pair_score ///< The aligned_pair_score to clone
	                                      ) {
		return prm_aligned_pair_score.clone().release();
	}

} // namespace cath::score

#endif // CATH_TOOLS_SOURCE_CT_UNI_CATH_SCORE_ALIGNED_PAIR_SCORE_ALIGNED_PAIR_SCORE_HPP
