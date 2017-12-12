/// \file
/// \brief The superposition_outputter class header

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

#ifndef _CATH_TOOLS_SOURCE_UNI_OUTPUTTER_SUPERPOSITION_OUTPUTTER_SUPERPOSITION_OUTPUTTER_H
#define _CATH_TOOLS_SOURCE_UNI_OUTPUTTER_SUPERPOSITION_OUTPUTTER_SUPERPOSITION_OUTPUTTER_H

#include <boost/utility/string_ref_fwd.hpp>

#include <iosfwd>
#include <memory>

namespace cath { namespace sup { class superposition_context; } }

namespace cath {
	namespace opts {

		/// \brief TODOCUMENT
		class superposition_outputter {
		private:
			/// \brief Pure virtual method with which each concrete superposition_outputter must define how to create a clone of itself
			virtual std::unique_ptr<superposition_outputter> do_clone() const = 0;

			/// \brief TODOCUMENT
			virtual void do_output_superposition(const sup::superposition_context &,
			                                     std::ostream &,
			                                     const boost::string_ref &) const = 0;

			/// \brief TODOCUMENT
			virtual bool do_involves_display_spec() const = 0;

			/// \brief Pure virtual method with which each concrete superposition_outputter must define its name
			virtual std::string do_get_name() const = 0;

		public:
			superposition_outputter() = default;
			std::unique_ptr<superposition_outputter> clone() const;
			virtual ~superposition_outputter() noexcept = default;

			superposition_outputter(const superposition_outputter &) = default;
			superposition_outputter(superposition_outputter &&) noexcept = default;
			superposition_outputter & operator=(const superposition_outputter &) = default;
			superposition_outputter & operator=(superposition_outputter &&) noexcept = default;

			void output_superposition(const sup::superposition_context &,
			                          std::ostream &,
			                          const boost::string_ref &) const;
			bool involves_display_spec() const;

			std::string get_name() const;
		};

		/// \brief NVI pass-through to get the concrete superposition_outputter's name
		inline std::string superposition_outputter::get_name() const {
			return do_get_name();
		}

		/// \brief Function to make superposition_outputter meet the Clonable concept (used in ptr_container)
		///
		/// NOTE: Don't call this yourself. Call the object's clone() method instead because that returns a
		///       smart pointer rather than the raw pointer this has to return to meet the Clonable concept.
		///
		/// This gets the smart pointer from the clone() method and then calls release on it.
		///
		/// \returns A raw pointer to a new copy of the superposition_outputter argument, with the same dynamic type.
		///          The caller is responsible for deleting this new object.
		inline superposition_outputter * new_clone(const superposition_outputter &arg_superposition_outputter ///< The superposition_outputter to clone
		                                           ) {
			return arg_superposition_outputter.clone().release();
		}
	} // namespace opts
} // namespace cath

#endif
