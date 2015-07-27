/// \file
/// \brief The alignment_outputter class header

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

#ifndef ALIGNMENT_OUTPUTTER_H_INCLUDED
#define ALIGNMENT_OUTPUTTER_H_INCLUDED

#include <iosfwd>
#include <memory>

namespace cath { namespace align { class alignment_context; } }

namespace cath {
	namespace opts {

		/// \brief TODOCUMENT
		class alignment_outputter {
		private:
			/// \brief Pure virtual method with which each concrete alignment_outputter must define how to create a clone of itself
			virtual std::unique_ptr<alignment_outputter> do_clone() const = 0;
			
			/// \brief TODOCUMENT
			virtual void do_output_alignment(const align::alignment_context &,
			                                 std::ostream &) const = 0;
			/// \brief TODOCUMENT
			virtual bool do_involves_display_spec() const = 0;

		public:
			std::unique_ptr<alignment_outputter> clone() const;
			virtual ~alignment_outputter() noexcept = default;

			void output_alignment(const align::alignment_context &,
			                      std::ostream &) const;
			bool involves_display_spec() const;
		};

		/// \brief Function to make alignment_outputter meet the Clonable concept (used in ptr_container)
		///
		/// NOTE: Don't call this yourself. Call the object's clone() method instead because that returns a
		///       smart pointer rather than the raw pointer this has to return to meet the Clonable concept.
		///
		/// This gets the smart pointer from the clone() method and then calls release on it.
		///
		/// \returns A raw pointer to a new copy of the alignment_outputter argument, with the same dynamic type.
		///          The caller is responsible for deleting this new object.
		inline alignment_outputter * new_clone(const alignment_outputter &arg_alignment_outputter ///< The alignment_outputter to clone
		                                       ) {
			return arg_alignment_outputter.clone().release();
		}
	}
}

#endif
