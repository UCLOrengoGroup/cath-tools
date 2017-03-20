/// \file
/// \brief The clone_ptr test suite

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
///
/// It might be helpful to run this test suite under Valgrind after any changes
/// to check for memory leaks.  Try a command something like the following:
/// \verbatim valgrind --leak-check=full debug/test --run_test=clone_ptr_test_suite \endverbatim
///
/// When doing that sort of Valgrind check, it also helps to temporarily comment the:
/// \verbatim BOOST_GLOBAL_FIXTURE(prepare_for_test_global_fixture) \endverbatim
/// line in test.cpp so that the only other error seen is the standard deregister_test_unit one.
///
/// There is no testcase specifically for get() but it is used throughout the testcases
/// so a failure should show up somewhere

#include <boost/test/auto_unit_test.hpp>

#include <boost/serialization/export.hpp>

#include "common/clone/check_uptr_clone_against_this.hpp"
#include "common/clone/clone_ptr.hpp"
#include "common/clone/make_uptr_clone.hpp"
#include "common/size_t_literal.hpp"

namespace cath { namespace test { } }

using namespace boost::archive;
using namespace boost::test_tools;
using namespace cath::common;
using namespace cath::common::detail;
using namespace cath::test;
using namespace std;

BOOST_TEST_DONT_PRINT_LOG_VALUE( nullptr_t )

const auto CONCRETE1_METHOD_RESULT = 3984756_z;
const auto CONCRETE2_METHOD_RESULT = 836_z;

namespace cath {
	namespace test {

		class clone_ptr_test_abstract_base {
		private:
			friend class boost::serialization::access;
			template<typename archive> void serialize(archive & /*ar*/,
			                                          const size_t /*version*/
			                                          ) {
			}
			virtual unique_ptr<clone_ptr_test_abstract_base> do_clone() const = 0;

		public:
			clone_ptr_test_abstract_base() = default;
			unique_ptr<clone_ptr_test_abstract_base> clone() const {
				return check_uptr_clone_against_this( do_clone(), *this );
			}
			virtual ~clone_ptr_test_abstract_base() noexcept = default;

			clone_ptr_test_abstract_base(const clone_ptr_test_abstract_base &) = default;
			clone_ptr_test_abstract_base(clone_ptr_test_abstract_base &&) noexcept = default;
			clone_ptr_test_abstract_base & operator=(const clone_ptr_test_abstract_base &) = default;
			clone_ptr_test_abstract_base & operator=(clone_ptr_test_abstract_base &&) noexcept = default;

			virtual size_t method() const = 0;
		};


		class clone_ptr_test_concrete1 final : public clone_ptr_test_abstract_base {
		private:
			friend class boost::serialization::access;
			template<typename archive> void serialize(archive & ar,
			                                          const size_t /*version*/
			                                          ) {
				ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(clone_ptr_test_abstract_base);
			}
			unique_ptr<clone_ptr_test_abstract_base> do_clone() const final {
				return { make_uptr_clone( *this ) };
			}

		public:
			size_t method() const final {
				return CONCRETE1_METHOD_RESULT;
			}
		};


		class clone_ptr_test_concrete2 final : public clone_ptr_test_abstract_base {
		private:
			friend class boost::serialization::access;
			template<typename archive> void serialize(archive & ar,
			                                       const size_t /*version*/
			                                       ) {
				ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(clone_ptr_test_abstract_base);
			}
			unique_ptr<clone_ptr_test_abstract_base> do_clone() const final {
				return { make_uptr_clone( *this ) };
			}

		public:
			size_t method() const final {
				return CONCRETE2_METHOD_RESULT;
			}
		};

	}  // namespace test
}  // namespace cath

//bool operator==(const clone_ptr_test_abstract_base &arg_obj1, ///< TODOCUMENT
//                const clone_ptr_test_abstract_base &arg_obj2  ///< TODOCUMENT
//                ) {
//	if (typeid(arg_obj1) != typeid(arg_obj2)) {
//		return false;
//	}
//	return (arg_obj1.method() == arg_obj2.method());
//}

BOOST_CLASS_EXPORT(clone_ptr_test_concrete1)
BOOST_CLASS_EXPORT(clone_ptr_test_concrete2)

BOOST_AUTO_TEST_SUITE(clone_ptr_test_suite)

//class has_clone_that_returns_int {
//public:
//	int clone() const {return 0;}
//};
//BOOST_AUTO_TEST_CASE(commented_out_test_for_compiler_errors) {
//	const clone_ptr<has_clone_that_returns_int> fail_ptr1;
//}

BOOST_AUTO_TEST_CASE(concrete_methods_return_different_values) {
	clone_ptr_test_concrete1 the_clone_ptr_test_concrete1;
	clone_ptr_test_concrete2 the_clone_ptr_test_concrete2;
	BOOST_CHECK_NE(the_clone_ptr_test_concrete1.method(), the_clone_ptr_test_concrete2.method());
}

BOOST_AUTO_TEST_CASE(default_constructs) {
	const clone_ptr<const clone_ptr_test_abstract_base> const_ptr;
	BOOST_CHECK( ! const_ptr );
	BOOST_CHECK_EQUAL( const_ptr.get(), nullptr );
}

BOOST_AUTO_TEST_CASE(constructs_from_raw_pointer) {
	clone_ptr_test_abstract_base * raw_ptr = new clone_ptr_test_concrete1;
	clone_ptr<const clone_ptr_test_abstract_base> ptr(raw_ptr);
	BOOST_CHECK( ptr );
	BOOST_CHECK_EQUAL(ptr.get(), raw_ptr);
}

BOOST_AUTO_TEST_CASE(constructs_from_lvalue_ref_to_unique_ptr) {
	clone_ptr_test_abstract_base * raw_ptr = new clone_ptr_test_concrete1;
	auto       the_unique_ptr = unique_ptr<clone_ptr_test_abstract_base>{ raw_ptr };
	const auto ptr            = clone_ptr<clone_ptr_test_abstract_base>{ the_unique_ptr };
	BOOST_CHECK      ( the_unique_ptr );
	BOOST_CHECK      ( ptr            );
	BOOST_CHECK_EQUAL( the_unique_ptr.get(), raw_ptr   );
	BOOST_CHECK_NE   ( the_unique_ptr.get(), ptr.get() );
	BOOST_CHECK_EQUAL( the_unique_ptr->method(), CONCRETE1_METHOD_RESULT );
	BOOST_CHECK_EQUAL( ptr->method(),            CONCRETE1_METHOD_RESULT );
}

BOOST_AUTO_TEST_CASE(constructs_from_rvalue_ref_to_unique_ptr) {
	clone_ptr_test_abstract_base * raw_ptr = new clone_ptr_test_concrete1;
	const auto ptr = clone_ptr<clone_ptr_test_abstract_base>{ unique_ptr<clone_ptr_test_abstract_base>{ raw_ptr } };
	BOOST_CHECK      ( ptr                  );
	BOOST_CHECK_EQUAL( ptr.get(), raw_ptr   );
	BOOST_CHECK_EQUAL( ptr->method(), CONCRETE1_METHOD_RESULT );
}

BOOST_AUTO_TEST_CASE(constructs_from_lvalue_ref_to_unique_ptr_to_const) {
	const clone_ptr_test_abstract_base * raw_ptr = new clone_ptr_test_concrete1;
	auto       the_unique_ptr = unique_ptr<const clone_ptr_test_abstract_base>{ raw_ptr };
	const auto ptr            = clone_ptr<const clone_ptr_test_abstract_base>{ the_unique_ptr };
	BOOST_CHECK      ( the_unique_ptr );
	BOOST_CHECK      ( ptr            );
	BOOST_CHECK_EQUAL( the_unique_ptr.get(), raw_ptr   );
	BOOST_CHECK_NE   ( the_unique_ptr.get(), ptr.get() );
	BOOST_CHECK_EQUAL( the_unique_ptr->method(), CONCRETE1_METHOD_RESULT );
	BOOST_CHECK_EQUAL( ptr->method(),            CONCRETE1_METHOD_RESULT );
}

BOOST_AUTO_TEST_CASE(constructs_from_rvalue_ref_to_unique_ptr_to_const) {
	const clone_ptr_test_abstract_base * raw_ptr = new clone_ptr_test_concrete1;
	const auto ptr = clone_ptr<const clone_ptr_test_abstract_base>{ unique_ptr<const clone_ptr_test_abstract_base>{ raw_ptr } };
	BOOST_CHECK      ( ptr                  );
	BOOST_CHECK_EQUAL( ptr.get(), raw_ptr   );
	BOOST_CHECK_EQUAL( ptr->method(), CONCRETE1_METHOD_RESULT );
}

/// \todo Should separate this into two tests: move-constructs (ie from rvalue) and copy-constructs
BOOST_AUTO_TEST_CASE(copy_constructs) {
	clone_ptr_test_abstract_base * raw_ptr = new clone_ptr_test_concrete1;
	const clone_ptr<clone_ptr_test_abstract_base> source_ptr( raw_ptr    );
	const clone_ptr<clone_ptr_test_abstract_base> dest_ptr  ( source_ptr );
	BOOST_CHECK( source_ptr );
	BOOST_CHECK( dest_ptr );
	BOOST_CHECK_EQUAL( source_ptr.get(), raw_ptr );
	BOOST_CHECK_NE   ( dest_ptr.get(),   raw_ptr );
	BOOST_CHECK_EQUAL( source_ptr->method(), CONCRETE1_METHOD_RESULT );
	BOOST_CHECK_EQUAL( dest_ptr->method(),   CONCRETE1_METHOD_RESULT );
}

/// \todo Should separate this into two tests: move-constructs (ie from rvalue) and copy-constructs
///
/// \todo Consider reducing redundancy by making template testcases for everything that should be repeated for:
///        * clone_ptr<      clone_ptr_test_abstract_base> and
///        * clone_ptr<const clone_ptr_test_abstract_base>
BOOST_AUTO_TEST_CASE(copy_constructs_with_const) {
	const clone_ptr_test_abstract_base * raw_ptr = new clone_ptr_test_concrete1;
	const clone_ptr<const clone_ptr_test_abstract_base> source_ptr(raw_ptr);
	const clone_ptr<const clone_ptr_test_abstract_base> dest_ptr(source_ptr);
	BOOST_CHECK( source_ptr );
	BOOST_CHECK( dest_ptr   );
	BOOST_CHECK_EQUAL( source_ptr.get(), raw_ptr);
	BOOST_CHECK_NE   ( dest_ptr.get(),   raw_ptr);
	BOOST_CHECK_EQUAL( source_ptr->method(), CONCRETE1_METHOD_RESULT );
	BOOST_CHECK_EQUAL( dest_ptr->method(),   CONCRETE1_METHOD_RESULT );
}

/// \todo Should separate this into two tests: move-assignment (ie from rvalue) and copy-assignment
BOOST_AUTO_TEST_CASE(assignment_operator_works) {
	clone_ptr_test_abstract_base * raw_ptr = new clone_ptr_test_concrete1;
	clone_ptr<clone_ptr_test_abstract_base> source_ptr(raw_ptr);

	clone_ptr<clone_ptr_test_abstract_base> previously_null_dest_ptr;
	clone_ptr<clone_ptr_test_abstract_base> previously_concrete2_dest_ptr(new clone_ptr_test_concrete2);

	BOOST_CHECK_EQUAL( previously_null_dest_ptr.get(), nullptr );
	BOOST_CHECK_EQUAL( previously_concrete2_dest_ptr->method(), CONCRETE2_METHOD_RESULT );

	previously_null_dest_ptr = source_ptr;
	previously_concrete2_dest_ptr = source_ptr;

	BOOST_CHECK_EQUAL( source_ptr.get(), raw_ptr );
	BOOST_CHECK_NE( previously_null_dest_ptr.get(),      raw_ptr );
	BOOST_CHECK_NE( previously_concrete2_dest_ptr.get(), raw_ptr );
	BOOST_CHECK_NE( previously_concrete2_dest_ptr.get(), previously_null_dest_ptr.get() );
	BOOST_CHECK_EQUAL( previously_null_dest_ptr->method(),      CONCRETE1_METHOD_RESULT );
	BOOST_CHECK_EQUAL( previously_concrete2_dest_ptr->method(), CONCRETE1_METHOD_RESULT );
}

BOOST_AUTO_TEST_CASE(lvalue_ref_to_unique_ptr_assignment_operator) {
	clone_ptr_test_abstract_base * raw_ptr = new clone_ptr_test_concrete1;
	auto ptr            = clone_ptr<clone_ptr_test_abstract_base>{};
	auto the_unique_ptr = unique_ptr<clone_ptr_test_abstract_base>{ raw_ptr };
	ptr                 = the_unique_ptr;
	BOOST_CHECK      ( the_unique_ptr );
	BOOST_CHECK      ( ptr            );
	BOOST_CHECK_EQUAL( the_unique_ptr.get(), raw_ptr   );
	BOOST_CHECK_NE   ( the_unique_ptr.get(), ptr.get() );
	BOOST_CHECK_EQUAL( the_unique_ptr->method(), CONCRETE1_METHOD_RESULT );
	BOOST_CHECK_EQUAL( ptr->method(),            CONCRETE1_METHOD_RESULT );
}

BOOST_AUTO_TEST_CASE(rvalue_ref_to_unique_ptr_assignment_operator) {
	clone_ptr_test_abstract_base * raw_ptr = new clone_ptr_test_concrete1;
	auto ptr = clone_ptr<clone_ptr_test_abstract_base>{};
	ptr      = unique_ptr<clone_ptr_test_abstract_base>{ raw_ptr };
	BOOST_CHECK      ( ptr                  );
	BOOST_CHECK_EQUAL( ptr.get(), raw_ptr   );
	BOOST_CHECK_EQUAL( ptr->method(), CONCRETE1_METHOD_RESULT );
}

BOOST_AUTO_TEST_CASE(pass_through_reset) {
	clone_ptr_test_abstract_base * raw_ptr = new clone_ptr_test_concrete1;
	auto ptr = clone_ptr<clone_ptr_test_abstract_base>{ raw_ptr };
	ptr.reset( new clone_ptr_test_concrete2 );
	BOOST_CHECK_NE(ptr.get(), raw_ptr);
	BOOST_CHECK_EQUAL(ptr->method(), CONCRETE2_METHOD_RESULT);
	ptr.reset( nullptr );
	BOOST_CHECK_EQUAL(ptr.get(), nullptr );
}

BOOST_AUTO_TEST_CASE(pass_through_dereference_operator) {
	clone_ptr<clone_ptr_test_abstract_base> ptr(new clone_ptr_test_concrete1);
	clone_ptr_test_abstract_base & ref = *ptr;
	BOOST_CHECK_EQUAL(ref.method(), CONCRETE1_METHOD_RESULT);
}

BOOST_AUTO_TEST_CASE(pass_through_member_by_pointer_operator) {
	clone_ptr<clone_ptr_test_abstract_base> ptr(new clone_ptr_test_concrete1);
	BOOST_CHECK_EQUAL(ptr->method(), CONCRETE1_METHOD_RESULT);
}

BOOST_AUTO_TEST_CASE(pass_through_swap) {
	clone_ptr<clone_ptr_test_abstract_base> ptr1(new clone_ptr_test_concrete1);
	clone_ptr<clone_ptr_test_abstract_base> ptr2(new clone_ptr_test_concrete2);
	BOOST_CHECK_EQUAL(ptr1->method(), CONCRETE1_METHOD_RESULT);
	BOOST_CHECK_EQUAL(ptr2->method(), CONCRETE2_METHOD_RESULT);
	ptr1.swap( ptr2 );
	BOOST_CHECK_EQUAL(ptr1->method(), CONCRETE2_METHOD_RESULT);
	BOOST_CHECK_EQUAL(ptr2->method(), CONCRETE1_METHOD_RESULT);
}

BOOST_AUTO_TEST_SUITE_END()
