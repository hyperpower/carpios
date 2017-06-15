#ifndef __TEST_ANY_HPP_
#define __TEST_ANY_HPP_

#include "gtest/gtest.h"
#include "../../lib/utility/any.hpp"

#include <math.h>
#include <string>
#include <memory>

namespace carpio {
template<size_t N>
struct words {
	void* w[N];
};

struct big_type {
	char i_wanna_be_big[256];
	std::string value;

	big_type() :
			value(std::string(300, 'b')) {
		i_wanna_be_big[0] = i_wanna_be_big[50] = 'k';
	}

	bool check() {
		EXPECT_EQ(value.size(), 300);
		//ASSERT_TRUE(value.front() == 'b' && value.back() == 'b');
		//ASSERT_TRUE(i_wanna_be_big[0] == 'k' && i_wanna_be_big[50] == 'k');
		return true;
	}
};

// small type which has nothrow move ctor but throw copy ctor
struct regression1_type {
	const void* confuse_stack_storage = (void*) (0);
	regression1_type() {
	}
	regression1_type(const regression1_type&) {
	}
	regression1_type(regression1_type&&) noexcept {
	}
	regression1_type& operator=(const regression1_type&) {
		return *this;
	}
	regression1_type& operator=(regression1_type&&) {
		return *this;
	}
};

TEST(Any, empty) {
	Any x = 4;
	Any y = big_type();
	Any z = 6;

	EXPECT_TRUE(Any().empty());
	EXPECT_TRUE(!Any(1).empty());
	EXPECT_TRUE(!Any(big_type()).empty());

	EXPECT_TRUE(!x.empty() && !y.empty() && !z.empty());
	y.clear();
	EXPECT_TRUE(!x.empty() && y.empty() && !z.empty());
	x = y;
	EXPECT_TRUE(x.empty() && y.empty() && !z.empty());
	z = Any();
	EXPECT_TRUE(x.empty() && y.empty() && z.empty());
}

TEST(Any, type) {
	EXPECT_TRUE(Any().type() == typeid(void));
	EXPECT_TRUE(Any(4).type() == typeid(int));
	EXPECT_TRUE(Any(big_type()).type() == typeid(big_type));
	EXPECT_TRUE(Any(1.5f).type() == typeid(float));
}

TEST(Any, except) {
	bool except0 = false;
	bool except1 = false, except2 = false;
	bool except3 = false, except4 = false;

	try {
		any_cast<int>(Any());
	} catch (const bad_any_cast&) {
		except0 = true;
	}

	try {
		any_cast<int>(Any(4.0f));
	} catch (const bad_any_cast&) {
		except1 = true;
	}

	try {
		any_cast<float>(Any(4.0f));
	} catch (const bad_any_cast&) {
		except2 = true;
	}

	try {
		any_cast<float>(Any(big_type()));
	} catch (const bad_any_cast&) {
		except3 = true;
	}

	try {
		any_cast<big_type>(Any(big_type()));
	} catch (const bad_any_cast&) {
		except4 = true;
	}

	EXPECT_TRUE(except0 == true);
	EXPECT_TRUE(except1 == true && except2 == false);
	EXPECT_TRUE(except3 == true && except4 == false);
}

TEST(Any, cast) {
	Any i4 = 4;
	Any i5 = 5;
	Any f6 = 6.0f;
	Any big1 = big_type();
	Any big2 = big_type();
	Any big3 = big_type();

	EXPECT_TRUE(any_cast<int>(&i4) != nullptr);
	EXPECT_TRUE(any_cast<float>(&i4) == nullptr);
	EXPECT_TRUE(any_cast<int>(i5) == 5);
	EXPECT_TRUE(any_cast<float>(f6) == 6.0f);
	EXPECT_TRUE(
			any_cast<big_type>(big1).check() && any_cast<big_type>(big2).check()
					&& any_cast<big_type>(big3).check());
}

TEST(Any, shared_pointer) {
	std::shared_ptr<int> ptr_count(new int);
	std::weak_ptr<int> weak = ptr_count;
	Any p0 = 0;

	EXPECT_TRUE(weak.use_count() == 1);
	Any p1 = ptr_count;
	EXPECT_TRUE(weak.use_count() == 2);
	Any p2 = p1;
	EXPECT_TRUE(weak.use_count() == 3);
	p0 = p1;
	EXPECT_TRUE(weak.use_count() == 4);
	p0 = 0;
	EXPECT_TRUE(weak.use_count() == 3);
	p0 = std::move(p1);
	EXPECT_TRUE(weak.use_count() == 3);
	p0.swap(p1);
	EXPECT_TRUE(weak.use_count() == 3);
	p0 = 0;
	EXPECT_TRUE(weak.use_count() == 3);
	p1.clear();
	EXPECT_TRUE(weak.use_count() == 2);
	p2 = Any(big_type());
	EXPECT_TRUE(weak.use_count() == 1);
	p1 = ptr_count;
	EXPECT_TRUE(weak.use_count() == 2);
	ptr_count = nullptr;
	EXPECT_TRUE(weak.use_count() == 1);
	p1 = Any();
	EXPECT_TRUE(weak.use_count() == 0);
}

TEST(Any, stack) {
	auto is_stack_allocated = [](const Any& a, const void* obj1) {
		uintptr_t a_ptr = (uintptr_t)(&a);
		uintptr_t obj = (uintptr_t)(obj1);
		return (obj >= a_ptr && obj < a_ptr + sizeof(Any));
	};

	//static_assert(sizeof(std::unique_ptr<big_type>) <= sizeof(void*) * 1, "unique_ptr too big");
	static_assert(sizeof(std::shared_ptr<big_type>) <= sizeof(void*) * 2, "shared_ptr too big");

	Any i = 400;
	Any f = 400.0f;
	//any unique = std::unique_ptr<big_type>(); -- must be copy constructible
	Any shared = std::shared_ptr<big_type>();
	Any rawptr = (void*) (nullptr);
	Any big = big_type();
	Any w2 = words<2>();
	Any w3 = words<3>();

	EXPECT_TRUE(is_stack_allocated(i, any_cast<int>(&i)));
	EXPECT_TRUE(is_stack_allocated(f, any_cast<float>(&f)));
	EXPECT_TRUE(is_stack_allocated(rawptr, any_cast<void*>(&rawptr)));
	//EXPECT_TRUE(is_stack_allocated(unique, any_cast<std::unique_ptr<big_type>>(&unique)));
	EXPECT_TRUE(
			is_stack_allocated(shared,
					any_cast<std::shared_ptr<big_type>>(&shared)));
	EXPECT_TRUE(!is_stack_allocated(big, any_cast<big_type>(&big)));
	EXPECT_TRUE(is_stack_allocated(w2, any_cast<words<2>>(&w2)));
	EXPECT_TRUE(!is_stack_allocated(w3, any_cast<words<3>>(&w3)));

	// Regression test for GitHub Issue #1
	Any r1 = regression1_type();
	EXPECT_TRUE(is_stack_allocated(r1, any_cast<const regression1_type>(&r1)));
}

}

#endif
