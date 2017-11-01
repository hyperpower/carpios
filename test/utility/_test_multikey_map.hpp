/*
 * _test_multikey_map.hpp
 *
 *  Created on: Oct 3, 2017
 *      Author: zhou
 */

#ifndef _TEST_MULTIKEY_MAP_HPP_
#define _TEST_MULTIKEY_MAP_HPP_

#include "gtest/gtest.h"
#include "../../lib/utility/multikey_map.hpp"
#include <iostream>
#include <assert.h>

namespace carpio {

typedef MultikeyMap<int, int, int> MMiii;

void print_map(MMiii& mmap) {
	for (auto& p : mmap) {
		std::cout << p.second->key1 << " " << p.second->key2 << " "
				<< p.second->val << std::endl;
	}
	std::cout << "-----\n";
}

TEST(multikey, insert) {
	MMiii mmap = { { 1, 2, 3 }, { 1, 2, 4 }, { 4, 5, 6 }, { 4, 5, 7 },
			{ 7, 8, 9 } };

	mmap.insert(1, 12, 13);
	mmap.insert(14, 5, 16);
	mmap.insert(17, 18, 19);
	mmap.insert( { 21, 22, 23 });

	std::cout << "count1 = 1 --> ";
	std::cout << mmap.count1(1) << std::endl;
	std::cout << "count1 = 2 --> ";
	std::cout << mmap.count1(2) << std::endl;
	std::cout << "count2 = 1 --> ";
	std::cout << mmap.count2(1) << std::endl;
	std::cout << "count2 = 5 --> ";
	std::cout << mmap.count2(5) << std::endl;

	ASSERT_EQ(mmap.count1(1), 3);

	print_map(mmap);

	std::cout << "insert same \n";
	mmap.insert(1, 12, 13);
	print_map(mmap);
}

TEST(multikey, erase) {
	MMiii mmap = { { 1, 2, 3 }, { 1, 2, 4 }, { 4, 5, 6 }, { 4, 5, 7 },
			{ 7, 8, 9 } };

	mmap.insert(1, 12, 13);
	mmap.insert(14, 5, 16);
	mmap.insert(17, 18, 19);
	mmap.insert( { 21, 22, 23 });
	std::cout << "The map :\n";
	print_map(mmap);

	std::cout << "size = "<< mmap.size() << std::endl;
	std::cout << "erase key1 =4 \n";
	mmap.erase1(4);
	std::cout << "size = "<< mmap.size() << std::endl;
	std::cout << "erase key2 =8\n";
	mmap.erase2(8);
	std::cout << "size = "<< mmap.size() << std::endl;

	std::cout << "erase key1 =17 key2 = 18\n";
	mmap.erase(17, 18);
	std::cout << "size = "<< mmap.size() << std::endl;

	print_map(mmap);
}

TEST(multikey, get) {
	MMiii mmap = { { 1, 2, 3 }, { 1, 2, 4 }, { 4, 5, 6 }, { 4, 5, 7 },
			{ 7, 8, 9 } };

	mmap.insert(1, 12, 13);
	mmap.insert(14, 5, 16);
	mmap.insert(17, 18, 19);
	mmap.insert( { 21, 22, 23 });

	std::vector<MMiii::EntryPtr> vec1 = mmap.get1(1);
	std::vector<MMiii::EntryPtr> vec2 = mmap.get2(5);
	std::cout << vec1.size() << " " << vec2.size() << std::endl;
}

}

#endif /* TEST_ALGEBRA__TEST_MULTIKEY_MAP_HPP_ */
