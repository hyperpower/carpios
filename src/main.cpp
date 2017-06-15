//============================================================================
// Name        : carpio.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include "../test/_test.hpp"


using namespace std;

int main(int argc, char **argv) {
	int res = carpio::RunTests(argc, argv);
    std::cout<<" ====  end ====";
    return res;
}
