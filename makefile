
default: deltagtest.cpp HopscotchDelta.hpp HopscotchSegment.hpp
	g++ deltagtest.cpp -lgtest -lgtest_main -pthread -o test
