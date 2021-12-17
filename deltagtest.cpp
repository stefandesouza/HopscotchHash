#include "HopscotchSegment.hpp"
#include <stdlib.h>
#include <future>
#include <vector>
#include <gtest/gtest.h>
using namespace GaspiLS;
const int GLOBAL_SIZE = 1000;
const int GLOBAL_HOP_RANGE = 16;
const double GLOBAL_CAPACITY = 0.80;

/*
TEST(SegmentHopscotchTest, insertTestSeq)
{   
    int size = GLOBAL_SIZE;
    int hop_range = GLOBAL_HOP_RANGE;
    GaspiLS::HopscotchSegment<int, int> table(size, hop_range);
    std::vector<int> elements;
    int numberOfElements = (int)(GLOBAL_CAPACITY * size);

    for (int i = 0; i < numberOfElements; ++i)
    {   
        elements.push_back(rand());
        EXPECT_EQ(true, table.insert(elements[i], elements[i]).second);
    }

    for (int i = 0; i < numberOfElements; ++i)
    {   
        EXPECT_EQ(true, table.contains(elements[i]));
    }
}

TEST(SegmentHopscotchTest, removeTestSeq)
{   
    int size = GLOBAL_SIZE;
    int hop_range = GLOBAL_HOP_RANGE;
    GaspiLS::HopscotchSegment<int, int> table(size, hop_range);
    std::vector<int> elements;
    int numberOfElements = (int)(GLOBAL_CAPACITY * size);

    for (int i = 0; i < numberOfElements; ++i)
    {   
        elements.push_back(rand());
        table.insert(elements[i], elements[i]);
    }

    for (int i = 0; i < numberOfElements; ++i)
    {   
        EXPECT_EQ(true, table.remove(elements[i]));
    }

    for (int i = 0; i < numberOfElements; ++i)
    {   
        EXPECT_EQ(false, table.contains(elements[i]));
    }

}
*/

TEST(SegmentHopscotchTest, insertTest)
{   
    int size = GLOBAL_SIZE;
    int hop_range = GLOBAL_HOP_RANGE;
    GaspiLS::HopscotchSegment<int, int> table(size, hop_range);
    std::vector<int> elements;
    int numberOfElements = (int)(GLOBAL_CAPACITY * size);
    for (int i = 0; i <numberOfElements; ++i)
    {
        elements.push_back(rand());
    }

    int numthreads = 2;
    std::future<std::pair<int,bool>> thread[numthreads];
    std::future<bool> checks[numthreads];
    std::pair<int,bool> end[numthreads];
    bool endCheck[numberOfElements];

    for (int j = 0; j < numberOfElements / numthreads; ++j)
    {
        for (int i = 0; i < numthreads; ++i)
        {   
            thread[i] = std::async(std::launch::async, &HopscotchSegment<int,int>::insert, &table, elements[i + j * numthreads], elements[i + j * numthreads]);
        }
        for (int i = 0; i < numthreads; ++i)
        {
            end[i] = thread[i].get();
            EXPECT_EQ(true, end[i].second);
        }
    }

    for (int j = 0; j < numberOfElements / numthreads; ++j)
    {
        for (int i = 0; i < numthreads; ++i)
        {   
            checks[i] = std::async(std::launch::async, &HopscotchSegment<int,int>::contains, &table, elements[i + j * numthreads]);
        }
        for (int i = 0; i < numthreads; ++i)
        {
            endCheck[i] = checks[i].get();
            EXPECT_EQ(true, endCheck[i]);
        }
    }
}


TEST(SegmentHopscotchTest, removeTest)
{   
    int size = GLOBAL_SIZE;
    int hop_range = GLOBAL_HOP_RANGE;
    GaspiLS::HopscotchSegment<int, int> table(size, hop_range);
    std::vector<int> elements;
    int numberOfElements = (int)(GLOBAL_CAPACITY * size);
    for (int i = 0; i <numberOfElements; ++i)
    {
        elements.push_back(rand());
    }

    int numthreads = 2;
    std::future<std::pair<int,bool>> thread[numthreads];
    std::future<bool> removeThread[numthreads];
    std::future<bool> checks[numthreads];
    std::pair<int,bool> end[numthreads];
    bool endCheck[numberOfElements];
    
    // add half the number of elements
    for (int j = 0; j < numberOfElements / numthreads; ++j)
    {
        for (int i = 0; i < numthreads / 2; ++i)
        {   
            thread[i] = std::async(std::launch::async, &HopscotchSegment<int,int>::insert, &table, elements[i + j * numthreads / 2], elements[i + j * numthreads / 2]);
            end[i] = thread[i].get();
            EXPECT_EQ(true, end[i].second);
        }
    }

    // remove the previous elements, add the next half
    for (int j = 0; j < numberOfElements /  numthreads; ++j)
    {
        for (int i = 0; i < numthreads / 2; ++i)
        {   
            removeThread[i] = std::async(std::launch::async, &HopscotchSegment<int,int>::remove, &table, elements[i + j * numthreads / 2]);

            thread[numthreads / 2 + i] = std::async(std::launch::async, &HopscotchSegment<int,int>::insert, &table, elements[numberOfElements / 2 + i + j * numthreads / 2], elements[numberOfElements / 2 + i + j * numthreads / 2]);
        }

        for (int i = 0; i < numthreads / 2; ++i)
        {
            endCheck[i] = removeThread[i].get();
            EXPECT_EQ(true, endCheck[i]);
            end[numthreads / 2 + i] = thread[numthreads / 2 + i].get();
            EXPECT_EQ(true, end[numthreads / 2 + i].second);
        }

        for (int i = 0; i < numthreads / 2; ++i)
        {   
            checks[i] = std::async(std::launch::async, &HopscotchSegment<int,int>::contains, &table, elements[i + j * numthreads / 2]);
            checks[numthreads / 2 + i] = std::async(std::launch::async, &HopscotchSegment<int,int>::contains, &table, elements[numberOfElements / 2 + i + j * numthreads / 2]);
        }
        for (int i = 0; i < numthreads / 2; ++i)
        {
            endCheck[i] = checks[i].get();
            EXPECT_EQ(false, endCheck[i]);
            endCheck[numthreads / 2 + i] = checks[numthreads / 2 + i].get();
            EXPECT_EQ(true, endCheck[numthreads / 2 + i]);
        }
    }
}


int main(int argc, char *argv[])
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}