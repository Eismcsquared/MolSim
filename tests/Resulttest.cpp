#include <gtest/gtest.h>

TEST(Test3, Test3) {
    std::cout << "Test1" << std::endl;
    ASSERT_TRUE(true);
}

TEST(Test4, Test4) {
    std::cout << "Test2" << std::endl;
    EXPECT_FALSE(false);
}