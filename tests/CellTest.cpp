#include <vector>
#include <gtest/gtest.h>
#include <spdlog/spdlog.h>
#include "container/Cell.h"
#include "Logger.h"

class CellTest: public ::testing::Test {
protected:
    Cell cell = Cell(std::array<double, 3>{5, 10, 0}, std::array<double, 3>{3, 1, 1});
};

// Test whether adding and removing indices from a cell works correctly.
TEST_F(CellTest, AddAndRemove) {
    test_logger->info("Cell - Add and remove test");
    for(unsigned int i = 0; i < 10; ++i) {
        cell.addIndex(i);
    }
    cell.removeIndex(0);
    cell.removeIndex(9);
    cell.removeIndex(3);
    cell.removeIndex(7);
    // remove indicis that are not in the cell
    cell.removeIndex(-1);
    cell.removeIndex(10);
    EXPECT_EQ(6, cell.getParticleIndices().size());
    std::set<unsigned int> ref = {1, 2, 4, 5, 6, 8};
    std::set<unsigned int> indexSet(cell.getParticleIndices().begin(), cell.getParticleIndices().end());
    EXPECT_EQ(ref, indexSet);
    if (::testing::Test::HasFailure()) {
        test_logger->info("Cell - Add and remove test failed");
    } else {
        test_logger->info("Cell - Add and remove test passed");
    }
}

// Test whether the contains function correctly determines whether a position is in the cell
TEST_F(CellTest, Contains) {
    test_logger->info("Cell - contains test");
    // In the cell
    EXPECT_TRUE(cell.contains(std::array<double, 3>{6, 10, 0.5}));
    // On boundary
    EXPECT_TRUE(cell.contains(std::array<double, 3>{5, 10.5, 0.3}));
    EXPECT_TRUE(cell.contains(std::array<double, 3>{5, 10, 0.4}));
    EXPECT_TRUE(cell.contains(std::array<double, 3>{8, 11, 0}));
    // ourside
    EXPECT_FALSE(cell.contains(std::array<double, 3>{4, 10.4, 0.8}));
    EXPECT_FALSE(cell.contains(std::array<double, 3>{11, 10.4, 0.8}));
    EXPECT_FALSE(cell.contains(std::array<double, 3>{6, 8.5, 0.8}));
    EXPECT_FALSE(cell.contains(std::array<double, 3>{6, 11.5, 0.8}));
    EXPECT_FALSE(cell.contains(std::array<double, 3>{6, 10.5, -1}));
    EXPECT_FALSE(cell.contains(std::array<double, 3>{6, 10.5, 3}));
    if (::testing::Test::HasFailure()) {
        test_logger->info("Cell - contains test failed");
    } else {
        test_logger->info("Cell - contains test passed");
    }
}

