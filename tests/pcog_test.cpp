//
// Created by rolf on 18-11-22.
//

#include <pcog/pcog.hpp>
#include <gtest/gtest.h>

TEST(add_test,add_1_2){
   EXPECT_EQ(pcog::add(1,2),3);
}
