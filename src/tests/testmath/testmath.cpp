#include "testmath.h"
#include "mathlib/math.h"

MathTest::MathTest(){};

MathTest::~MathTest() {};

void MathTest::SetUp() {};

void MathTest::TearDown() {};


TEST_F(MathTest, Max) {
	float v = math::max(1.0f, 3.0f);
    EXPECT_GT(v, 1.0f);
    EXPECT_GT(v, 3.0f);
}
