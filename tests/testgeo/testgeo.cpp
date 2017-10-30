#include "testgeo.h"

GeoTest::GeoTest(){};

GeoTest::~GeoTest() {};

void GeoTest::SetUp() {};

void GeoTest::TearDown() {};


TEST_F(GeoTest, ByDefaultBazTrueIsTrue) {
    EXPECT_EQ(true, true);
}
