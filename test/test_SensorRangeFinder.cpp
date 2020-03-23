/****************************************************************************
 *
 *   Copyright (c) 2020 ECL Development Team. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the
 *    distribution.
 * 3. Neither the name PX4 nor the names of its contributors may be
 *    used to endorse or promote products derived from this software
 *    without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 * FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 * COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
 * OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 * AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 * ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 ****************************************************************************/

#include <gtest/gtest.h>
#include <math.h>
#include "EKF/common.h"
#include "EKF/SensorRangeFinder.hpp"
#include <matrix/math.hpp>

using estimator::rangeSample;
using matrix::Dcmf;
using namespace estimator::sensor;

class SensorRangeFinderTest : public ::testing::Test {
public:
	// Setup the Ekf with synthetic measurements
	void SetUp() override
	{
	}

	// Use this method to clean up any memory, network etc. after each test
	void TearDown() override
	{
	}

protected:
	SensorRangeFinder _range_finder{};
};


TEST_F(SensorRangeFinderTest, setData)
{
	rangeSample sample{};
	sample.rng = 1.f;
	sample.quality = 9.f;
	sample.time_us = 1e6f;
	Dcmf attitude{};

	_range_finder.setDelayedRng(sample.rng);
	_range_finder.setValidity(true);
	/* _range_finder_data.runChecks(1e6f, attitude); */
	EXPECT_TRUE(_range_finder.isHealthy());
}
