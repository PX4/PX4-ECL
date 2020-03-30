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
using matrix::Eulerf;
using namespace estimator::sensor;

class SensorRangeFinderTest : public ::testing::Test {
public:
	// Setup the Ekf with synthetic measurements
	void SetUp() override
	{
		_range_finder.setSensorTilt(0.f);
		_range_finder.setCosMaxTilt(0.707f);
		_range_finder.setLimits(_min_range, _max_range);
	}

	// Use this method to clean up any memory, network etc. after each test
	void TearDown() override
	{
	}

protected:
	SensorRangeFinder _range_finder{};
	const rangeSample _good_sample{1.f, (uint64_t)1e6, 100}; // {range, time_us, quality}
	const float _min_range{0.5f};
	const float _max_range{10.f};
};


TEST_F(SensorRangeFinderTest, setRange)
{
	rangeSample sample{};
	sample.rng = 1.f;
	sample.time_us = 1e6;
	sample.quality = 9;

	_range_finder.setDelayedRng(sample.rng);
	_range_finder.setValidity(true);
	EXPECT_TRUE(_range_finder.isHealthy());
}

TEST_F(SensorRangeFinderTest, goodData)
{
	// WHEN: the drone is leveled and the data is good
	Dcmf attitude{Eulerf(0.f, 0.f, 0.f)};
	_range_finder.setDelayedSample(_good_sample);
	_range_finder.runChecks(1e6f, attitude);

	// THEN: the data can be used for aiding
	EXPECT_TRUE(_range_finder.isDelayedDataHealthy());
}

TEST_F(SensorRangeFinderTest, tiltExceeded)
{
	// WHEN: the drone is excessively tilted
	Dcmf attitude{Eulerf(0.f, 1.f, 0.f)};

	_range_finder.setDelayedSample(_good_sample);
	_range_finder.runChecks(1e6f, attitude);

	// THEN: the data should be marked as unhealthy
	EXPECT_FALSE(_range_finder.isDelayedDataHealthy());
}

TEST_F(SensorRangeFinderTest, rangeExceeded)
{
	// WHEN: the measured range is larger than the maximum
	Dcmf attitude{Eulerf(0.f, 0.f, 0.f)};

	rangeSample bad_sample = _good_sample;
	bad_sample. rng = 100.f;
	_range_finder.setDelayedSample(bad_sample);
	_range_finder.runChecks(1e6f, attitude);

	// THEN: the data should be marked as unhealthy
	EXPECT_FALSE(_range_finder.isDelayedDataHealthy());
}
