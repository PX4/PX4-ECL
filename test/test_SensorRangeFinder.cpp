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
	const rangeSample _good_sample{1.f, (uint64_t)2e6, 100}; // {range, time_us, quality}
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
	_range_finder.runChecks(_good_sample.time_us, attitude);

	// THEN: the data can be used for aiding
	EXPECT_TRUE(_range_finder.isDelayedDataHealthy());
}

TEST_F(SensorRangeFinderTest, tiltExceeded)
{
	// WHEN: the drone is excessively tilted
	Dcmf attitude{Eulerf(0.f, 1.f, 0.f)};

	_range_finder.setDelayedSample(_good_sample);
	_range_finder.runChecks(_good_sample.time_us, attitude);

	// THEN: the data should be marked as unhealthy
	EXPECT_FALSE(_range_finder.isDelayedDataHealthy());
}

TEST_F(SensorRangeFinderTest, rangeMaxExceeded)
{
	Dcmf attitude{Eulerf(0.f, 0.f, 0.f)};

	// WHEN: the measured range is larger than the maximum
	rangeSample bad_sample = _good_sample;
	bad_sample.rng = _max_range + 0.01f;
	_range_finder.setDelayedSample(bad_sample);
	_range_finder.runChecks(bad_sample.time_us, attitude);

	// THEN: the data should be marked as unhealthy
	EXPECT_FALSE(_range_finder.isDelayedDataHealthy());
}

TEST_F(SensorRangeFinderTest, rangeMinExceeded)
{
	Dcmf attitude{Eulerf(0.f, 0.f, 0.f)};

	// WHEN: the measured range is shorter than the minimum
	rangeSample bad_sample = _good_sample;
	bad_sample.rng = _min_range - 0.01f;
	_range_finder.setDelayedSample(bad_sample);
	_range_finder.runChecks(bad_sample.time_us, attitude);

	// THEN: the data should be marked as unhealthy
	EXPECT_FALSE(_range_finder.isDelayedDataHealthy());
}

TEST_F(SensorRangeFinderTest, outOfDate)
{
	Dcmf attitude{Eulerf(0.f, 0.f, 0.f)};

	// WHEN: the data is outdated
	rangeSample outdated_sample = _good_sample;
	outdated_sample.time_us = 0;
	uint64_t t_now = _good_sample.time_us;
	_range_finder.setDelayedSample(outdated_sample);
	_range_finder.runChecks(t_now, attitude);

	// THEN: the data should be marked as unhealthy
	EXPECT_FALSE(_range_finder.isDelayedDataHealthy());
}

TEST_F(SensorRangeFinderTest, rangeStuck)
{
	Dcmf attitude{Eulerf(0.f, 0.f, 0.f)};

	// WHEN: the data is constantly the same
	rangeSample new_sample = _good_sample;
	const float dt = 3e5;
	const float stuck_timeout = 10e6;
	for (int i = 0; i < (stuck_timeout / dt) + 2; i++) {
		_range_finder.setDelayedSample(new_sample);
		_range_finder.runChecks(new_sample.time_us, attitude);
		new_sample.time_us += float(i) * dt;
	}

	// THEN: the data should be marked as unhealthy
	EXPECT_FALSE(_range_finder.isDelayedDataHealthy());
}

TEST_F(SensorRangeFinderTest, qualityHysteresis)
{
	Dcmf attitude{Eulerf(0.f, 0.f, 0.f)};

	// WHEN: the data is first bad and then good
	rangeSample new_sample = _good_sample;

	new_sample.quality = 0;
	_range_finder.setDelayedSample(new_sample);
	_range_finder.runChecks(new_sample.time_us, attitude);
	EXPECT_FALSE(_range_finder.isDelayedDataHealthy());

	new_sample.quality = _good_sample.quality;
	_range_finder.setDelayedSample(new_sample);
	_range_finder.runChecks(new_sample.time_us, attitude);
	EXPECT_FALSE(_range_finder.isDelayedDataHealthy());

	// AND: we need to put enough good data to pass the hysteresis
	const float dt = 3e5;
	const float hyst_time = 1e6;
	for (int i = 0; i < (hyst_time / dt) + 1; i++) {
		_range_finder.setDelayedSample(new_sample);
		_range_finder.runChecks(new_sample.time_us, attitude);
		new_sample.time_us += float(i) * dt;
	}

	// THEN: the data is again declared healthy
	EXPECT_TRUE(_range_finder.isDelayedDataHealthy());
}

/* TEST_F(SensorRangeFinderTest, continuity) */
/* { */
/* 	Dcmf attitude{Eulerf(0.f, 0.f, 0.f)}; */

/* 	// WHEN: the data too slow */
/* 	rangeSample new_sample = _good_sample; */
/* 	float dt = 2e6; */
/* 	for (int i = 0; i < 10; i++) { */
/* 		_range_finder.setDelayedSample(new_sample); */
/* 		_range_finder.runChecks(new_sample.time_us, attitude); */
/* 		new_sample.time_us += float(i) * dt; */
/* 	} */

/* 	// THEN: the data should be marked as unhealthy */
/* 	EXPECT_TRUE(_range_finder.isDelayedDataHealthy()); */
/* } */
