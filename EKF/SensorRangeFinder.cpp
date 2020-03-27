/****************************************************************************
 *
 *   Copyright (c) 2020 Estimation and Control Library (ECL). All rights reserved.
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
 * 3. Neither the name ECL nor the names of its contributors may be
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

/**
 * @file SensorRangeFinder.cpp
 *
 * @author Mathieu Bresciani <brescianimathieu@gmail.com>
 *
 */

#include "SensorRangeFinder.hpp"

namespace estimator
{
namespace sensor
{

void SensorRangeFinder::runChecks(const uint64_t time_delayed_us, const Dcmf &R_to_earth)
{
	// calculate 2,2 element of rotation matrix from sensor frame to earth frame
	// this is required for use of range finder and flow data
	_R_rng_to_earth_2_2 = R_to_earth(2, 0) * _sin_tilt_rng + R_to_earth(2, 2) * _cos_tilt_rng;

	updateRangeDataValidity(time_delayed_us);
}

void SensorRangeFinder::updateRangeDataValidity(uint64_t time_delayed_us)
{
	updateRangeDataContinuity(time_delayed_us);

	// check if out of date
	if (((time_delayed_us - _range_sample_delayed.time_us) > 2 * RNG_MAX_INTERVAL)
	    || !isRangeDataContinuous()) {
		_rng_hgt_valid = false;
		return;
	}

	// Don't run the checks unless we have retrieved new data from the buffer
	if (_range_data_ready) {
		_rng_hgt_valid = false;

		if (_range_sample_delayed.quality == 0) {
			_time_bad_rng_signal_quality = time_delayed_us;

		} else if (time_delayed_us - _time_bad_rng_signal_quality > _range_signal_hysteresis_ms) {
			const bool is_in_range = ((_range_sample_delayed.rng >= _rng_valid_min_val)
						  && (_range_sample_delayed.rng <= _rng_valid_max_val));

			if (isTiltOk() || is_in_range) {
				updateRangeDataStuck();

				if (!_is_stuck) {
					_rng_hgt_valid = true;
				}
			}
		}
	}
}

void SensorRangeFinder::updateRangeDataContinuity(uint64_t time_delayed_us)
{
	// Calculate a first order IIR low-pass filtered time of arrival between samples using a 2 second time constant.
	float alpha = 0.5f * _dt_update;
	_dt_last_range_update_filt_us = _dt_last_range_update_filt_us * (1.0f - alpha) + alpha *
					(time_delayed_us - _range_sample_delayed.time_us);

	// Apply spike protection to the filter state.
	_dt_last_range_update_filt_us = fminf(_dt_last_range_update_filt_us, 4e6f);
}

void SensorRangeFinder::updateRangeDataStuck()
{
	// Check for "stuck" range finder measurements when range was not valid for certain period
	// This handles a failure mode observed with some lidar sensors
	if (((_range_sample_delayed.time_us - _time_last_rng_ready) > (uint64_t)10e6)) {

		// require a variance of rangefinder values to check for "stuck" measurements
		if (_rng_stuck_max_val - _rng_stuck_min_val > _range_stuck_threshold) {
			_time_last_rng_ready = _range_sample_delayed.time_us;
			_rng_stuck_min_val = 0.0f;
			_rng_stuck_max_val = 0.0f;
			_is_stuck = false;

		} else {
			if (_range_sample_delayed.rng > _rng_stuck_max_val) {
				_rng_stuck_max_val = _range_sample_delayed.rng;
			}

			if (_rng_stuck_min_val < 0.1f || _range_sample_delayed.rng < _rng_stuck_min_val) {
				_rng_stuck_min_val = _range_sample_delayed.rng;
			}

			_is_stuck = true;
		}

	} else {
		_time_last_rng_ready = _range_sample_delayed.time_us;
	}
}

bool SensorRangeFinder::canBeusedAsFailover() const
{
	return false;
}

bool SensorRangeFinder::canResetOnSensor() const
{
	return false;
}

} // namespace sensor
} // namespace estimator
