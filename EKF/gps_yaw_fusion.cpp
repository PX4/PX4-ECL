/****************************************************************************
 *
 *   Copyright (c) 2018 Estimation and Control Library (ECL). All rights reserved.
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
 * @file gps_yaw_fusion.cpp
 * Definition of functions required to use yaw obtained from GPS dual antenna measurements.
 *
 * @author Paul Riseborough <p_riseborough@live.com.au>
 *
 */

#include "ekf.h"

#include <ecl.h>
#include <mathlib/mathlib.h>
#include <cstdlib>

void Ekf::fuseGpsYaw()
{
	// calculate the observed yaw angle of antenna array, converting a from body to antenna yaw measurement
	const float measured_hdg = wrap_pi(_gps_sample_delayed.yaw);

	// define the predicted antenna array vector and rotate into earth frame
	const Vector3f ant_vec_bf = {cosf(_gps_yaw_offset), sinf(_gps_yaw_offset), 0.0f};
	const Vector3f ant_vec_ef = _R_to_earth * ant_vec_bf;

	// check if antenna array vector is within 30 degrees of vertical and therefore unable to provide a reliable heading
	if (fabsf(ant_vec_ef(2)) > cosf(math::radians(30.0f)))  {
		return;
	}

	// using magnetic heading tuning parameter
	const float R_YAW = sq(fmaxf(_params.mag_heading_noise, 1.0e-2f));

	// determine if a 321 or 312 Euler sequence is best
	bool success;
	if (fabsf(_R_to_earth(2, 0)) < fabsf(_R_to_earth(2, 1))) {
		success = fuseYaw321(measured_hdg, R_YAW);

	} else {
		success = fuseYaw312(measured_hdg, R_YAW);
	}

	if (success) {
		_time_last_gps_yaw_fuse = _time_last_imu;
	}
}

bool Ekf::resetYawToGps()
{
	// define the predicted antenna array vector and rotate into earth frame
	const Vector3f ant_vec_bf = {cosf(_gps_yaw_offset), sinf(_gps_yaw_offset), 0.0f};
	const Vector3f ant_vec_ef = _R_to_earth * ant_vec_bf;

	// check if antenna array vector is within 30 degrees of vertical and therefore unable to provide a reliable heading
	if (fabsf(ant_vec_ef(2)) > cosf(math::radians(30.0f)))  {
		return false;
	}

	// GPS yaw measurement is alreday compensated for antenna offset in the driver
	const float measured_yaw = _gps_sample_delayed.yaw;

	const float yaw_variance = sq(fmaxf(_params.mag_heading_noise, 1.0e-2f));
	resetQuatStateYaw(measured_yaw, yaw_variance, true);

	_time_last_gps_yaw_fuse = _time_last_imu;

	return true;
}
