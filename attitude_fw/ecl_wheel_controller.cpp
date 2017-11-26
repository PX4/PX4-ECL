/****************************************************************************
 *
 *   Copyright (c) 2013-2018 Estimation and Control Library (ECL). All rights reserved.
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
 * @file ecl_wheel_controller.cpp
 * Implementation of a simple PID wheel controller for heading tracking.
 *
 * Authors and acknowledgments in header.
 */

#include "ecl_wheel_controller.h"

#include <geo/geo.h>

using math::constrain;
using matrix::wrap_pi;

float ECL_WheelController::control_attitude(const ECL_ControlData &ctl_data)
{
	// Do not calculate control signal with bad inputs
	if (!(ISFINITE(ctl_data.yaw_setpoint) &&
	      ISFINITE(ctl_data.yaw))) {

		ECL_WARN("not controlling wheel");
		return _rate_setpoint;
	}

	// Calculate error
	const float error = wrap_pi(ctl_data.yaw_setpoint - ctl_data.yaw);

	// Apply P controller
	_rate_setpoint = constrain(error / _tc, -_max_rate, _max_rate);

	return _rate_setpoint;
}

float ECL_WheelController::control_bodyrate(const ECL_ControlData &ctl_data)
{
	// Do not calculate control signal with bad inputs
	if (!(ISFINITE(_rate_setpoint) &&
	      ISFINITE(ctl_data.body_z_rate) &&
	      ISFINITE(ctl_data.groundspeed) &&
	      ISFINITE(ctl_data.groundspeed_scaler))) {

		ECL_WARN("not controlling wheel body");
		return constrain(_last_output, -1.0f, 1.0f);
	}

	// Calculate body angular rate error
	const float rate_error = _rate_setpoint - ctl_data.body_z_rate;

	const float min_speed = 1.0f;

	if (!ctl_data.lock_integrator && (ctl_data.groundspeed > min_speed)) {
		integrate(rate_error);
	}

	// Apply PI rate controller and store non-limited output
	_last_output = (_rate_setpoint * _k_ff + rate_error * _k_p + _integrator) * ctl_data.groundspeed_scaler *
		       ctl_data.groundspeed_scaler;

	return constrain(_last_output, -1.0f, 1.0f);
}
