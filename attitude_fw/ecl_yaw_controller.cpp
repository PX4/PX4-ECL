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
 * @file ecl_yaw_controller.cpp
 * Implementation of a simple orthogonal coordinated turn yaw PID controller.
 *
 * Authors and acknowledgments in header.
 */

#include "ecl_yaw_controller.h"

#include <geo/geo.h>

using math::constrain;
using math::radians;

float ECL_YawController::control_attitude(const ECL_ControlData &ctl_data)
{
	// Do not calculate control signal with bad inputs
	if (!(ISFINITE(ctl_data.roll) &&
	      ISFINITE(ctl_data.pitch) &&
	      ISFINITE(ctl_data.roll_rate_setpoint) &&
	      ISFINITE(ctl_data.pitch_rate_setpoint))) {

		ECL_WARN("not controlling yaw");
		return _rate_setpoint;
	}

	float constrained_roll = 0.0f;
	bool inverted = false;

	// roll is used as feedforward term and inverted flight needs to be considered
	if (fabsf(ctl_data.roll) < math::radians(90.0f)) {
		// not inverted, but numerically still potentially close to infinity
		constrained_roll = constrain(ctl_data.roll, radians(-80.0f), radians(80.0f));

	} else {
		inverted = true;

		// inverted flight, constrain on the two extremes of -pi..+pi to avoid infinity
		// note: the ranges are extended by 10 deg here to avoid numeric resolution effects
		if (ctl_data.roll > 0.0f) {
			// right hemisphere
			constrained_roll = constrain(ctl_data.roll, radians(100.0f), radians(180.0f));

		} else {
			// left hemisphere
			constrained_roll = constrain(ctl_data.roll, radians(-180.0f), radians(-100.0f));
		}
	}

	constrained_roll = constrain(constrained_roll, -fabsf(ctl_data.roll_setpoint), fabsf(ctl_data.roll_setpoint));

	if (!inverted) {
		// Calculate desired yaw rate from coordinated turn constraint / (no side forces)
		const float airspeed = constrain(ctl_data.airspeed, ctl_data.airspeed_min, ctl_data.airspeed_max);
		_rate_setpoint = tanf(constrained_roll) * cosf(ctl_data.pitch) * CONSTANTS_ONE_G / airspeed;
	}

	if (!ISFINITE(_rate_setpoint)) {
		ECL_WARN("yaw rate sepoint not finite");
		_rate_setpoint = 0.0f;
	}

	return _rate_setpoint;
}

float ECL_YawController::control_bodyrate(const ECL_ControlData &ctl_data)
{
	// Do not calculate control signal with bad inputs
	if (!(ISFINITE(_bodyrate_setpoint) &&
	      ISFINITE(ctl_data.body_z_rate) &&
	      ISFINITE(ctl_data.airspeed) &&
	      ISFINITE(ctl_data.airspeed_min) &&
	      ISFINITE(ctl_data.scaler))) {

		ECL_WARN("not controlling roll body");
		return constrain(_last_output, -1.0f, 1.0f);
	}

	// Calculate body angular rate error
	const float rate_error = _bodyrate_setpoint - ctl_data.body_z_rate;

	if (!ctl_data.lock_integrator) {
		integrate(rate_error);
	}

	// Apply PI rate controller and store non-limited output
	_last_output = (_bodyrate_setpoint * _k_ff + rate_error * _k_p + _integrator) * ctl_data.scaler * ctl_data.scaler;

	return constrain(_last_output, -1.0f, 1.0f);
}

float ECL_YawController::control_euler_rate(const ECL_ControlData &ctl_data)
{
	// Do not calculate control signal with bad inputs
	if (!(ISFINITE(ctl_data.roll) &&
	      ISFINITE(ctl_data.pitch_rate_setpoint) &&
	      ISFINITE(ctl_data.pitch) &&
	      ISFINITE(ctl_data.yaw_rate_setpoint))) {

		ECL_WARN("not controlling yaw euler rate");
		return constrain(_last_output, -1.0f, 1.0f);
	}

	// Transform setpoint to body angular rates (jacobian)
	set_bodyrate_setpoint(-sinf(ctl_data.roll) * ctl_data.pitch_rate_setpoint +
			      cosf(ctl_data.roll) * cosf(ctl_data.pitch) * ctl_data.yaw_rate_setpoint);

	return control_bodyrate(ctl_data);
}
