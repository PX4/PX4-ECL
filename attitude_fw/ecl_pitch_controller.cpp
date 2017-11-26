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
 * @file ecl_pitch_controller.cpp
 * Implementation of a simple orthogonal pitch PID controller.
 *
 * Authors and acknowledgments in header.
 */

#include "ecl_pitch_controller.h"

using math::constrain;

float ECL_PitchController::control_attitude(const ECL_ControlData &ctl_data)
{
	// Do not calculate control signal with bad inputs
	if (!(ISFINITE(ctl_data.pitch_setpoint) &&
	      ISFINITE(ctl_data.pitch))) {

		ECL_WARN("not controlling pitch");
		return _rate_setpoint;
	}

	// Calculate error
	const float error = ctl_data.pitch_setpoint - ctl_data.pitch;

	// Apply P controller
	_rate_setpoint = error / _tc;

	return _rate_setpoint;
}

float ECL_PitchController::control_bodyrate(const ECL_ControlData &ctl_data)
{
	// Do not calculate control signal with bad inputs
	if (!(ISFINITE(_bodyrate_setpoint) &&
	      ISFINITE(ctl_data.body_y_rate) &&
	      ISFINITE(ctl_data.scaler))) {

		ECL_WARN("not controlling pitch body");
		return constrain(_last_output, -1.0f, 1.0f);
	}

	// Calculate body angular rate error
	const float rate_error = _bodyrate_setpoint - ctl_data.body_y_rate;

	if (!ctl_data.lock_integrator) {
		integrate(rate_error);
	}

	// Apply PI rate controller and store non-limited output
	_last_output = (_bodyrate_setpoint * _k_ff + rate_error * _k_p + _integrator) * ctl_data.scaler * ctl_data.scaler;

	return constrain(_last_output, -1.0f, 1.0f);
}

float ECL_PitchController::control_euler_rate(const ECL_ControlData &ctl_data)
{
	// Do not calculate control signal with bad inputs
	if (!(ISFINITE(ctl_data.roll) &&
	      ISFINITE(ctl_data.pitch) &&
	      ISFINITE(ctl_data.pitch_rate_setpoint) &&
	      ISFINITE(ctl_data.yaw_rate_setpoint))) {

		ECL_WARN("not controlling pitch euler rate");
		return constrain(_last_output, -1.0f, 1.0f);
	}

	// Transform setpoint to body angular rates (jacobian)
	set_bodyrate_setpoint(cosf(ctl_data.roll) * ctl_data.pitch_rate_setpoint +
			      cosf(ctl_data.pitch) * sinf(ctl_data.roll) * ctl_data.yaw_rate_setpoint);

	return control_bodyrate(ctl_data);
}
