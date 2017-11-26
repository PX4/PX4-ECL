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
 * @file ecl_controller.cpp
 * Definition of base class for other controllers
 *
 * @author Lorenz Meier <lm@inf.ethz.ch>
 * @author Thomas Gubler <thomasgubler@gmail.com>
 *
 * Acknowledgements:
 *
 *   The control design is based on a design
 *   by Paul Riseborough and Andrew Tridgell, 2013,
 *   which in turn is based on initial work of
 *   Jonathan Challinger, 2012.
 */

#include "ecl_controller.h"

using math::constrain;
using math::max;
using math::min;

void ECL_Controller::integrate(const float rate_error)
{
	const float dt = ecl_elapsed_time(&_last_run) * 1e-6f;
	_last_run = ecl_absolute_time();

	// lock integral for long intervals
	if ((_k_i > 0.0f) && (dt < 0.5f)) {
		float id = rate_error * dt;

		// anti-windup: do not allow integrator to increase if actuator is at limit
		if (_last_output < -1.0f) {
			// only allow motion to center: increase value
			id = max(id, 0.0f);

		} else if (_last_output > 1.0f) {
			// only allow motion to center: decrease value
			id = min(id, 0.0f);
		}

		// add and constrain
		_integrator = constrain(_integrator + id * _k_i, -_integrator_max, _integrator_max);
	}
}
