/****************************************************************************
 *
 *   Copyright (c) 2017 PX4 Development Team. All rights reserved.
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

/**
 * @file mag_validator.cpp
 *
 * A data validation class to identify anomalies in magnetometer data streams
 *
 * @author Stephan Brown <Stephan.Brown.07@gmail.com>
 */

#include "mag_validator.h"

#include <uORB/topics/vehicle_attitude.h>
#include <uORB/topics/vehicle_global_position.h>

#include <geo_lookup/geo_mag_declination.h>

#include <matrix/matrix/math.hpp>
#include <ecl/ecl.h>
#include <math.h>

MagValidator::MagValidator(DataValidator *prev_sibling) :
		DataValidator(prev_sibling),
		_att_sub(-1),
		_gpos_sub(-1),
		_error_filter(150, 1),
		_local_inc(0.0f),
		_mag_inc(0.0f),
		_error_filtered(0.0f)
	{ }

void MagValidator::complex_error_update()
{
	if (_att_sub < 0) {
		ECL_INFO("Subscribing to vehicle attitude");
		_att_sub = orb_subscribe(ORB_ID(vehicle_attitude));
	}

	if (_gpos_sub < 0) {
		ECL_INFO("Subscribing to global position");
		_gpos_sub = orb_subscribe(ORB_ID(vehicle_global_position));
	}

	bool updated = false;

	static struct vehicle_global_position_s gpos = {};
	int ret = orb_check(_gpos_sub, &updated);
	if (ret == OK && updated) {
		orb_copy(ORB_ID(vehicle_global_position), _gpos_sub, &gpos);
		if (gpos.timestamp > 0) {
			_local_inc = get_mag_inclination(gpos.lat, gpos.lon);
		}
	}

	struct vehicle_attitude_s att = {};
	ret = orb_check(_att_sub, &updated);
	if (ret == OK && updated) {
		orb_copy(ORB_ID(vehicle_attitude), _att_sub, &att);

		matrix::Vector<float, 3> mag_vect_bf(_value);
		matrix::Quaternionf att_q(att.q);
		matrix::Vector<float, 3> mag_vect_wf = att_q.conjugate(mag_vect_bf);

		_mag_inc = 90.0f - acosf((mag_vect_wf(2)/mag_vect_wf.length())) * (180.0f / M_PI_F);
		if (gpos.timestamp > 0) {
			_error_filtered = _error_filter.apply(fabsf(_local_inc - _mag_inc));
		}
	}
}

bool MagValidator::complex_error_state()
{
	return _error_filtered > 20.0f;
}

void
MagValidator::print()
{
	DataValidator::print();
	ECL_INFO("\tLcl: %.2f, Act: %.2f Err: %.2f",
		(double)_local_inc, (double)_mag_inc, (double)(_error_filtered));
}
