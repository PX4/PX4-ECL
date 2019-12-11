/****************************************************************************
 *
 *   Copyright (c) 2013 Estimation and Control Library (ECL). All rights reserved.
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
 * @file gps_checks.cpp
 * Perform pre-flight and in-flight GPS quality checks
 *
 * @author Paul Riseborough <p_riseborough@live.com.au>
 *
 */

#include "ekf.h"

#include <ecl.h>
#include <geo_lookup/geo_mag_declination.h>
#include <mathlib/mathlib.h>

// GPS pre-flight check bit locations
#define MASK_GPS_NSATS  (1<<0)
#define MASK_GPS_PDOP   (1<<1)
#define MASK_GPS_HACC   (1<<2)
#define MASK_GPS_VACC   (1<<3)
#define MASK_GPS_SACC   (1<<4)
#define MASK_GPS_HDRIFT (1<<5)
#define MASK_GPS_VDRIFT (1<<6)
#define MASK_GPS_HSPD   (1<<7)
#define MASK_GPS_VSPD   (1<<8)

bool Ekf::collect_gps(const gps_message &gps)
{
	// Run GPS checks always
	_gps_checks_passed = gps_is_good(gps);
	if (!_NED_origin_initialised && _gps_checks_passed) {
		// If we have good GPS data set the origin's WGS-84 position to the last gps fix
		double lat = gps.lat / 1.0e7;
		double lon = gps.lon / 1.0e7;
		map_projection_init_timestamped(&_pos_ref, lat, lon, _time_last_imu);

		// if we are already doing aiding, correct for the change in position since the EKF started navigationg
		if (_control_status.flags.opt_flow || _control_status.flags.gps || _control_status.flags.ev_pos || _control_status.flags.ev_vel) {
			double est_lat, est_lon;
			map_projection_reproject(&_pos_ref, -_state.pos(0), -_state.pos(1), &est_lat, &est_lon);
			map_projection_init_timestamped(&_pos_ref, est_lat, est_lon, _time_last_imu);
		}

		// Take the current GPS height and subtract the filter height above origin to estimate the GPS height of the origin
		_gps_alt_ref = 1e-3f * (float)gps.alt + _state.pos(2);
		_NED_origin_initialised = true;
		_last_gps_origin_time_us = _time_last_imu;

		// set the magnetic field data returned by the geo library using the current GPS position
		_mag_declination_gps = math::radians(get_mag_declination(lat, lon));
		_mag_inclination_gps = math::radians(get_mag_inclination(lat, lon));
		_mag_strength_gps = 0.01f * get_mag_strength(lat, lon);

		// request a reset of the yaw using the new declination
		_mag_yaw_reset_req = true;
		// save the horizontal and vertical position uncertainty of the origin
		_gps_origin_eph = gps.eph;
		_gps_origin_epv = gps.epv;

		// if the user has selected GPS as the primary height source, switch across to using it
		if (_params.vdist_sensor_type == VDIST_SENSOR_GPS) {
			ECL_INFO_TIMESTAMPED("EKF GPS checks passed (WGS-84 origin set, using GPS height)");
			_control_status.flags.baro_hgt = false;
			_control_status.flags.gps_hgt = true;
			_control_status.flags.rng_hgt = false;
			// zero the sensor offset
			_hgt_sensor_offset = 0.0f;
		} else {
			ECL_INFO_TIMESTAMPED("EKF GPS checks passed (WGS-84 origin set)");
		}
	}

	// start collecting GPS if there is a 3D fix and the NED origin has been set
	return _NED_origin_initialised && (gps.fix_type >= 3);
}

/*
 * Return true if the GPS solution quality is adequate to set an origin for the EKF
 * and start GPS aiding.
 * All activated checks must pass for 10 seconds.
 * Checks are activated using the EKF2_GPS_CHECK bitmask parameter
 * Checks are adjusted using the EKF2_REQ_* parameters
*/
bool Ekf::gps_is_good(const gps_message &gps)
{
	// Check the fix type
	_gps_check_fail_status.flags.fix = (gps.fix_type < 3);

	// Check the number of satellites
	_gps_check_fail_status.flags.nsats = (gps.nsats < _params.req_nsats);

	// Check the position dilution of precision
	_gps_check_fail_status.flags.pdop = (gps.pdop > _params.req_pdop);

	// Check the reported horizontal and vertical position accuracy
	_gps_check_fail_status.flags.hacc = (gps.eph > _params.req_hacc);
	_gps_check_fail_status.flags.vacc = (gps.epv > _params.req_vacc);

	// Check the reported speed accuracy
	_gps_check_fail_status.flags.sacc = (gps.sacc > _params.req_sacc);

	// check if GPS quality is degraded
	_gps_error_norm = fmaxf((gps.eph / _params.req_hacc) , (gps.epv / _params.req_vacc));
	_gps_error_norm = fmaxf(_gps_error_norm , (gps.sacc / _params.req_sacc));

	// Calculate time lapsed since last update, limit to prevent numerical errors and calculate a lowpass filter coefficient
	const float filt_time_const = 10.0f;
	float dt = fminf(fmaxf(float(_time_last_imu - _gps_pos_prev.timestamp) * 1e-6f, 0.001f), filt_time_const);
	float filter_coef = dt / filt_time_const;

	// The following checks are only valid when the vehicle is at rest
	double lat = gps.lat * 1.0e-7;
	double lon = gps.lon * 1.0e-7;
	if (!_control_status.flags.in_air && _vehicle_at_rest) {
		// Calculate position movement since last measurement
		float delta_posN = 0.0f;
		float delta_PosE = 0.0f;

		// calculate position movement since last GPS fix
		if (_gps_pos_prev.timestamp > 0) {
			map_projection_project(&_gps_pos_prev, lat, lon, &delta_posN, &delta_PosE);

		} else {
			// no previous position has been set
			map_projection_init_timestamped(&_gps_pos_prev, lat, lon, _time_last_imu);
			_gps_alt_prev = 1e-3f * (float)gps.alt;

		}

		// Calculate the horizontal drift velocity components and limit to 10x the threshold
		float vel_limit = 10.0f * _params.req_hdrift;
		float velN = fminf(fmaxf(delta_posN / dt, -vel_limit), vel_limit);
		float velE = fminf(fmaxf(delta_PosE / dt, -vel_limit), vel_limit);

		// Apply a low pass filter
		_gpsDriftVelN = velN * filter_coef + _gpsDriftVelN * (1.0f - filter_coef);
		_gpsDriftVelE = velE * filter_coef + _gpsDriftVelE * (1.0f - filter_coef);

		// Calculate the horizontal drift speed and fail if too high
		_gps_drift_metrics[0] = sqrtf(_gpsDriftVelN * _gpsDriftVelN + _gpsDriftVelE * _gpsDriftVelE);
		_gps_check_fail_status.flags.hdrift = (_gps_drift_metrics[0] > _params.req_hdrift);

		// Calculate the vertical drift velocity and limit to 10x the threshold
		float vz_drift_limit = 10.0f * _params.req_vdrift;
		float gps_alt_m = 1e-3f * (float)gps.alt;
		float velD = math::constrain(((_gps_alt_prev - gps_alt_m) / dt), -vz_drift_limit, vz_drift_limit);

		// Apply a low pass filter to the vertical velocity
		_gps_drift_velD = velD * filter_coef + _gps_drift_velD * (1.0f - filter_coef);

		// Fail if the vertical drift speed is too high
		_gps_drift_metrics[1] = fabsf(_gps_drift_velD);
		_gps_check_fail_status.flags.vdrift = (_gps_drift_metrics[1] > _params.req_vdrift);

		// Check the magnitude of the filtered horizontal GPS velocity
		float vxy_drift_limit = 10.0f * _params.req_hdrift;
		float gps_velN = fminf(fmaxf(gps.vel_ned[0], -vxy_drift_limit), vxy_drift_limit);
		float gps_velE = fminf(fmaxf(gps.vel_ned[1], -vxy_drift_limit), vxy_drift_limit);
		_gps_velN_filt = gps_velN * filter_coef + _gps_velN_filt * (1.0f - filter_coef);
		_gps_velE_filt  = gps_velE * filter_coef + _gps_velE_filt  * (1.0f - filter_coef);
		_gps_drift_metrics[2] = sqrtf(_gps_velN_filt * _gps_velN_filt + _gps_velE_filt * _gps_velE_filt);
		_gps_check_fail_status.flags.hspeed = (_gps_drift_metrics[2] > _params.req_hdrift);

		_gps_drift_updated = true;

	} else if (_control_status.flags.in_air) {
		// These checks are always declared as passed when flying
		// If on ground and moving, the last result before movement commenced is kept
		_gps_check_fail_status.flags.hdrift = false;
		_gps_check_fail_status.flags.vdrift = false;
		_gps_check_fail_status.flags.hspeed = false;
		_gps_drift_updated = false;

	} else {
		// This is the case where the vehicle is on ground and IMU movement is blocking the drift calculation
		_gps_drift_updated = true;

	}

	// save GPS fix for next time
	map_projection_init_timestamped(&_gps_pos_prev, lat, lon, _time_last_imu);
	_gps_alt_prev = 1e-3f * (float)gps.alt;

	// Check  the filtered difference between GPS and EKF vertical velocity
	float vz_diff_limit = 10.0f * _params.req_vdrift;
	float vertVel = fminf(fmaxf((gps.vel_ned[2] - _state.vel(2)), -vz_diff_limit), vz_diff_limit);
	_gps_velD_diff_filt = vertVel * filter_coef + _gps_velD_diff_filt * (1.0f - filter_coef);
	_gps_check_fail_status.flags.vspeed = (fabsf(_gps_velD_diff_filt) > _params.req_vdrift);

	// assume failed first time through
	if (_last_gps_fail_us == 0) {
		_last_gps_fail_us = _time_last_imu;
	}

	// if any user selected checks have failed, record the fail time
	if (
		_gps_check_fail_status.flags.fix ||
		(_gps_check_fail_status.flags.nsats   && (_params.gps_check_mask & MASK_GPS_NSATS)) ||
		(_gps_check_fail_status.flags.pdop    && (_params.gps_check_mask & MASK_GPS_PDOP)) ||
		(_gps_check_fail_status.flags.hacc    && (_params.gps_check_mask & MASK_GPS_HACC)) ||
		(_gps_check_fail_status.flags.vacc    && (_params.gps_check_mask & MASK_GPS_VACC)) ||
		(_gps_check_fail_status.flags.sacc    && (_params.gps_check_mask & MASK_GPS_SACC)) ||
		(_gps_check_fail_status.flags.hdrift  && (_params.gps_check_mask & MASK_GPS_HDRIFT)) ||
		(_gps_check_fail_status.flags.vdrift  && (_params.gps_check_mask & MASK_GPS_VDRIFT)) ||
		(_gps_check_fail_status.flags.hspeed  && (_params.gps_check_mask & MASK_GPS_HSPD)) ||
		(_gps_check_fail_status.flags.vspeed  && (_params.gps_check_mask & MASK_GPS_VSPD))
	) {
		_last_gps_fail_us = _time_last_imu;
	} else {
		_last_gps_pass_us = _time_last_imu;
	}

	// continuous period without fail of x seconds required to return a healthy status
	return _time_last_imu - _last_gps_fail_us > (uint64_t)_min_gps_health_time_us;
}
