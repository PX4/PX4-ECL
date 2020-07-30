/****************************************************************************
 *
 *   Copyright (c) 2012-2014 PX4 Development Team. All rights reserved.
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
 * @file map_projection.cpp
 *
 * Geo / math functions to perform geodesic calculations
 *
 * @author Thomas Gubler <thomasgubler@student.ethz.ch>
 * @author Julian Oes <joes@student.ethz.ch>
 * @author Lorenz Meier <lm@inf.ethz.ch>
 * @author Anton Babushkin <anton.babushkin@me.com>
 */

#include "map_projection.h"
#include <ecl.h>

#include <mathlib/mathlib.h>
#include <matrix/math.hpp>
#include <float.h>

using matrix::wrap_pi;
using matrix::wrap_2pi;

static constexpr double CONSTANTS_RADIUS_OF_EARTH = 6371000;					// meters (m)

/*
 * Azimuthal Equidistant Projection
 * formulas according to: http://mathworld.wolfram.com/AzimuthalEquidistantProjection.html
 */

static struct map_projection_reference_s mp_ref;
static struct globallocal_converter_reference_s gl_ref = {0.0f, false};

bool map_projection_global_initialized()
{
	return map_projection_initialized(&mp_ref);
}

bool map_projection_initialized(const struct map_projection_reference_s *ref)
{
	return ref->init_done;
}

uint64_t map_projection_global_timestamp()
{
	return map_projection_timestamp(&mp_ref);
}

uint64_t map_projection_timestamp(const struct map_projection_reference_s *ref)
{
	return ref->timestamp;
}

// lat_0, lon_0 are expected to be in correct format: -> 47.1234567 and not 471234567
int map_projection_global_init(double lat_0, double lon_0, uint64_t timestamp)
{
	return map_projection_init_timestamped(&mp_ref, lat_0, lon_0, timestamp);
}

// lat_0, lon_0 are expected to be in correct format: -> 47.1234567 and not 471234567
int map_projection_init_timestamped(struct map_projection_reference_s *ref, double lat_0, double lon_0, uint64_t timestamp)
{

	ref->lat_rad = math::radians(lat_0);
	ref->lon_rad = math::radians(lon_0);
	ref->sin_lat = sin(ref->lat_rad);
	ref->cos_lat = cos(ref->lat_rad);

	ref->timestamp = timestamp;
	ref->init_done = true;

	return 0;
}

//lat_0, lon_0 are expected to be in correct format: -> 47.1234567 and not 471234567
int map_projection_init(struct map_projection_reference_s *ref, double lat_0, double lon_0)
{
	return map_projection_init_timestamped(ref, lat_0, lon_0, ecl_absolute_time());
}

int map_projection_global_reference(double *ref_lat_rad, double *ref_lon_rad)
{
	return map_projection_reference(&mp_ref, ref_lat_rad, ref_lon_rad);
}

int map_projection_reference(const struct map_projection_reference_s *ref, double *ref_lat_rad, double *ref_lon_rad)
{
	if (!map_projection_initialized(ref)) {
		return -1;
	}

	*ref_lat_rad = ref->lat_rad;
	*ref_lon_rad = ref->lon_rad;

	return 0;
}

int map_projection_global_project(double lat, double lon, float *x, float *y)
{
	return map_projection_project(&mp_ref, lat, lon, x, y);
}

int map_projection_project(const struct map_projection_reference_s *ref, double lat, double lon, float *x, float *y)
{
	if (!map_projection_initialized(ref)) {
		return -1;
	}

	const double lat_rad = math::radians(lat);
	const double lon_rad = math::radians(lon);

	const double sin_lat = sin(lat_rad);
	const double cos_lat = cos(lat_rad);

	const double cos_d_lon = cos(lon_rad - ref->lon_rad);

	const double arg = math::constrain(ref->sin_lat * sin_lat + ref->cos_lat * cos_lat * cos_d_lon, -1.0,  1.0);
	const double c = acos(arg);

	double k = 1.0;

	if (fabs(c) > 0) {
		k = (c / sin(c));
	}

	*x = static_cast<float>(k * (ref->cos_lat * sin_lat - ref->sin_lat * cos_lat * cos_d_lon) * CONSTANTS_RADIUS_OF_EARTH);
	*y = static_cast<float>(k * cos_lat * sin(lon_rad - ref->lon_rad) * CONSTANTS_RADIUS_OF_EARTH);

	return 0;
}

int map_projection_global_reproject(float x, float y, double *lat, double *lon)
{
	return map_projection_reproject(&mp_ref, x, y, lat, lon);
}

int map_projection_reproject(const struct map_projection_reference_s *ref, float x, float y, double *lat, double *lon)
{
	if (!map_projection_initialized(ref)) {
		return -1;
	}

	const double x_rad = (double)x / CONSTANTS_RADIUS_OF_EARTH;
	const double y_rad = (double)y / CONSTANTS_RADIUS_OF_EARTH;
	const double c = sqrt(x_rad * x_rad + y_rad * y_rad);

	if (fabs(c) > 0) {
		const double sin_c = sin(c);
		const double cos_c = cos(c);

		const double lat_rad = asin(cos_c * ref->sin_lat + (x_rad * sin_c * ref->cos_lat) / c);
		const double lon_rad = (ref->lon_rad + atan2(y_rad * sin_c, c * ref->cos_lat * cos_c - x_rad * ref->sin_lat * sin_c));

		*lat = math::degrees(lat_rad);
		*lon = math::degrees(lon_rad);

	} else {
		*lat = math::degrees(ref->lat_rad);
		*lon = math::degrees(ref->lon_rad);
	}

	return 0;
}

int map_projection_global_getref(double *lat_0, double *lon_0)
{
	if (!map_projection_global_initialized()) {
		return -1;
	}

	if (lat_0 != nullptr) {
		*lat_0 = math::degrees(mp_ref.lat_rad);
	}

	if (lon_0 != nullptr) {
		*lon_0 = math::degrees(mp_ref.lon_rad);
	}

	return 0;

}
int globallocalconverter_init(double lat_0, double lon_0, float alt_0, uint64_t timestamp)
{
	gl_ref.alt = alt_0;

	if (!map_projection_global_init(lat_0, lon_0, timestamp)) {
		gl_ref.init_done = true;
		return 0;
	}

	gl_ref.init_done = false;
	return -1;
}

bool globallocalconverter_initialized()
{
	return gl_ref.init_done && map_projection_global_initialized();
}

int globallocalconverter_tolocal(double lat, double lon, float alt, float *x, float *y, float *z)
{
	if (!map_projection_global_initialized()) {
		return -1;
	}

	map_projection_global_project(lat, lon, x, y);
	*z = gl_ref.alt - alt;

	return 0;
}

int globallocalconverter_toglobal(float x, float y, float z,  double *lat, double *lon, float *alt)
{
	if (!map_projection_global_initialized()) {
		return -1;
	}

	map_projection_global_reproject(x, y, lat, lon);
	*alt = gl_ref.alt - z;

	return 0;
}

int globallocalconverter_getref(double *lat_0, double *lon_0, float *alt_0)
{
	if (map_projection_global_initialized() != 0) {
		return -1;
	}

	if (map_projection_global_getref(lat_0, lon_0)) {
		return -1;
	}

	if (alt_0 != nullptr) {
		*alt_0 = gl_ref.alt;
	}

	return 0;
}

