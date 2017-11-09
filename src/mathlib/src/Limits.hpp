/****************************************************************************
 *
 *   Copyright (c) 2013 PX4 Development Team. All rights reserved.
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
 * @file Limits.hpp
 *
 * Limiting / constrain helper functions
 */

#pragma once

#include <stdio.h>
#include <cmath>

namespace math
{

template<typename _Tp>
bool is_finite(_Tp x) {
#if defined (__PX4_NUTTX)
    return PX4_ISFINITE(x);
#elif defined (__PX4_QURT)
    return __builtin_isfinite(x);
#else
    return std::isfinite(x);
#endif
}


template<typename _Tp>
inline _Tp min(const _Tp a, const _Tp b)
{
	return (a < b) ? a : b;
}

template<typename _Tp>
inline _Tp max(const _Tp a, const _Tp b)
{
	return (a > b) ? a : b;
}

template<typename _Tp>
inline _Tp constrain(const _Tp val, const _Tp min_val, const _Tp max_val)
{
	return (val < min_val) ? min_val : ((val > max_val) ? max_val : val);
}

template<typename _Tp>
inline _Tp radians(const _Tp degrees)
{
	return (degrees / ((_Tp)180) * (_Tp)3.14159265358979323);
}

template<typename _Tp>
inline _Tp degrees(const _Tp radians)
{
	return (radians / ((_Tp)3.14159265358979323) * (_Tp)180);
}

template<typename _Tp>
_Tp wrap_pi(_Tp radians)
{
	/* value is inf or NaN */
	if (!is_finite(radians)) {
		return radians;
	}

	/* define pi */
	_Tp Pi = (_Tp)3.14159265358979323;

	int c = 0;

	while (radians >= Pi) {

		radians -= (_Tp)2 * Pi;

		if (c++ > 3) {
			return NAN;
		}
	}

	c = 0;

	while (radians < Pi) {
		radians += (_Tp)2 * Pi;

		if (c++ > 3) {
			return NAN;
		}
	}

	return radians;
}

template<typename _Tp>
_Tp wrap_2pi(_Tp radians)
{
	/* value is inf or NaN */
	if (!is_finite(radians)) {
		return radians;
	}

	/* define pi */
	_Tp Pi = (_Tp)3.14159265358979323;

	int c = 0;

	while (radians >= (_Tp)2 * Pi) {
		radians -= (_Tp)2 * Pi;

		if (c++ > 3) {
			return NAN;
		}
	}

	c = 0;

	while (radians < (_Tp)0) {
		radians += (_Tp)2 * Pi;

		if (c++ > 3) {
			return NAN;
		}
	}

	return radians;
}

template<typename _Tp>
_Tp wrap_180(_Tp degrees)
{
	/* value is inf or NaN */
	if (!is_finite(degrees)) {
		return degrees;
	}

	int c = 0;

	while (degrees >= (_Tp)180) {
		degrees -= (_Tp)360;

		if (c++ > 3) {
			return NAN;
		}
	}

	c = 0;

	while (degrees < (_Tp)-180) {
		degrees += (_Tp)360;

		if (c++ > 3) {
			return NAN;
		}
	}

	return degrees;
}

template<typename _Tp>
_Tp wrap_360(_Tp degrees)
{
	/* value is inf or NaN */
	if (!is_finite(degrees)) {
		return degrees;
	}

	int c = 0;

	while (degrees >= (_Tp)360) {
		degrees -= (_Tp)360;

		if (c++ > 3) {
			return NAN;
		}
	}

	c = 0;

	while (degrees < (_Tp)0) {
		degrees += (_Tp)360;

		if (c++ > 3) {
			return NAN;
		}
	}

	return degrees;
}

}
