
#include <ecl.h>

#include <matrix/math.hpp>

#pragma once

// return the square of two floating point numbers - used in auto coded sections
template<typename T> static constexpr T sq(T var) { return var * var; }

// converts Tait-Bryan 312 sequence of rotations from frame 1 to frame 2
// to the corresponding rotation matrix that rotates from frame 2 to frame 1
// rot312(0) - First rotation is a RH rotation about the Z axis (rad)
// rot312(1) - Second rotation is a RH rotation about the X axis (rad)
// rot312(2) - Third rotation is a RH rotation about the Y axis (rad)
// See http://www.atacolorado.com/eulersequences.doc
matrix::Dcm<ecl_float_t> taitBryan312ToRotMat(const matrix::Vector3<ecl_float_t> &rot312);

// Use Kahan summation algorithm to get the sum of "sum_previous" and "input".
// This function relies on the caller to be responsible for keeping a copy of
// "accumulator" and passing this value at the next iteration.
// Ref: https://en.wikipedia.org/wiki/Kahan_summation_algorithm
ecl_float_t kahanSummation(ecl_float_t sum_previous, ecl_float_t input, ecl_float_t &accumulator);

// calculate the inverse rotation matrix from a quaternion rotation
// this produces the inverse rotation to that produced by the math library quaternion to matrix::Dcm<ecl_float_t> operator
matrix::Dcm<ecl_float_t> quatToInverseRotMat(const matrix::Quaternion<ecl_float_t> &quat);

// We should use a 3-2-1 Tait-Bryan (yaw-pitch-roll) rotation sequence
// when there is more roll than pitch tilt and a 3-1-2 rotation sequence
// when there is more pitch than roll tilt to avoid gimbal lock.
bool shouldUse321RotationSequence(const matrix::Dcm<ecl_float_t> &R);

ecl_float_t getEuler321Yaw(const matrix::Quaternion<ecl_float_t> &q);
ecl_float_t getEuler321Yaw(const matrix::Dcm<ecl_float_t> &R);

ecl_float_t getEuler312Yaw(const matrix::Quaternion<ecl_float_t> &q);
ecl_float_t getEuler312Yaw(const matrix::Dcm<ecl_float_t> &R);

matrix::Dcm<ecl_float_t> updateEuler321YawInRotMat(ecl_float_t yaw, const matrix::Dcm<ecl_float_t> &rot_in);
matrix::Dcm<ecl_float_t> updateEuler312YawInRotMat(ecl_float_t yaw, const matrix::Dcm<ecl_float_t> &rot_in);

// Checks which euler rotation sequence to use and update yaw in rotation matrix
matrix::Dcm<ecl_float_t> updateYawInRotMat(ecl_float_t yaw, const matrix::Dcm<ecl_float_t> &rot_in);

namespace ecl
{
inline float powf(float x, int exp)
{
	float ret;

	if (exp > 0) {
		ret = x;

		for (int count = 1; count < exp; count++) {
			ret *= x;
		}

		return ret;

	} else if (exp < 0) {
		return 1.0f / ecl::powf(x, -exp);
	}

	return 1.0f;
}
}
