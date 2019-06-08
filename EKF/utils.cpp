

#include <ecl.h>

#include "utils.hpp"

matrix::Dcm<ecl_float_t> taitBryan312ToRotMat(const matrix::Vector3<ecl_float_t> &rot312)
{
	// Calculate the frame2 to frame 1 rotation matrix from a 312 Tait-Bryan rotation sequence
	const ecl_float_t c2 = cosf(rot312(2)); // third rotation is pitch
	const ecl_float_t s2 = sinf(rot312(2));
	const ecl_float_t s1 = sinf(rot312(1)); // second rotation is roll
	const ecl_float_t c1 = cosf(rot312(1));
	const ecl_float_t s0 = sinf(rot312(0)); // first rotation is yaw
	const ecl_float_t c0 = cosf(rot312(0));

	matrix::Dcm<ecl_float_t> R;
	R(0, 0) = c0 * c2 - s0 * s1 * s2;
	R(1, 1) = c0 * c1;
	R(2, 2) = c2 * c1;
	R(0, 1) = -c1 * s0;
	R(0, 2) = s2 * c0 + c2 * s1 * s0;
	R(1, 0) = c2 * s0 + s2 * s1 * c0;
	R(1, 2) = s0 * s2 - s1 * c0 * c2;
	R(2, 0) = -s2 * c1;
	R(2, 1) = s1;

	return R;
}

ecl_float_t kahanSummation(ecl_float_t sum_previous, ecl_float_t input, ecl_float_t &accumulator)
{
	const ecl_float_t y = input - accumulator;
	const ecl_float_t t = sum_previous + y;
	accumulator = (t - sum_previous) - y;
	return t;
}

matrix::Dcm<ecl_float_t> quatToInverseRotMat(const matrix::Quaternion<ecl_float_t> &quat)
{
	const auto q00 = quat(0) * quat(0);
	const auto q11 = quat(1) * quat(1);
	const auto q22 = quat(2) * quat(2);
	const auto q33 = quat(3) * quat(3);
	const auto q01 = quat(0) * quat(1);
	const auto q02 = quat(0) * quat(2);
	const auto q03 = quat(0) * quat(3);
	const auto q12 = quat(1) * quat(2);
	const auto q13 = quat(1) * quat(3);
	const auto q23 = quat(2) * quat(3);

	matrix::Dcm<ecl_float_t> dcm;
	dcm(0, 0) = q00 + q11 - q22 - q33;
	dcm(1, 1) = q00 - q11 + q22 - q33;
	dcm(2, 2) = q00 - q11 - q22 + q33;
	dcm(1, 0) = 2.0 * (q12 - q03);
	dcm(2, 0) = 2.0 * (q13 + q02);
	dcm(0, 1) = 2.0 * (q12 + q03);
	dcm(2, 1) = 2.0 * (q23 - q01);
	dcm(0, 2) = 2.0 * (q13 - q02);
	dcm(1, 2) = 2.0 * (q23 + q01);

	return dcm;
}

bool shouldUse321RotationSequence(const matrix::Dcm<ecl_float_t> &R)
{
	return fabsf(R(2, 0)) < fabsf(R(2, 1));
}

ecl_float_t getEuler321Yaw(const matrix::Quaternion<ecl_float_t> &q)
{
	// Values from yaw_input_321.c file produced by
	// https://github.com/PX4/ecl/blob/master/matlab/scripts/Inertial%20Nav%20EKF/quat2yaw321.m
	const float a = 2.f * (q(0) * q(3) + q(1) * q(2));
	const float b = sq(q(0)) + sq(q(1)) - sq(q(2)) - sq(q(3));
	return atan2f(a, b);
}

ecl_float_t getEuler321Yaw(const matrix::Dcm<ecl_float_t> &R)
{
	return atan2f(R(1, 0), R(0, 0));
}

ecl_float_t getEuler312Yaw(const matrix::Quaternion<ecl_float_t> &q)
{
	// Values from yaw_input_312.c file produced by
	// https://github.com/PX4/ecl/blob/master/matlab/scripts/Inertial%20Nav%20EKF/quat2yaw312.m
	const auto a = 2 * (q(0) * q(3) - q(1) * q(2));
	const auto b = sq(q(0)) - sq(q(1)) + sq(q(2)) - sq(q(3));
	return atan2f(a, b);
}

ecl_float_t getEuler312Yaw(const matrix::Dcm<ecl_float_t> &R)
{
	return atan2f(-R(0, 1), R(1, 1));
}

matrix::Dcm<ecl_float_t> updateEuler321YawInRotMat(ecl_float_t yaw, const matrix::Dcm<ecl_float_t> &rot_in)
{
	matrix::Euler<ecl_float_t> euler321(rot_in);
	euler321(2) = yaw;
	return matrix::Dcm<ecl_float_t>(euler321);
}

matrix::Dcm<ecl_float_t> updateEuler312YawInRotMat(ecl_float_t yaw, const matrix::Dcm<ecl_float_t> &rot_in)
{
	const matrix::Vector3<ecl_float_t> rotVec312(yaw,  // yaw
			asinf(rot_in(2, 1)),  // roll
			atan2f(-rot_in(2, 0), rot_in(2, 2)));  // pitch
	return taitBryan312ToRotMat(rotVec312);
}

matrix::Dcm<ecl_float_t> updateYawInRotMat(ecl_float_t yaw, const matrix::Dcm<ecl_float_t> &rot_in)
{
	return shouldUse321RotationSequence(rot_in) ?
	       updateEuler321YawInRotMat(yaw, rot_in) :
	       updateEuler312YawInRotMat(yaw, rot_in);
}
