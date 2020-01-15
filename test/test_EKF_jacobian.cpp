/****************************************************************************
 *
 *   Copyright (c) 2019 ECL Development Team. All rights reserved.
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

#include <gtest/gtest.h>
#include <cmath>
#include <matrix/math.hpp>
#include "EKF/ekf.h"

using matrix::Vector3f;
using matrix::Quaternion;
using matrix::Dual;
using matrix::Eulerf;

class EkfJacobianTest : public ::testing::Test {
 public:
	Dual<float, 4> yaw321FromQuaterion(const Quaternion<Dual<float,4 >>& q)
	{
		return atan2(Dual<float, 4>(2) * (q(0) * q(3) + q(1) * q(2)),
			      q(0) * q(0) + q(1) * q(1) - q(2) * q(2) - q(3) * q(3));
		// matrix::Dcm<Dual<float, 4>> R_to_earth(q); // this is slower
		// return atan2(R_to_earth(1, 0), R_to_earth(0, 0));
	}

	Dual<float, 4> yaw312FromQuaterion(const Quaternion<Dual<float,4 >>& q)
	{
		return atan2(Dual<float, 4>(-2) * (q(1) * q(2) - q(0) * q(3)),
			      q(0) * q(0) - q(1) * q(1) + q(2) * q(2) - q(3) * q(3));
		// matrix::Dcm<Dual<float, 4>> R_to_earth(q); // this is slower
		// return atan2(-R_to_earth(0, 1), R_to_earth(1, 1));
	}
};

TEST_F(EkfJacobianTest, jacobianWith321EulerAngles)
{

	for(float roll = -89.0f; roll <= 90.0f; roll+=45.0f)
	{
		for(float pitch = -89.0f; pitch <= 90.0f; pitch+=45.0f)
		{
			for(float heading = -180.0f; heading <= 180.0f; heading+=45.0f)
			{
				Quaternion<float> quat(Eulerf(math::radians(roll), math::radians(pitch), math::radians(heading)));

				const float q0 = quat(0);
				const float q1 = quat(1);
				const float q2 = quat(2);
				const float q3 = quat(3);

				float t9 = q0*q3;
				float t10 = q1*q2;
				float t2 = t9+t10;
				float t3 = q0*q0;
				float t4 = q1*q1;
				float t5 = q2*q2;
				float t6 = q3*q3;
				float t7 = t3+t4-t5-t6;
				float t8 = t7*t7;
				if (t8 > 1e-6f) {
					t8 = 1.0f/t8;
				} else {
					// return;
					t8 = 1.0f/t8;
				}
				float t11 = t2*t2;
				float t12 = t8*t11*4.0f;
				float t13 = t12+1.0f;
				float t14;
				if (fabsf(t13) > 1e-6f) {
					t14 = 1.0f/t13;
				} else {
					// return;
					t14 = 1.0f/t13;
				}

				matrix::Vector<float, 4> H_YAW;
				H_YAW(0) = t8*t14*(q3*t3-q3*t4+q3*t5+q3*t6+q0*q1*q2*2.0f)*-2.0f;
				H_YAW(1) = t8*t14*(-q2*t3+q2*t4+q2*t5+q2*t6+q0*q1*q3*2.0f)*-2.0f;
				H_YAW(2) = t8*t14*(q1*t3+q1*t4+q1*t5-q1*t6+q0*q2*q3*2.0f)*2.0f;
				H_YAW(3) = t8*t14*(q0*t3+q0*t4-q0*t5+q0*t6+q1*q2*q3*2.0f)*2.0f;


				// alternative way
				using D4 = Dual<float, 4>;
				Quaternion<Dual<float, 4>> q(D4(q0, 0),D4(q1, 1),D4(q2, 2),D4(q3, 3));
				Dual<float, 4> yaw = yaw321FromQuaterion(q);

				EXPECT_TRUE(matrix::isEqual(H_YAW, yaw.derivative)) << "yaw " << heading << std::endl << H_YAW(0) << " "
														<< H_YAW(1) << " "
														<< H_YAW(2) << " "
														<< H_YAW(3) << std::endl
														<< yaw.derivative(0) << " "
														<< yaw.derivative(1) << " "
														<< yaw.derivative(2) << " "
														<< yaw.derivative(3);

				EXPECT_NEAR(math::degrees(yaw.value), heading, 0.1f) << "roll: " << roll << ", pitch" << pitch << std::endl;
			}
		}
	}
}

TEST_F(EkfJacobianTest, jacobianWith312EulerAngles)
{


	for(float roll = -89.0f; roll <= 90.0f; roll+=45.0f)
	{
		for(float pitch = -89.0f; pitch <= 90.0f; pitch+=45.0f)
		{
			for(float heading = -180.0f; heading <= 180.0f; heading+=45.0f)
			{
				Quaternion<float> quat(Eulerf(math::radians(roll), math::radians(pitch), math::radians(heading)));

				const float q0 = quat(0);
				const float q1 = quat(1);
				const float q2 = quat(2);
				const float q3 = quat(3);

				// calculate observation jacobian when we are observing a rotation in a 312 sequence
				float t9 = q0*q3;
				float t10 = q1*q2;
				float t2 = t9-t10;
				float t3 = q0*q0;
				float t4 = q1*q1;
				float t5 = q2*q2;
				float t6 = q3*q3;
				float t7 = t3-t4+t5-t6;
				float t8 = t7*t7;
				if (t8 > 1e-6f) {
					t8 = 1.0f/t8;
				} else {
					t8 = 1.0f/t8;
				}
				float t11 = t2*t2;
				float t12 = t8*t11*4.0f;
				float t13 = t12+1.0f;
				float t14;

				if (fabsf(t13) > 1e-6f) {
					t14 = 1.0f/t13;
				} else {
					t14 = 1.0f/t13;
				}

				matrix::Vector<float, 4> H_YAW;
				H_YAW(0) = t8*t14*(q3*t3+q3*t4-q3*t5+q3*t6-q0*q1*q2*2.0f)*-2.0f;
				H_YAW(1) = t8*t14*(q2*t3+q2*t4+q2*t5-q2*t6-q0*q1*q3*2.0f)*-2.0f;
				H_YAW(2) = t8*t14*(-q1*t3+q1*t4+q1*t5+q1*t6-q0*q2*q3*2.0f)*2.0f;
				H_YAW(3) = t8*t14*(q0*t3-q0*t4+q0*t5+q0*t6-q1*q2*q3*2.0f)*2.0f;

				// alternative way
				using D4 = Dual<float, 4>;
				Quaternion<Dual<float, 4>> q(D4(q0, 0),D4(q1, 1),D4(q2, 2),D4(q3, 3));
				Dual<float, 4> yaw = yaw312FromQuaterion(q);

				EXPECT_TRUE(matrix::isEqual(H_YAW, yaw.derivative)) << "yaw " << heading << std::endl << H_YAW(0) << " "
														<< H_YAW(1) << " "
														<< H_YAW(2) << " "
														<< H_YAW(3) << std::endl
														<< yaw.derivative(0) << " "
														<< yaw.derivative(1) << " "
														<< yaw.derivative(2) << " "
														<< yaw.derivative(3);

				// This test will fail very likely as we swap the euler angle sequence.
				// EXPECT_NEAR(math::degrees(yaw.value), heading, 0.1f) << "roll: " << roll << ", pitch" << pitch << std::endl;
			}
		}
	}
}
