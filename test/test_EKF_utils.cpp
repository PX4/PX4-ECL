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

/**
 * @file test_EKF_utils.cpp
 *
 * @brief Unit tests for the miscellaneous EKF utilities
 */

#include <gtest/gtest.h>
#include <cmath>

#include "EKF/utils.hpp"

TEST(eclPowfTest, checkAccuracy)
{
	int exponent[7] = {-3,-2,-1,0,1,2,3};
	float x[2] = {-0.5f,0.5f};
	float expected_result[2][7] = {{-1.0f/(0.5f*0.5f*0.5f), 1.0f/(0.5f*0.5f), -1.0f/0.5f, 1.0f, -0.5f, 0.5f*0.5f, -0.5f*0.5f*0.5f},
				       { 1.0f/(0.5f*0.5f*0.5f), 1.0f/(0.5f*0.5f),  1.0f/0.5f, 1.0f,  0.5f, 0.5f*0.5f,  0.5f*0.5f*0.5f}};

	for (uint8_t exp_index = 0; exp_index < 7; exp_index++) {
		for (uint8_t x_index = 0; x_index < 2; x_index++) {
			const float test_result = ecl::powf(x[x_index],exponent[exp_index]);
			ASSERT_EQ(test_result, expected_result[x_index][exp_index]);
		}
	}
}
