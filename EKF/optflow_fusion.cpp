/****************************************************************************
 *
 *   Copyright (c) 2015 Estimation and Control Library (ECL). All rights reserved.
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
 * @file vel_pos_fusion.cpp
 * Function for fusing gps and baro measurements/
 *
 * @author Paul Riseborough <p_riseborough@live.com.au>
 * @author Siddharth Bharat Purohit <siddharthbharatpurohit@gmail.com>
 *
 */

#include "ekf.h"
#include <mathlib/mathlib.h>

void Ekf::fuseOptFlow()
{
    float terrainState = 0.0f;
    float gndclearance = _params.rng_gnd_clearance;
	// get latest estimated orientation
	float q0 = _state.quat_nominal(0);
	float q1 = _state.quat_nominal(1);
	float q2 = _state.quat_nominal(2);
	float q3 = _state.quat_nominal(3);

	// get latest velocity in earth frame
	float vn = _state.vel(0);
	float ve = _state.vel(1);
	float vd = _state.vel(2);
	float pd = _state.pos(2);
	Vector3f vel_body;

	float R_LOS = sq(_params.flow_noise);

	matrix::Dcm<float> earth_to_body(_state.quat_nominal);	//convert quat to DCM
	earth_to_body = earth_to_body.transpose();	// calculate earth to body rot mat

	vel_body = earth_to_body*_state.vel;	// rotate vel in earth frame to sensor/body frame
	Vector2f _flow_innov;

    // calculate innovations
	float range = -pd/earth_to_body(2,2); // absolute distance to the frame region in view 

	// range cannot be less than gnd clearance
	if(range < gndclearance) {
		range = gndclearance;
	}

	if (fabsf(range) < 1e-6f) {
		return;
	}

    _flow_innov(0) = vel_body(1)/range - _flow_sample_delayed.flowRadXYcomp(1); //innov in X direction
	_flow_innov(1) = -vel_body(0)/range + _flow_sample_delayed.flowRadXYcomp(0); //innov in Y direction

	// TODO: work out the process to take terrain changes into account
	// constrain height above ground to be above range measured on ground
    float heightAboveGndEst = math::max((terrainState - pd), gndclearance);
    float ptd = pd + heightAboveGndEst;

	// intermediate variable from algebraic optimisation
	float SH_LOS[7];
	memset(SH_LOS,0,sizeof(SH_LOS));
	SH_LOS[0] = sq(q0) - sq(q1) - sq(q2) + sq(q3);
	SH_LOS[1] = vn*(sq(q0) + sq(q1) - sq(q2) - sq(q3)) - vd*(2*q0*q2 - 2*q1*q3) + ve*(2*q0*q3 + 2*q1*q2);
	SH_LOS[2] = ve*(sq(q0) - sq(q1) + sq(q2) - sq(q3)) + vd*(2*q0*q1 + 2*q2*q3) - vn*(2*q0*q3 - 2*q1*q2);
	SH_LOS[3] = 1/(pd - ptd);
	SH_LOS[4] = vd*SH_LOS[0] - ve*(2*q0*q1 - 2*q2*q3) + vn*(2*q0*q2 + 2*q1*q3);
	SH_LOS[5] = 2*q0*q2 - 2*q1*q3;
	SH_LOS[6] = 2*q0*q1 + 2*q2*q3;

	// Intermediate variables used to calculate the Kalman gain matrices
	float SK_LOS[5];
	memset(SK_LOS,0,sizeof(SK_LOS));
	// this is 1/(innovation variance) for the X axis measurement
	SK_LOS[0] = 1/(R_LOS - (SH_LOS[0]*SH_LOS[3]*SH_LOS[4] - SH_LOS[2]*SH_LOS[3]*SH_LOS[6])*(P[8][0]*SH_LOS[0]*SH_LOS[2]*sq(SH_LOS[3]) - P[0][0]*(SH_LOS[0]*SH_LOS[3]*SH_LOS[4] - SH_LOS[2]*SH_LOS[3]*SH_LOS[6]) + P[3][0]*SH_LOS[0]*SH_LOS[3]*(2*q0*q3 - 2*q1*q2) + P[1][0]*SH_LOS[2]*SH_LOS[3]*SH_LOS[5] + P[2][0]*SH_LOS[0]*SH_LOS[1]*SH_LOS[3] - P[5][0]*SH_LOS[0]*SH_LOS[3]*SH_LOS[6] - P[4][0]*SH_LOS[0]*SH_LOS[3]*(sq(q0) - sq(q1) + sq(q2) - sq(q3))) + SH_LOS[2]*SH_LOS[3]*SH_LOS[5]*(P[8][1]*SH_LOS[0]*SH_LOS[2]*sq(SH_LOS[3]) - P[0][1]*(SH_LOS[0]*SH_LOS[3]*SH_LOS[4] - SH_LOS[2]*SH_LOS[3]*SH_LOS[6]) + P[3][1]*SH_LOS[0]*SH_LOS[3]*(2*q0*q3 - 2*q1*q2) + P[1][1]*SH_LOS[2]*SH_LOS[3]*SH_LOS[5] + P[2][1]*SH_LOS[0]*SH_LOS[1]*SH_LOS[3] - P[5][1]*SH_LOS[0]*SH_LOS[3]*SH_LOS[6] - P[4][1]*SH_LOS[0]*SH_LOS[3]*(sq(q0) - sq(q1) + sq(q2) - sq(q3))) + SH_LOS[0]*SH_LOS[1]*SH_LOS[3]*(P[8][2]*SH_LOS[0]*SH_LOS[2]*sq(SH_LOS[3]) - P[0][2]*(SH_LOS[0]*SH_LOS[3]*SH_LOS[4] - SH_LOS[2]*SH_LOS[3]*SH_LOS[6]) + P[3][2]*SH_LOS[0]*SH_LOS[3]*(2*q0*q3 - 2*q1*q2) + P[1][2]*SH_LOS[2]*SH_LOS[3]*SH_LOS[5] + P[2][2]*SH_LOS[0]*SH_LOS[1]*SH_LOS[3] - P[5][2]*SH_LOS[0]*SH_LOS[3]*SH_LOS[6] - P[4][2]*SH_LOS[0]*SH_LOS[3]*(sq(q0) - sq(q1) + sq(q2) - sq(q3))) - SH_LOS[0]*SH_LOS[3]*SH_LOS[6]*(P[8][5]*SH_LOS[0]*SH_LOS[2]*sq(SH_LOS[3]) - P[0][5]*(SH_LOS[0]*SH_LOS[3]*SH_LOS[4] - SH_LOS[2]*SH_LOS[3]*SH_LOS[6]) + P[3][5]*SH_LOS[0]*SH_LOS[3]*(2*q0*q3 - 2*q1*q2) + P[1][5]*SH_LOS[2]*SH_LOS[3]*SH_LOS[5] + P[2][5]*SH_LOS[0]*SH_LOS[1]*SH_LOS[3] - P[5][5]*SH_LOS[0]*SH_LOS[3]*SH_LOS[6] - P[4][5]*SH_LOS[0]*SH_LOS[3]*(sq(q0) - sq(q1) + sq(q2) - sq(q3))) - SH_LOS[0]*SH_LOS[3]*(sq(q0) - sq(q1) + sq(q2) - sq(q3))*(P[8][4]*SH_LOS[0]*SH_LOS[2]*sq(SH_LOS[3]) - P[0][4]*(SH_LOS[0]*SH_LOS[3]*SH_LOS[4] - SH_LOS[2]*SH_LOS[3]*SH_LOS[6]) + P[3][4]*SH_LOS[0]*SH_LOS[3]*(2*q0*q3 - 2*q1*q2) + P[1][4]*SH_LOS[2]*SH_LOS[3]*SH_LOS[5] + P[2][4]*SH_LOS[0]*SH_LOS[1]*SH_LOS[3] - P[5][4]*SH_LOS[0]*SH_LOS[3]*SH_LOS[6] - P[4][4]*SH_LOS[0]*SH_LOS[3]*(sq(q0) - sq(q1) + sq(q2) - sq(q3))) + SH_LOS[0]*SH_LOS[2]*sq(SH_LOS[3])*(P[8][8]*SH_LOS[0]*SH_LOS[2]*sq(SH_LOS[3]) - P[0][8]*(SH_LOS[0]*SH_LOS[3]*SH_LOS[4] - SH_LOS[2]*SH_LOS[3]*SH_LOS[6]) + P[3][8]*SH_LOS[0]*SH_LOS[3]*(2*q0*q3 - 2*q1*q2) + P[1][8]*SH_LOS[2]*SH_LOS[3]*SH_LOS[5] + P[2][8]*SH_LOS[0]*SH_LOS[1]*SH_LOS[3] - P[5][8]*SH_LOS[0]*SH_LOS[3]*SH_LOS[6] - P[4][8]*SH_LOS[0]*SH_LOS[3]*(sq(q0) - sq(q1) + sq(q2) - sq(q3))) + SH_LOS[0]*SH_LOS[3]*(2*q0*q3 - 2*q1*q2)*(P[8][3]*SH_LOS[0]*SH_LOS[2]*sq(SH_LOS[3]) - P[0][3]*(SH_LOS[0]*SH_LOS[3]*SH_LOS[4] - SH_LOS[2]*SH_LOS[3]*SH_LOS[6]) + P[3][3]*SH_LOS[0]*SH_LOS[3]*(2*q0*q3 - 2*q1*q2) + P[1][3]*SH_LOS[2]*SH_LOS[3]*SH_LOS[5] + P[2][3]*SH_LOS[0]*SH_LOS[1]*SH_LOS[3] - P[5][3]*SH_LOS[0]*SH_LOS[3]*SH_LOS[6] - P[4][3]*SH_LOS[0]*SH_LOS[3]*(sq(q0) - sq(q1) + sq(q2) - sq(q3))));
	// this is 1/(innovation variance) for the Y axis measurement
	SK_LOS[1] = 1/(R_LOS + (SH_LOS[0]*SH_LOS[3]*SH_LOS[4] + SH_LOS[1]*SH_LOS[3]*SH_LOS[5])*(P[1][1]*(SH_LOS[0]*SH_LOS[3]*SH_LOS[4] + SH_LOS[1]*SH_LOS[3]*SH_LOS[5]) + P[8][1]*SH_LOS[0]*SH_LOS[1]*sq(SH_LOS[3]) - P[4][1]*SH_LOS[0]*SH_LOS[3]*(2*q0*q3 + 2*q1*q2) + P[0][1]*SH_LOS[1]*SH_LOS[3]*SH_LOS[6] - P[2][1]*SH_LOS[0]*SH_LOS[2]*SH_LOS[3] + P[5][1]*SH_LOS[0]*SH_LOS[3]*SH_LOS[5] - P[3][1]*SH_LOS[0]*SH_LOS[3]*(sq(q0) + sq(q1) - sq(q2) - sq(q3))) + SH_LOS[1]*SH_LOS[3]*SH_LOS[6]*(P[1][0]*(SH_LOS[0]*SH_LOS[3]*SH_LOS[4] + SH_LOS[1]*SH_LOS[3]*SH_LOS[5]) + P[8][0]*SH_LOS[0]*SH_LOS[1]*sq(SH_LOS[3]) - P[4][0]*SH_LOS[0]*SH_LOS[3]*(2*q0*q3 + 2*q1*q2) + P[0][0]*SH_LOS[1]*SH_LOS[3]*SH_LOS[6] - P[2][0]*SH_LOS[0]*SH_LOS[2]*SH_LOS[3] + P[5][0]*SH_LOS[0]*SH_LOS[3]*SH_LOS[5] - P[3][0]*SH_LOS[0]*SH_LOS[3]*(sq(q0) + sq(q1) - sq(q2) - sq(q3))) - SH_LOS[0]*SH_LOS[2]*SH_LOS[3]*(P[1][2]*(SH_LOS[0]*SH_LOS[3]*SH_LOS[4] + SH_LOS[1]*SH_LOS[3]*SH_LOS[5]) + P[8][2]*SH_LOS[0]*SH_LOS[1]*sq(SH_LOS[3]) - P[4][2]*SH_LOS[0]*SH_LOS[3]*(2*q0*q3 + 2*q1*q2) + P[0][2]*SH_LOS[1]*SH_LOS[3]*SH_LOS[6] - P[2][2]*SH_LOS[0]*SH_LOS[2]*SH_LOS[3] + P[5][2]*SH_LOS[0]*SH_LOS[3]*SH_LOS[5] - P[3][2]*SH_LOS[0]*SH_LOS[3]*(sq(q0) + sq(q1) - sq(q2) - sq(q3))) + SH_LOS[0]*SH_LOS[3]*SH_LOS[5]*(P[1][5]*(SH_LOS[0]*SH_LOS[3]*SH_LOS[4] + SH_LOS[1]*SH_LOS[3]*SH_LOS[5]) + P[8][5]*SH_LOS[0]*SH_LOS[1]*sq(SH_LOS[3]) - P[4][5]*SH_LOS[0]*SH_LOS[3]*(2*q0*q3 + 2*q1*q2) + P[0][5]*SH_LOS[1]*SH_LOS[3]*SH_LOS[6] - P[2][5]*SH_LOS[0]*SH_LOS[2]*SH_LOS[3] + P[5][5]*SH_LOS[0]*SH_LOS[3]*SH_LOS[5] - P[3][5]*SH_LOS[0]*SH_LOS[3]*(sq(q0) + sq(q1) - sq(q2) - sq(q3))) - SH_LOS[0]*SH_LOS[3]*(sq(q0) + sq(q1) - sq(q2) - sq(q3))*(P[1][3]*(SH_LOS[0]*SH_LOS[3]*SH_LOS[4] + SH_LOS[1]*SH_LOS[3]*SH_LOS[5]) + P[8][3]*SH_LOS[0]*SH_LOS[1]*sq(SH_LOS[3]) - P[4][3]*SH_LOS[0]*SH_LOS[3]*(2*q0*q3 + 2*q1*q2) + P[0][3]*SH_LOS[1]*SH_LOS[3]*SH_LOS[6] - P[2][3]*SH_LOS[0]*SH_LOS[2]*SH_LOS[3] + P[5][3]*SH_LOS[0]*SH_LOS[3]*SH_LOS[5] - P[3][3]*SH_LOS[0]*SH_LOS[3]*(sq(q0) + sq(q1) - sq(q2) - sq(q3))) + SH_LOS[0]*SH_LOS[1]*sq(SH_LOS[3])*(P[1][8]*(SH_LOS[0]*SH_LOS[3]*SH_LOS[4] + SH_LOS[1]*SH_LOS[3]*SH_LOS[5]) + P[8][8]*SH_LOS[0]*SH_LOS[1]*sq(SH_LOS[3]) - P[4][8]*SH_LOS[0]*SH_LOS[3]*(2*q0*q3 + 2*q1*q2) + P[0][8]*SH_LOS[1]*SH_LOS[3]*SH_LOS[6] - P[2][8]*SH_LOS[0]*SH_LOS[2]*SH_LOS[3] + P[5][8]*SH_LOS[0]*SH_LOS[3]*SH_LOS[5] - P[3][8]*SH_LOS[0]*SH_LOS[3]*(sq(q0) + sq(q1) - sq(q2) - sq(q3))) - SH_LOS[0]*SH_LOS[3]*(2*q0*q3 + 2*q1*q2)*(P[1][4]*(SH_LOS[0]*SH_LOS[3]*SH_LOS[4] + SH_LOS[1]*SH_LOS[3]*SH_LOS[5]) + P[8][4]*SH_LOS[0]*SH_LOS[1]*sq(SH_LOS[3]) - P[4][4]*SH_LOS[0]*SH_LOS[3]*(2*q0*q3 + 2*q1*q2) + P[0][4]*SH_LOS[1]*SH_LOS[3]*SH_LOS[6] - P[2][4]*SH_LOS[0]*SH_LOS[2]*SH_LOS[3] + P[5][4]*SH_LOS[0]*SH_LOS[3]*SH_LOS[5] - P[3][4]*SH_LOS[0]*SH_LOS[3]*(sq(q0) + sq(q1) - sq(q2) - sq(q3))));
	SK_LOS[2] = sq(q0) + sq(q1) - sq(q2) - sq(q3);
	SK_LOS[3] = sq(q0) - sq(q1) + sq(q2) - sq(q3);
	SK_LOS[4] = SH_LOS[3];

    // run innovation variance checks
	float optflow_test_ratio_x = sq(_flow_innov(0)) / (sq(math::max(_params.flow_innov_gate, 3.0f)) * (1/SK_LOS[0]));
	float optflow_test_ratio_y = sq(_flow_innov(1)) / (sq(math::max(_params.flow_innov_gate, 3.0f)) * (1/SK_LOS[1]));
	if(optflow_test_ratio_x > 1.0f || optflow_test_ratio_y > 1.0f) {
        return;
	}

	for(uint8_t index = 0; index < 2; index++ ) {

		float Kfusion[24];
		float H_LOS[24];

        //initialise observation matrix
		memset(H_LOS,0,sizeof(H_LOS));

		if(index == 0) {
			// Calculate the observation jacobians for the LOS rate about the X body axis
			H_LOS[0] = SH_LOS[2]*SH_LOS[3]*SH_LOS[6] - SH_LOS[0]*SH_LOS[3]*SH_LOS[4];
			H_LOS[1] = SH_LOS[2]*SH_LOS[3]*SH_LOS[5];
			H_LOS[2] = SH_LOS[0]*SH_LOS[1]*SH_LOS[3];
			H_LOS[3] = SH_LOS[0]*SH_LOS[3]*(2*q0*q3 - 2*q1*q2);
			H_LOS[4] = -SH_LOS[0]*SH_LOS[3]*(sq(q0) - sq(q1) + sq(q2) - sq(q3));
			H_LOS[5] = -SH_LOS[0]*SH_LOS[3]*SH_LOS[6];
			H_LOS[8] = SH_LOS[0]*SH_LOS[2]*sq(SH_LOS[3]);

			// Calculate the Kalman gain matrix for the X axis measurement
			Kfusion[0] = SK_LOS[0]*(P[0][8]*SH_LOS[0]*SH_LOS[2]*sq(SK_LOS[4]) - P[0][0]*(SH_LOS[0]*SH_LOS[4]*SK_LOS[4] - SH_LOS[2]*SH_LOS[6]*SK_LOS[4]) + P[0][3]*SH_LOS[0]*SK_LOS[4]*(2*q0*q3 - 2*q1*q2) + P[0][2]*SH_LOS[0]*SH_LOS[1]*SK_LOS[4] + P[0][1]*SH_LOS[2]*SH_LOS[5]*SK_LOS[4] - P[0][5]*SH_LOS[0]*SH_LOS[6]*SK_LOS[4] - P[0][4]*SH_LOS[0]*SK_LOS[3]*SK_LOS[4]);
			Kfusion[1] = SK_LOS[0]*(P[1][8]*SH_LOS[0]*SH_LOS[2]*sq(SK_LOS[4]) - P[1][0]*(SH_LOS[0]*SH_LOS[4]*SK_LOS[4] - SH_LOS[2]*SH_LOS[6]*SK_LOS[4]) + P[1][3]*SH_LOS[0]*SK_LOS[4]*(2*q0*q3 - 2*q1*q2) + P[1][2]*SH_LOS[0]*SH_LOS[1]*SK_LOS[4] + P[1][1]*SH_LOS[2]*SH_LOS[5]*SK_LOS[4] - P[1][5]*SH_LOS[0]*SH_LOS[6]*SK_LOS[4] - P[1][4]*SH_LOS[0]*SK_LOS[3]*SK_LOS[4]);
			Kfusion[2] = SK_LOS[0]*(P[2][8]*SH_LOS[0]*SH_LOS[2]*sq(SK_LOS[4]) - P[2][0]*(SH_LOS[0]*SH_LOS[4]*SK_LOS[4] - SH_LOS[2]*SH_LOS[6]*SK_LOS[4]) + P[2][3]*SH_LOS[0]*SK_LOS[4]*(2*q0*q3 - 2*q1*q2) + P[2][2]*SH_LOS[0]*SH_LOS[1]*SK_LOS[4] + P[2][1]*SH_LOS[2]*SH_LOS[5]*SK_LOS[4] - P[2][5]*SH_LOS[0]*SH_LOS[6]*SK_LOS[4] - P[2][4]*SH_LOS[0]*SK_LOS[3]*SK_LOS[4]);
			Kfusion[3] = SK_LOS[0]*(P[3][8]*SH_LOS[0]*SH_LOS[2]*sq(SK_LOS[4]) - P[3][0]*(SH_LOS[0]*SH_LOS[4]*SK_LOS[4] - SH_LOS[2]*SH_LOS[6]*SK_LOS[4]) + P[3][3]*SH_LOS[0]*SK_LOS[4]*(2*q0*q3 - 2*q1*q2) + P[3][2]*SH_LOS[0]*SH_LOS[1]*SK_LOS[4] + P[3][1]*SH_LOS[2]*SH_LOS[5]*SK_LOS[4] - P[3][5]*SH_LOS[0]*SH_LOS[6]*SK_LOS[4] - P[3][4]*SH_LOS[0]*SK_LOS[3]*SK_LOS[4]);
			Kfusion[4] = SK_LOS[0]*(P[4][8]*SH_LOS[0]*SH_LOS[2]*sq(SK_LOS[4]) - P[4][0]*(SH_LOS[0]*SH_LOS[4]*SK_LOS[4] - SH_LOS[2]*SH_LOS[6]*SK_LOS[4]) + P[4][3]*SH_LOS[0]*SK_LOS[4]*(2*q0*q3 - 2*q1*q2) + P[4][2]*SH_LOS[0]*SH_LOS[1]*SK_LOS[4] + P[4][1]*SH_LOS[2]*SH_LOS[5]*SK_LOS[4] - P[4][5]*SH_LOS[0]*SH_LOS[6]*SK_LOS[4] - P[4][4]*SH_LOS[0]*SK_LOS[3]*SK_LOS[4]);
			Kfusion[5] = SK_LOS[0]*(P[5][8]*SH_LOS[0]*SH_LOS[2]*sq(SK_LOS[4]) - P[5][0]*(SH_LOS[0]*SH_LOS[4]*SK_LOS[4] - SH_LOS[2]*SH_LOS[6]*SK_LOS[4]) + P[5][3]*SH_LOS[0]*SK_LOS[4]*(2*q0*q3 - 2*q1*q2) + P[5][2]*SH_LOS[0]*SH_LOS[1]*SK_LOS[4] + P[5][1]*SH_LOS[2]*SH_LOS[5]*SK_LOS[4] - P[5][5]*SH_LOS[0]*SH_LOS[6]*SK_LOS[4] - P[5][4]*SH_LOS[0]*SK_LOS[3]*SK_LOS[4]);
			Kfusion[6] = SK_LOS[0]*(P[6][8]*SH_LOS[0]*SH_LOS[2]*sq(SK_LOS[4]) - P[6][0]*(SH_LOS[0]*SH_LOS[4]*SK_LOS[4] - SH_LOS[2]*SH_LOS[6]*SK_LOS[4]) + P[6][3]*SH_LOS[0]*SK_LOS[4]*(2*q0*q3 - 2*q1*q2) + P[6][2]*SH_LOS[0]*SH_LOS[1]*SK_LOS[4] + P[6][1]*SH_LOS[2]*SH_LOS[5]*SK_LOS[4] - P[6][5]*SH_LOS[0]*SH_LOS[6]*SK_LOS[4] - P[6][4]*SH_LOS[0]*SK_LOS[3]*SK_LOS[4]);
			Kfusion[7] = SK_LOS[0]*(P[7][8]*SH_LOS[0]*SH_LOS[2]*sq(SK_LOS[4]) - P[7][0]*(SH_LOS[0]*SH_LOS[4]*SK_LOS[4] - SH_LOS[2]*SH_LOS[6]*SK_LOS[4]) + P[7][3]*SH_LOS[0]*SK_LOS[4]*(2*q0*q3 - 2*q1*q2) + P[7][2]*SH_LOS[0]*SH_LOS[1]*SK_LOS[4] + P[7][1]*SH_LOS[2]*SH_LOS[5]*SK_LOS[4] - P[7][5]*SH_LOS[0]*SH_LOS[6]*SK_LOS[4] - P[7][4]*SH_LOS[0]*SK_LOS[3]*SK_LOS[4]);
			Kfusion[8] = SK_LOS[0]*(P[8][8]*SH_LOS[0]*SH_LOS[2]*sq(SK_LOS[4]) - P[8][0]*(SH_LOS[0]*SH_LOS[4]*SK_LOS[4] - SH_LOS[2]*SH_LOS[6]*SK_LOS[4]) + P[8][3]*SH_LOS[0]*SK_LOS[4]*(2*q0*q3 - 2*q1*q2) + P[8][2]*SH_LOS[0]*SH_LOS[1]*SK_LOS[4] + P[8][1]*SH_LOS[2]*SH_LOS[5]*SK_LOS[4] - P[8][5]*SH_LOS[0]*SH_LOS[6]*SK_LOS[4] - P[8][4]*SH_LOS[0]*SK_LOS[3]*SK_LOS[4]);
			Kfusion[9] = SK_LOS[0]*(P[9][8]*SH_LOS[0]*SH_LOS[2]*sq(SK_LOS[4]) - P[9][0]*(SH_LOS[0]*SH_LOS[4]*SK_LOS[4] - SH_LOS[2]*SH_LOS[6]*SK_LOS[4]) + P[9][3]*SH_LOS[0]*SK_LOS[4]*(2*q0*q3 - 2*q1*q2) + P[9][2]*SH_LOS[0]*SH_LOS[1]*SK_LOS[4] + P[9][1]*SH_LOS[2]*SH_LOS[5]*SK_LOS[4] - P[9][5]*SH_LOS[0]*SH_LOS[6]*SK_LOS[4] - P[9][4]*SH_LOS[0]*SK_LOS[3]*SK_LOS[4]);
			Kfusion[10] = SK_LOS[0]*(P[10][8]*SH_LOS[0]*SH_LOS[2]*sq(SK_LOS[4]) - P[10][0]*(SH_LOS[0]*SH_LOS[4]*SK_LOS[4] - SH_LOS[2]*SH_LOS[6]*SK_LOS[4]) + P[10][3]*SH_LOS[0]*SK_LOS[4]*(2*q0*q3 - 2*q1*q2) + P[10][2]*SH_LOS[0]*SH_LOS[1]*SK_LOS[4] + P[10][1]*SH_LOS[2]*SH_LOS[5]*SK_LOS[4] - P[10][5]*SH_LOS[0]*SH_LOS[6]*SK_LOS[4] - P[10][4]*SH_LOS[0]*SK_LOS[3]*SK_LOS[4]);
			Kfusion[11] = SK_LOS[0]*(P[11][8]*SH_LOS[0]*SH_LOS[2]*sq(SK_LOS[4]) - P[11][0]*(SH_LOS[0]*SH_LOS[4]*SK_LOS[4] - SH_LOS[2]*SH_LOS[6]*SK_LOS[4]) + P[11][3]*SH_LOS[0]*SK_LOS[4]*(2*q0*q3 - 2*q1*q2) + P[11][2]*SH_LOS[0]*SH_LOS[1]*SK_LOS[4] + P[11][1]*SH_LOS[2]*SH_LOS[5]*SK_LOS[4] - P[11][5]*SH_LOS[0]*SH_LOS[6]*SK_LOS[4] - P[11][4]*SH_LOS[0]*SK_LOS[3]*SK_LOS[4]);
			Kfusion[12] = SK_LOS[0]*(P[12][8]*SH_LOS[0]*SH_LOS[2]*sq(SK_LOS[4]) - P[12][0]*(SH_LOS[0]*SH_LOS[4]*SK_LOS[4] - SH_LOS[2]*SH_LOS[6]*SK_LOS[4]) + P[12][3]*SH_LOS[0]*SK_LOS[4]*(2*q0*q3 - 2*q1*q2) + P[12][2]*SH_LOS[0]*SH_LOS[1]*SK_LOS[4] + P[12][1]*SH_LOS[2]*SH_LOS[5]*SK_LOS[4] - P[12][5]*SH_LOS[0]*SH_LOS[6]*SK_LOS[4] - P[12][4]*SH_LOS[0]*SK_LOS[3]*SK_LOS[4]);
			Kfusion[13] = SK_LOS[0]*(P[13][8]*SH_LOS[0]*SH_LOS[2]*sq(SK_LOS[4]) - P[13][0]*(SH_LOS[0]*SH_LOS[4]*SK_LOS[4] - SH_LOS[2]*SH_LOS[6]*SK_LOS[4]) + P[13][3]*SH_LOS[0]*SK_LOS[4]*(2*q0*q3 - 2*q1*q2) + P[13][2]*SH_LOS[0]*SH_LOS[1]*SK_LOS[4] + P[13][1]*SH_LOS[2]*SH_LOS[5]*SK_LOS[4] - P[13][5]*SH_LOS[0]*SH_LOS[6]*SK_LOS[4] - P[13][4]*SH_LOS[0]*SK_LOS[3]*SK_LOS[4]);
			Kfusion[14] = SK_LOS[0]*(P[14][8]*SH_LOS[0]*SH_LOS[2]*sq(SK_LOS[4]) - P[14][0]*(SH_LOS[0]*SH_LOS[4]*SK_LOS[4] - SH_LOS[2]*SH_LOS[6]*SK_LOS[4]) + P[14][3]*SH_LOS[0]*SK_LOS[4]*(2*q0*q3 - 2*q1*q2) + P[14][2]*SH_LOS[0]*SH_LOS[1]*SK_LOS[4] + P[14][1]*SH_LOS[2]*SH_LOS[5]*SK_LOS[4] - P[14][5]*SH_LOS[0]*SH_LOS[6]*SK_LOS[4] - P[14][4]*SH_LOS[0]*SK_LOS[3]*SK_LOS[4]);
			Kfusion[15] = SK_LOS[0]*(P[15][8]*SH_LOS[0]*SH_LOS[2]*sq(SK_LOS[4]) - P[15][0]*(SH_LOS[0]*SH_LOS[4]*SK_LOS[4] - SH_LOS[2]*SH_LOS[6]*SK_LOS[4]) + P[15][3]*SH_LOS[0]*SK_LOS[4]*(2*q0*q3 - 2*q1*q2) + P[15][2]*SH_LOS[0]*SH_LOS[1]*SK_LOS[4] + P[15][1]*SH_LOS[2]*SH_LOS[5]*SK_LOS[4] - P[15][5]*SH_LOS[0]*SH_LOS[6]*SK_LOS[4] - P[15][4]*SH_LOS[0]*SK_LOS[3]*SK_LOS[4]);
			Kfusion[16] = SK_LOS[0]*(P[16][8]*SH_LOS[0]*SH_LOS[2]*sq(SK_LOS[4]) - P[16][0]*(SH_LOS[0]*SH_LOS[4]*SK_LOS[4] - SH_LOS[2]*SH_LOS[6]*SK_LOS[4]) + P[16][3]*SH_LOS[0]*SK_LOS[4]*(2*q0*q3 - 2*q1*q2) + P[16][2]*SH_LOS[0]*SH_LOS[1]*SK_LOS[4] + P[16][1]*SH_LOS[2]*SH_LOS[5]*SK_LOS[4] - P[16][5]*SH_LOS[0]*SH_LOS[6]*SK_LOS[4] - P[16][4]*SH_LOS[0]*SK_LOS[3]*SK_LOS[4]);
			Kfusion[17] = SK_LOS[0]*(P[17][8]*SH_LOS[0]*SH_LOS[2]*sq(SK_LOS[4]) - P[17][0]*(SH_LOS[0]*SH_LOS[4]*SK_LOS[4] - SH_LOS[2]*SH_LOS[6]*SK_LOS[4]) + P[17][3]*SH_LOS[0]*SK_LOS[4]*(2*q0*q3 - 2*q1*q2) + P[17][2]*SH_LOS[0]*SH_LOS[1]*SK_LOS[4] + P[17][1]*SH_LOS[2]*SH_LOS[5]*SK_LOS[4] - P[17][5]*SH_LOS[0]*SH_LOS[6]*SK_LOS[4] - P[17][4]*SH_LOS[0]*SK_LOS[3]*SK_LOS[4]);
			Kfusion[18] = SK_LOS[0]*(P[18][8]*SH_LOS[0]*SH_LOS[2]*sq(SK_LOS[4]) - P[18][0]*(SH_LOS[0]*SH_LOS[4]*SK_LOS[4] - SH_LOS[2]*SH_LOS[6]*SK_LOS[4]) + P[18][3]*SH_LOS[0]*SK_LOS[4]*(2*q0*q3 - 2*q1*q2) + P[18][2]*SH_LOS[0]*SH_LOS[1]*SK_LOS[4] + P[18][1]*SH_LOS[2]*SH_LOS[5]*SK_LOS[4] - P[18][5]*SH_LOS[0]*SH_LOS[6]*SK_LOS[4] - P[18][4]*SH_LOS[0]*SK_LOS[3]*SK_LOS[4]);
			Kfusion[19] = SK_LOS[0]*(P[19][8]*SH_LOS[0]*SH_LOS[2]*sq(SK_LOS[4]) - P[19][0]*(SH_LOS[0]*SH_LOS[4]*SK_LOS[4] - SH_LOS[2]*SH_LOS[6]*SK_LOS[4]) + P[19][3]*SH_LOS[0]*SK_LOS[4]*(2*q0*q3 - 2*q1*q2) + P[19][2]*SH_LOS[0]*SH_LOS[1]*SK_LOS[4] + P[19][1]*SH_LOS[2]*SH_LOS[5]*SK_LOS[4] - P[19][5]*SH_LOS[0]*SH_LOS[6]*SK_LOS[4] - P[19][4]*SH_LOS[0]*SK_LOS[3]*SK_LOS[4]);
			Kfusion[20] = SK_LOS[0]*(P[20][8]*SH_LOS[0]*SH_LOS[2]*sq(SK_LOS[4]) - P[20][0]*(SH_LOS[0]*SH_LOS[4]*SK_LOS[4] - SH_LOS[2]*SH_LOS[6]*SK_LOS[4]) + P[20][3]*SH_LOS[0]*SK_LOS[4]*(2*q0*q3 - 2*q1*q2) + P[20][2]*SH_LOS[0]*SH_LOS[1]*SK_LOS[4] + P[20][1]*SH_LOS[2]*SH_LOS[5]*SK_LOS[4] - P[20][5]*SH_LOS[0]*SH_LOS[6]*SK_LOS[4] - P[20][4]*SH_LOS[0]*SK_LOS[3]*SK_LOS[4]);
			Kfusion[21] = SK_LOS[0]*(P[21][8]*SH_LOS[0]*SH_LOS[2]*sq(SK_LOS[4]) - P[21][0]*(SH_LOS[0]*SH_LOS[4]*SK_LOS[4] - SH_LOS[2]*SH_LOS[6]*SK_LOS[4]) + P[21][3]*SH_LOS[0]*SK_LOS[4]*(2*q0*q3 - 2*q1*q2) + P[21][2]*SH_LOS[0]*SH_LOS[1]*SK_LOS[4] + P[21][1]*SH_LOS[2]*SH_LOS[5]*SK_LOS[4] - P[21][5]*SH_LOS[0]*SH_LOS[6]*SK_LOS[4] - P[21][4]*SH_LOS[0]*SK_LOS[3]*SK_LOS[4]);
			Kfusion[22] = SK_LOS[0]*(P[22][8]*SH_LOS[0]*SH_LOS[2]*sq(SK_LOS[4]) - P[22][0]*(SH_LOS[0]*SH_LOS[4]*SK_LOS[4] - SH_LOS[2]*SH_LOS[6]*SK_LOS[4]) + P[22][3]*SH_LOS[0]*SK_LOS[4]*(2*q0*q3 - 2*q1*q2) + P[22][2]*SH_LOS[0]*SH_LOS[1]*SK_LOS[4] + P[22][1]*SH_LOS[2]*SH_LOS[5]*SK_LOS[4] - P[22][5]*SH_LOS[0]*SH_LOS[6]*SK_LOS[4] - P[22][4]*SH_LOS[0]*SK_LOS[3]*SK_LOS[4]);
			Kfusion[23] = SK_LOS[0]*(P[23][8]*SH_LOS[0]*SH_LOS[2]*sq(SK_LOS[4]) - P[23][0]*(SH_LOS[0]*SH_LOS[4]*SK_LOS[4] - SH_LOS[2]*SH_LOS[6]*SK_LOS[4]) + P[23][3]*SH_LOS[0]*SK_LOS[4]*(2*q0*q3 - 2*q1*q2) + P[23][2]*SH_LOS[0]*SH_LOS[1]*SK_LOS[4] + P[23][1]*SH_LOS[2]*SH_LOS[5]*SK_LOS[4] - P[23][5]*SH_LOS[0]*SH_LOS[6]*SK_LOS[4] - P[23][4]*SH_LOS[0]*SK_LOS[3]*SK_LOS[4]);
		} else {
			// Calculate the observation jacobians for the LOS rate about the Y body axis
			H_LOS[0] = -SH_LOS[1]*SH_LOS[3]*SH_LOS[6];
			H_LOS[1] = - SH_LOS[0]*SH_LOS[3]*SH_LOS[4] - SH_LOS[1]*SH_LOS[3]*SH_LOS[5];
			H_LOS[2] = SH_LOS[0]*SH_LOS[2]*SH_LOS[3];
			H_LOS[3] = SH_LOS[0]*SH_LOS[3]*(sq(q0) + sq(q1) - sq(q2) - sq(q3));
			H_LOS[4] = SH_LOS[0]*SH_LOS[3]*(2.0f*q0*q3 + 2.0f*q1*q2);
			H_LOS[5] = -SH_LOS[0]*SH_LOS[3]*SH_LOS[5];
			H_LOS[8] = -SH_LOS[0]*SH_LOS[1]*sq(SH_LOS[3]);

			// Calculate the Kalman gain matrix for the Y axis measurement
			Kfusion[0] = -SK_LOS[1]*(P[0][1]*(SH_LOS[0]*SH_LOS[4]*SK_LOS[4] + SH_LOS[1]*SH_LOS[5]*SK_LOS[4]) + P[0][8]*SH_LOS[0]*SH_LOS[1]*sq(SK_LOS[4]) - P[0][4]*SH_LOS[0]*SK_LOS[4]*(2*q0*q3 + 2*q1*q2) - P[0][2]*SH_LOS[0]*SH_LOS[2]*SK_LOS[4] + P[0][0]*SH_LOS[1]*SH_LOS[6]*SK_LOS[4] + P[0][5]*SH_LOS[0]*SH_LOS[5]*SK_LOS[4] - P[0][3]*SH_LOS[0]*SK_LOS[2]*SK_LOS[4]);
			Kfusion[1] = -SK_LOS[1]*(P[1][1]*(SH_LOS[0]*SH_LOS[4]*SK_LOS[4] + SH_LOS[1]*SH_LOS[5]*SK_LOS[4]) + P[1][8]*SH_LOS[0]*SH_LOS[1]*sq(SK_LOS[4]) - P[1][4]*SH_LOS[0]*SK_LOS[4]*(2*q0*q3 + 2*q1*q2) - P[1][2]*SH_LOS[0]*SH_LOS[2]*SK_LOS[4] + P[1][0]*SH_LOS[1]*SH_LOS[6]*SK_LOS[4] + P[1][5]*SH_LOS[0]*SH_LOS[5]*SK_LOS[4] - P[1][3]*SH_LOS[0]*SK_LOS[2]*SK_LOS[4]);
			Kfusion[2] = -SK_LOS[1]*(P[2][1]*(SH_LOS[0]*SH_LOS[4]*SK_LOS[4] + SH_LOS[1]*SH_LOS[5]*SK_LOS[4]) + P[2][8]*SH_LOS[0]*SH_LOS[1]*sq(SK_LOS[4]) - P[2][4]*SH_LOS[0]*SK_LOS[4]*(2*q0*q3 + 2*q1*q2) - P[2][2]*SH_LOS[0]*SH_LOS[2]*SK_LOS[4] + P[2][0]*SH_LOS[1]*SH_LOS[6]*SK_LOS[4] + P[2][5]*SH_LOS[0]*SH_LOS[5]*SK_LOS[4] - P[2][3]*SH_LOS[0]*SK_LOS[2]*SK_LOS[4]);
			Kfusion[3] = -SK_LOS[1]*(P[3][1]*(SH_LOS[0]*SH_LOS[4]*SK_LOS[4] + SH_LOS[1]*SH_LOS[5]*SK_LOS[4]) + P[3][8]*SH_LOS[0]*SH_LOS[1]*sq(SK_LOS[4]) - P[3][4]*SH_LOS[0]*SK_LOS[4]*(2*q0*q3 + 2*q1*q2) - P[3][2]*SH_LOS[0]*SH_LOS[2]*SK_LOS[4] + P[3][0]*SH_LOS[1]*SH_LOS[6]*SK_LOS[4] + P[3][5]*SH_LOS[0]*SH_LOS[5]*SK_LOS[4] - P[3][3]*SH_LOS[0]*SK_LOS[2]*SK_LOS[4]);
			Kfusion[4] = -SK_LOS[1]*(P[4][1]*(SH_LOS[0]*SH_LOS[4]*SK_LOS[4] + SH_LOS[1]*SH_LOS[5]*SK_LOS[4]) + P[4][8]*SH_LOS[0]*SH_LOS[1]*sq(SK_LOS[4]) - P[4][4]*SH_LOS[0]*SK_LOS[4]*(2*q0*q3 + 2*q1*q2) - P[4][2]*SH_LOS[0]*SH_LOS[2]*SK_LOS[4] + P[4][0]*SH_LOS[1]*SH_LOS[6]*SK_LOS[4] + P[4][5]*SH_LOS[0]*SH_LOS[5]*SK_LOS[4] - P[4][3]*SH_LOS[0]*SK_LOS[2]*SK_LOS[4]);
			Kfusion[5] = -SK_LOS[1]*(P[5][1]*(SH_LOS[0]*SH_LOS[4]*SK_LOS[4] + SH_LOS[1]*SH_LOS[5]*SK_LOS[4]) + P[5][8]*SH_LOS[0]*SH_LOS[1]*sq(SK_LOS[4]) - P[5][4]*SH_LOS[0]*SK_LOS[4]*(2*q0*q3 + 2*q1*q2) - P[5][2]*SH_LOS[0]*SH_LOS[2]*SK_LOS[4] + P[5][0]*SH_LOS[1]*SH_LOS[6]*SK_LOS[4] + P[5][5]*SH_LOS[0]*SH_LOS[5]*SK_LOS[4] - P[5][3]*SH_LOS[0]*SK_LOS[2]*SK_LOS[4]);
			Kfusion[6] = -SK_LOS[1]*(P[6][1]*(SH_LOS[0]*SH_LOS[4]*SK_LOS[4] + SH_LOS[1]*SH_LOS[5]*SK_LOS[4]) + P[6][8]*SH_LOS[0]*SH_LOS[1]*sq(SK_LOS[4]) - P[6][4]*SH_LOS[0]*SK_LOS[4]*(2*q0*q3 + 2*q1*q2) - P[6][2]*SH_LOS[0]*SH_LOS[2]*SK_LOS[4] + P[6][0]*SH_LOS[1]*SH_LOS[6]*SK_LOS[4] + P[6][5]*SH_LOS[0]*SH_LOS[5]*SK_LOS[4] - P[6][3]*SH_LOS[0]*SK_LOS[2]*SK_LOS[4]);
			Kfusion[7] = -SK_LOS[1]*(P[7][1]*(SH_LOS[0]*SH_LOS[4]*SK_LOS[4] + SH_LOS[1]*SH_LOS[5]*SK_LOS[4]) + P[7][8]*SH_LOS[0]*SH_LOS[1]*sq(SK_LOS[4]) - P[7][4]*SH_LOS[0]*SK_LOS[4]*(2*q0*q3 + 2*q1*q2) - P[7][2]*SH_LOS[0]*SH_LOS[2]*SK_LOS[4] + P[7][0]*SH_LOS[1]*SH_LOS[6]*SK_LOS[4] + P[7][5]*SH_LOS[0]*SH_LOS[5]*SK_LOS[4] - P[7][3]*SH_LOS[0]*SK_LOS[2]*SK_LOS[4]);
			Kfusion[8] = -SK_LOS[1]*(P[8][1]*(SH_LOS[0]*SH_LOS[4]*SK_LOS[4] + SH_LOS[1]*SH_LOS[5]*SK_LOS[4]) + P[8][8]*SH_LOS[0]*SH_LOS[1]*sq(SK_LOS[4]) - P[8][4]*SH_LOS[0]*SK_LOS[4]*(2*q0*q3 + 2*q1*q2) - P[8][2]*SH_LOS[0]*SH_LOS[2]*SK_LOS[4] + P[8][0]*SH_LOS[1]*SH_LOS[6]*SK_LOS[4] + P[8][5]*SH_LOS[0]*SH_LOS[5]*SK_LOS[4] - P[8][3]*SH_LOS[0]*SK_LOS[2]*SK_LOS[4]);
			Kfusion[9] = -SK_LOS[1]*(P[9][1]*(SH_LOS[0]*SH_LOS[4]*SK_LOS[4] + SH_LOS[1]*SH_LOS[5]*SK_LOS[4]) + P[9][8]*SH_LOS[0]*SH_LOS[1]*sq(SK_LOS[4]) - P[9][4]*SH_LOS[0]*SK_LOS[4]*(2*q0*q3 + 2*q1*q2) - P[9][2]*SH_LOS[0]*SH_LOS[2]*SK_LOS[4] + P[9][0]*SH_LOS[1]*SH_LOS[6]*SK_LOS[4] + P[9][5]*SH_LOS[0]*SH_LOS[5]*SK_LOS[4] - P[9][3]*SH_LOS[0]*SK_LOS[2]*SK_LOS[4]);
			Kfusion[10] = -SK_LOS[1]*(P[10][1]*(SH_LOS[0]*SH_LOS[4]*SK_LOS[4] + SH_LOS[1]*SH_LOS[5]*SK_LOS[4]) + P[10][8]*SH_LOS[0]*SH_LOS[1]*sq(SK_LOS[4]) - P[10][4]*SH_LOS[0]*SK_LOS[4]*(2*q0*q3 + 2*q1*q2) - P[10][2]*SH_LOS[0]*SH_LOS[2]*SK_LOS[4] + P[10][0]*SH_LOS[1]*SH_LOS[6]*SK_LOS[4] + P[10][5]*SH_LOS[0]*SH_LOS[5]*SK_LOS[4] - P[10][3]*SH_LOS[0]*SK_LOS[2]*SK_LOS[4]);
			Kfusion[11] = -SK_LOS[1]*(P[11][1]*(SH_LOS[0]*SH_LOS[4]*SK_LOS[4] + SH_LOS[1]*SH_LOS[5]*SK_LOS[4]) + P[11][8]*SH_LOS[0]*SH_LOS[1]*sq(SK_LOS[4]) - P[11][4]*SH_LOS[0]*SK_LOS[4]*(2*q0*q3 + 2*q1*q2) - P[11][2]*SH_LOS[0]*SH_LOS[2]*SK_LOS[4] + P[11][0]*SH_LOS[1]*SH_LOS[6]*SK_LOS[4] + P[11][5]*SH_LOS[0]*SH_LOS[5]*SK_LOS[4] - P[11][3]*SH_LOS[0]*SK_LOS[2]*SK_LOS[4]);
			Kfusion[12] = -SK_LOS[1]*(P[12][1]*(SH_LOS[0]*SH_LOS[4]*SK_LOS[4] + SH_LOS[1]*SH_LOS[5]*SK_LOS[4]) + P[12][8]*SH_LOS[0]*SH_LOS[1]*sq(SK_LOS[4]) - P[12][4]*SH_LOS[0]*SK_LOS[4]*(2*q0*q3 + 2*q1*q2) - P[12][2]*SH_LOS[0]*SH_LOS[2]*SK_LOS[4] + P[12][0]*SH_LOS[1]*SH_LOS[6]*SK_LOS[4] + P[12][5]*SH_LOS[0]*SH_LOS[5]*SK_LOS[4] - P[12][3]*SH_LOS[0]*SK_LOS[2]*SK_LOS[4]);
			Kfusion[13] = -SK_LOS[1]*(P[13][1]*(SH_LOS[0]*SH_LOS[4]*SK_LOS[4] + SH_LOS[1]*SH_LOS[5]*SK_LOS[4]) + P[13][8]*SH_LOS[0]*SH_LOS[1]*sq(SK_LOS[4]) - P[13][4]*SH_LOS[0]*SK_LOS[4]*(2*q0*q3 + 2*q1*q2) - P[13][2]*SH_LOS[0]*SH_LOS[2]*SK_LOS[4] + P[13][0]*SH_LOS[1]*SH_LOS[6]*SK_LOS[4] + P[13][5]*SH_LOS[0]*SH_LOS[5]*SK_LOS[4] - P[13][3]*SH_LOS[0]*SK_LOS[2]*SK_LOS[4]);
			Kfusion[14] = -SK_LOS[1]*(P[14][1]*(SH_LOS[0]*SH_LOS[4]*SK_LOS[4] + SH_LOS[1]*SH_LOS[5]*SK_LOS[4]) + P[14][8]*SH_LOS[0]*SH_LOS[1]*sq(SK_LOS[4]) - P[14][4]*SH_LOS[0]*SK_LOS[4]*(2*q0*q3 + 2*q1*q2) - P[14][2]*SH_LOS[0]*SH_LOS[2]*SK_LOS[4] + P[14][0]*SH_LOS[1]*SH_LOS[6]*SK_LOS[4] + P[14][5]*SH_LOS[0]*SH_LOS[5]*SK_LOS[4] - P[14][3]*SH_LOS[0]*SK_LOS[2]*SK_LOS[4]);
			Kfusion[15] = -SK_LOS[1]*(P[15][1]*(SH_LOS[0]*SH_LOS[4]*SK_LOS[4] + SH_LOS[1]*SH_LOS[5]*SK_LOS[4]) + P[15][8]*SH_LOS[0]*SH_LOS[1]*sq(SK_LOS[4]) - P[15][4]*SH_LOS[0]*SK_LOS[4]*(2*q0*q3 + 2*q1*q2) - P[15][2]*SH_LOS[0]*SH_LOS[2]*SK_LOS[4] + P[15][0]*SH_LOS[1]*SH_LOS[6]*SK_LOS[4] + P[15][5]*SH_LOS[0]*SH_LOS[5]*SK_LOS[4] - P[15][3]*SH_LOS[0]*SK_LOS[2]*SK_LOS[4]);
			Kfusion[16] = -SK_LOS[1]*(P[16][1]*(SH_LOS[0]*SH_LOS[4]*SK_LOS[4] + SH_LOS[1]*SH_LOS[5]*SK_LOS[4]) + P[16][8]*SH_LOS[0]*SH_LOS[1]*sq(SK_LOS[4]) - P[16][4]*SH_LOS[0]*SK_LOS[4]*(2*q0*q3 + 2*q1*q2) - P[16][2]*SH_LOS[0]*SH_LOS[2]*SK_LOS[4] + P[16][0]*SH_LOS[1]*SH_LOS[6]*SK_LOS[4] + P[16][5]*SH_LOS[0]*SH_LOS[5]*SK_LOS[4] - P[16][3]*SH_LOS[0]*SK_LOS[2]*SK_LOS[4]);
			Kfusion[17] = -SK_LOS[1]*(P[17][1]*(SH_LOS[0]*SH_LOS[4]*SK_LOS[4] + SH_LOS[1]*SH_LOS[5]*SK_LOS[4]) + P[17][8]*SH_LOS[0]*SH_LOS[1]*sq(SK_LOS[4]) - P[17][4]*SH_LOS[0]*SK_LOS[4]*(2*q0*q3 + 2*q1*q2) - P[17][2]*SH_LOS[0]*SH_LOS[2]*SK_LOS[4] + P[17][0]*SH_LOS[1]*SH_LOS[6]*SK_LOS[4] + P[17][5]*SH_LOS[0]*SH_LOS[5]*SK_LOS[4] - P[17][3]*SH_LOS[0]*SK_LOS[2]*SK_LOS[4]);
			Kfusion[18] = -SK_LOS[1]*(P[18][1]*(SH_LOS[0]*SH_LOS[4]*SK_LOS[4] + SH_LOS[1]*SH_LOS[5]*SK_LOS[4]) + P[18][8]*SH_LOS[0]*SH_LOS[1]*sq(SK_LOS[4]) - P[18][4]*SH_LOS[0]*SK_LOS[4]*(2*q0*q3 + 2*q1*q2) - P[18][2]*SH_LOS[0]*SH_LOS[2]*SK_LOS[4] + P[18][0]*SH_LOS[1]*SH_LOS[6]*SK_LOS[4] + P[18][5]*SH_LOS[0]*SH_LOS[5]*SK_LOS[4] - P[18][3]*SH_LOS[0]*SK_LOS[2]*SK_LOS[4]);
			Kfusion[19] = -SK_LOS[1]*(P[19][1]*(SH_LOS[0]*SH_LOS[4]*SK_LOS[4] + SH_LOS[1]*SH_LOS[5]*SK_LOS[4]) + P[19][8]*SH_LOS[0]*SH_LOS[1]*sq(SK_LOS[4]) - P[19][4]*SH_LOS[0]*SK_LOS[4]*(2*q0*q3 + 2*q1*q2) - P[19][2]*SH_LOS[0]*SH_LOS[2]*SK_LOS[4] + P[19][0]*SH_LOS[1]*SH_LOS[6]*SK_LOS[4] + P[19][5]*SH_LOS[0]*SH_LOS[5]*SK_LOS[4] - P[19][3]*SH_LOS[0]*SK_LOS[2]*SK_LOS[4]);
			Kfusion[20] = -SK_LOS[1]*(P[20][1]*(SH_LOS[0]*SH_LOS[4]*SK_LOS[4] + SH_LOS[1]*SH_LOS[5]*SK_LOS[4]) + P[20][8]*SH_LOS[0]*SH_LOS[1]*sq(SK_LOS[4]) - P[20][4]*SH_LOS[0]*SK_LOS[4]*(2*q0*q3 + 2*q1*q2) - P[20][2]*SH_LOS[0]*SH_LOS[2]*SK_LOS[4] + P[20][0]*SH_LOS[1]*SH_LOS[6]*SK_LOS[4] + P[20][5]*SH_LOS[0]*SH_LOS[5]*SK_LOS[4] - P[20][3]*SH_LOS[0]*SK_LOS[2]*SK_LOS[4]);
			Kfusion[21] = -SK_LOS[1]*(P[21][1]*(SH_LOS[0]*SH_LOS[4]*SK_LOS[4] + SH_LOS[1]*SH_LOS[5]*SK_LOS[4]) + P[21][8]*SH_LOS[0]*SH_LOS[1]*sq(SK_LOS[4]) - P[21][4]*SH_LOS[0]*SK_LOS[4]*(2*q0*q3 + 2*q1*q2) - P[21][2]*SH_LOS[0]*SH_LOS[2]*SK_LOS[4] + P[21][0]*SH_LOS[1]*SH_LOS[6]*SK_LOS[4] + P[21][5]*SH_LOS[0]*SH_LOS[5]*SK_LOS[4] - P[21][3]*SH_LOS[0]*SK_LOS[2]*SK_LOS[4]);
			Kfusion[22] = -SK_LOS[1]*(P[22][1]*(SH_LOS[0]*SH_LOS[4]*SK_LOS[4] + SH_LOS[1]*SH_LOS[5]*SK_LOS[4]) + P[22][8]*SH_LOS[0]*SH_LOS[1]*sq(SK_LOS[4]) - P[22][4]*SH_LOS[0]*SK_LOS[4]*(2*q0*q3 + 2*q1*q2) - P[22][2]*SH_LOS[0]*SH_LOS[2]*SK_LOS[4] + P[22][0]*SH_LOS[1]*SH_LOS[6]*SK_LOS[4] + P[22][5]*SH_LOS[0]*SH_LOS[5]*SK_LOS[4] - P[22][3]*SH_LOS[0]*SK_LOS[2]*SK_LOS[4]);
			Kfusion[23] = -SK_LOS[1]*(P[23][1]*(SH_LOS[0]*SH_LOS[4]*SK_LOS[4] + SH_LOS[1]*SH_LOS[5]*SK_LOS[4]) + P[23][8]*SH_LOS[0]*SH_LOS[1]*sq(SK_LOS[4]) - P[23][4]*SH_LOS[0]*SK_LOS[4]*(2*q0*q3 + 2*q1*q2) - P[23][2]*SH_LOS[0]*SH_LOS[2]*SK_LOS[4] + P[23][0]*SH_LOS[1]*SH_LOS[6]*SK_LOS[4] + P[23][5]*SH_LOS[0]*SH_LOS[5]*SK_LOS[4] - P[23][3]*SH_LOS[0]*SK_LOS[2]*SK_LOS[4]);
		}
		// by definition our error state is zero at the time of fusion
		_state.ang_error.setZero();
	    // only update magnetic field states if we are fusing 3-axis observations
        if (!_control_status.flags.mag_3D) {
            for (int row = 16; row <= 21; row++) {
                Kfusion[row] = 0.0f;
            }
        }
        // only update wind states if we are doing wind estimation
        if (!_control_status.flags.wind) {
            for (int row = 22; row <= 23; row++) {
                Kfusion[row] = 0.0f;
            }
        }	
        fuse(Kfusion,_flow_innov(index));

		Quaternion q_correction;
		q_correction.from_axis_angle(_state.ang_error);
		_state.quat_nominal = q_correction * _state.quat_nominal;
		_state.quat_nominal.normalize();
		_state.ang_error.setZero();

		// apply covariance correction via P_new = (I -K*H)*P
		// first calculate expression for KHP
		// then calculate P - KHP
		float KH[_k_num_states][_k_num_states] = {0};

		for (unsigned row = 0; row < _k_num_states; row++) {
			for (unsigned column = 0; column < 5; column++) {
				KH[row][column] = Kfusion[row] * H_LOS[column];
			}
			KH[row][8] = Kfusion[row] * H_LOS[8];
		}

		float KHP[_k_num_states][_k_num_states] = {0};

		for (unsigned row = 0; row < _k_num_states; row++) {
			for (unsigned column = 0; column < _k_num_states; column++) {
				float tmp = KH[row][0] * P[0][column];
				tmp += KH[row][1] * P[1][column];
				tmp += KH[row][2] * P[2][column];
				tmp += KH[row][3] * P[3][column];
				tmp += KH[row][4] * P[4][column];
				tmp += KH[row][5] * P[5][column];
				tmp += KH[row][8] * P[8][column];
				KHP[row][column] = tmp;
			}
		}

		for (unsigned row = 0; row < _k_num_states; row++) {
			for (unsigned column = 0; column < _k_num_states; column++) {
				P[row][column] -= KHP[row][column];
			}
		}
		_time_last_of_fuse = _time_last_imu;
		_gps_check_fail_status.value = 0;
		makeSymmetrical();
		limitCov();
	}
}