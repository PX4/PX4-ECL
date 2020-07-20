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
 * @file heading_fusion.cpp
 * Magnetometer fusion methods.
 * Equations generated using EKF/python/ekf_derivation/main.py
 *
 * @author Roman Bast <bapstroman@gmail.com>
 * @author Paul Riseborough <p_riseborough@live.com.au>
 *
 */

#include "ekf.h"
#include <ecl.h>
#include <mathlib/mathlib.h>

void Ekf::fuseMag()
{
	// assign intermediate variables
	const float &q0 = _state.quat_nominal(0);
	const float &q1 = _state.quat_nominal(1);
	const float &q2 = _state.quat_nominal(2);
	const float &q3 = _state.quat_nominal(3);

	const float &magN = _state.mag_I(0);
	const float &magE = _state.mag_I(1);
	const float &magD = _state.mag_I(2);

	// XYZ Measurement uncertainty. Need to consider timing errors for fast rotations
	const float R_MAG = sq(fmaxf(_params.mag_noise, 0.0f));

	// intermediate variables from algebraic optimisation
	const char* numerical_error_covariance_reset_string = "numerical error - covariance reset";
	const float HK0 = magE*q3;
	const float HK1 = magN*q0;
	const float HK2 = magD*q2;
	const float HK3 = HK0 + HK1 - HK2;
	const float HK4 = 2*HK3;
	const float HK5 = magD*q3 + magE*q2 + magN*q1;
	const float HK6 = 2*HK5;
	const float HK7 = magE*q1;
	const float HK8 = magD*q0;
	const float HK9 = magN*q2;
	const float HK10 = magD*q1;
	const float HK11 = magE*q0;
	const float HK12 = magN*q3;
	const float HK13 = HK10 + HK11 - HK12;
	const float HK14 = 2*HK13;
	const float HK15 = powf(q1, 2);
	const float HK16 = powf(q2, 2);
	const float HK17 = -HK16;
	const float HK18 = powf(q0, 2);
	const float HK19 = powf(q3, 2);
	const float HK20 = HK18 - HK19;
	const float HK21 = HK15 + HK17 + HK20;
	const float HK22 = q0*q3;
	const float HK23 = q1*q2;
	const float HK24 = HK22 + HK23;
	const float HK25 = q1*q3;
	const float HK26 = q0*q2;
	const float HK27 = 2*HK24;
	const float HK28 = -2*HK25 + 2*HK26;
	const float HK29 = 2*HK5;
	const float HK30 = 2*HK3;
	const float HK31 = -HK7 + HK8 + HK9;
	const float HK32 = 2*HK31;
	const float HK33 = HK32*P(0,2);
	const float HK34 = 2*HK13;
	const float HK35 = HK34*P(0,3);
	const float HK36 = HK21*P(0,16) + HK27*P(0,17) - HK28*P(0,18) + HK29*P(0,1) + HK30*P(0,0) - HK33 + HK35 + P(0,19);
	const float HK37 = HK21*P(16,16) + HK27*P(16,17) - HK28*P(16,18) + HK29*P(1,16) + HK30*P(0,16) - HK32*P(2,16) + HK34*P(3,16) + P(16,19);
	const float HK38 = HK21*P(16,18) + HK27*P(17,18) - HK28*P(18,18) + HK29*P(1,18) + HK30*P(0,18) - HK32*P(2,18) + HK34*P(3,18) + P(18,19);
	const float HK39 = HK29*P(1,2);
	const float HK40 = HK30*P(0,2);
	const float HK41 = HK21*P(2,16) + HK27*P(2,17) - HK28*P(2,18) - HK32*P(2,2) + HK34*P(2,3) + HK39 + HK40 + P(2,19);
	const float HK42 = HK21*P(16,17) + HK27*P(17,17) - HK28*P(17,18) + HK29*P(1,17) + HK30*P(0,17) - HK32*P(2,17) + HK34*P(3,17) + P(17,19);
	const float HK43 = HK29*P(1,3);
	const float HK44 = HK30*P(0,3);
	const float HK45 = HK21*P(3,16) + HK27*P(3,17) - HK28*P(3,18) - HK32*P(2,3) + HK34*P(3,3) + HK43 + HK44 + P(3,19);
	const float HK46 = HK32*P(1,2);
	const float HK47 = HK34*P(1,3);
	const float HK48 = HK21*P(1,16) + HK27*P(1,17) - HK28*P(1,18) + HK29*P(1,1) + HK30*P(0,1) - HK46 + HK47 + P(1,19);
	const float HK49 = HK21*P(16,19) + HK27*P(17,19) - HK28*P(18,19) + HK29*P(1,19) + HK30*P(0,19) - HK32*P(2,19) + HK34*P(3,19) + P(19,19);

	_mag_innov_var(0) = (HK21*HK37 + HK27*HK42 - HK28*HK38 + HK29*HK48 + HK30*HK36 - HK32*HK41 + HK34*HK45 + HK49 + R_MAG);
	float HK50;
	if (_mag_innov_var(0) >= R_MAG) {
		// the innovation variance contribution from the state covariances is non-negative - no fault
		_fault_status.flags.bad_mag_x = false;
		HK50 = 1.0F/_mag_innov_var(0);

	} else {
		// the innovation variance contribution from the state covariances is negative which means the covariance matrix is badly conditioned
		_fault_status.flags.bad_mag_x = true;

		// we need to re-initialise covariances and abort this fusion step
		resetMagRelatedCovariances();
		ECL_ERR("magX %s", numerical_error_covariance_reset_string);
		return;

	}

	const float HK51 = 2*HK31;
	const float HK52 = -HK15;
	const float HK53 = HK16 + HK20 + HK52;
	const float HK54 = q0*q1;
	const float HK55 = q2*q3;
	const float HK56 = HK54 + HK55;
	const float HK57 = 2*HK56;
	const float HK58 = 2*HK22 - 2*HK23;
	const float HK59 = HK32*P(0,1);
	const float HK60 = HK29*P(0,2) + HK34*P(0,0) - HK44 + HK53*P(0,17) + HK57*P(0,18) - HK58*P(0,16) + HK59 + P(0,20);
	const float HK61 = HK29*P(2,17) - HK30*P(3,17) + HK32*P(1,17) + HK34*P(0,17) + HK53*P(17,17) + HK57*P(17,18) - HK58*P(16,17) + P(17,20);
	const float HK62 = HK29*P(2,16) - HK30*P(3,16) + HK32*P(1,16) + HK34*P(0,16) + HK53*P(16,17) + HK57*P(16,18) - HK58*P(16,16) + P(16,20);
	const float HK63 = HK29*P(2,3);
	const float HK64 = -HK30*P(3,3) + HK32*P(1,3) + HK35 + HK53*P(3,17) + HK57*P(3,18) - HK58*P(3,16) + HK63 + P(3,20);
	const float HK65 = HK29*P(2,18) - HK30*P(3,18) + HK32*P(1,18) + HK34*P(0,18) + HK53*P(17,18) + HK57*P(18,18) - HK58*P(16,18) + P(18,20);
	const float HK66 = HK34*P(0,1);
	const float HK67 = -HK30*P(1,3) + HK32*P(1,1) + HK39 + HK53*P(1,17) + HK57*P(1,18) - HK58*P(1,16) + HK66 + P(1,20);
	const float HK68 = HK30*P(2,3);
	const float HK69 = HK29*P(2,2) + HK34*P(0,2) + HK46 + HK53*P(2,17) + HK57*P(2,18) - HK58*P(2,16) - HK68 + P(2,20);
	const float HK70 = HK29*P(2,20) - HK30*P(3,20) + HK32*P(1,20) + HK34*P(0,20) + HK53*P(17,20) + HK57*P(18,20) - HK58*P(16,20) + P(20,20);

	_mag_innov_var(1) = (HK29*HK69 - HK30*HK64 + HK32*HK67 + HK34*HK60 + HK53*HK61 + HK57*HK65 - HK58*HK62 + HK70 + R_MAG);
	float HK71;
	// check for a badly conditioned covariance matrix
	if (_mag_innov_var(1) >= R_MAG) {
		// the innovation variance contribution from the state covariances is non-negative - no fault
		_fault_status.flags.bad_mag_y = false;
		HK71 = 1.0F/_mag_innov_var(1);

	} else {
		// the innovation variance contribution from the state covariances is negtive which means the covariance matrix is badly conditioned
		_fault_status.flags.bad_mag_y = true;

		// we need to re-initialise covariances and abort this fusion step
		resetMagRelatedCovariances();
		ECL_ERR("magY %s", numerical_error_covariance_reset_string);
		return;

	}

	const float HK72 = HK25 + HK26;
	const float HK73 = HK17 + HK18 + HK19 + HK52;
	const float HK74 = 2*HK72;
	const float HK75 = 2*HK54 - 2*HK55;
	const float HK76 = HK29*P(0,3) + HK32*P(0,0) + HK40 - HK66 + HK73*P(0,18) + HK74*P(0,16) - HK75*P(0,17) + P(0,21);
	const float HK77 = HK29*P(3,18) + HK30*P(2,18) + HK32*P(0,18) - HK34*P(1,18) + HK73*P(18,18) + HK74*P(16,18) - HK75*P(17,18) + P(18,21);
	const float HK78 = HK29*P(3,17) + HK30*P(2,17) + HK32*P(0,17) - HK34*P(1,17) + HK73*P(17,18) + HK74*P(16,17) - HK75*P(17,17) + P(17,21);
	const float HK79 = HK30*P(1,2) - HK34*P(1,1) + HK43 + HK59 + HK73*P(1,18) + HK74*P(1,16) - HK75*P(1,17) + P(1,21);
	const float HK80 = HK29*P(3,16) + HK30*P(2,16) + HK32*P(0,16) - HK34*P(1,16) + HK73*P(16,18) + HK74*P(16,16) - HK75*P(16,17) + P(16,21);
	const float HK81 = HK29*P(3,3) + HK32*P(0,3) - HK47 + HK68 + HK73*P(3,18) + HK74*P(3,16) - HK75*P(3,17) + P(3,21);
	const float HK82 = HK30*P(2,2) + HK33 - HK34*P(1,2) + HK63 + HK73*P(2,18) + HK74*P(2,16) - HK75*P(2,17) + P(2,21);
	const float HK83 = HK29*P(3,21) + HK30*P(2,21) + HK32*P(0,21) - HK34*P(1,21) + HK73*P(18,21) + HK74*P(16,21) - HK75*P(17,21) + P(21,21);

	_mag_innov_var(2) = (HK29*HK81 + HK30*HK82 + HK32*HK76 - HK34*HK79 + HK73*HK77 + HK74*HK80 - HK75*HK78 + HK83 + R_MAG);
	float HK84;
	// check for a badly conditioned covariance matrix
	if (_mag_innov_var(2) >= R_MAG) {
		// the innovation variance contribution from the state covariances is non-negative - no fault
		_fault_status.flags.bad_mag_z = false;
		HK84 = 1.0F/_mag_innov_var(2);

	} else {
		// the innovation variance contribution from the state covariances is negative which means the covariance matrix is badly conditioned
		_fault_status.flags.bad_mag_z = true;

		// we need to re-initialise covariances and abort this fusion step
		resetMagRelatedCovariances();
		ECL_ERR("magZ %s", numerical_error_covariance_reset_string);
		return;

	}

	// rotate magnetometer earth field state into body frame
	const Dcmf R_to_body = quatToInverseRotMat(_state.quat_nominal);

	const Vector3f mag_I_rot = R_to_body * _state.mag_I;

	// compute magnetometer innovations
	_mag_innov = mag_I_rot + _state.mag_B - _mag_sample_delayed.mag;

	// do not use the synthesized measurement for the magnetomter Z component for 3D fusion
	if (_control_status.flags.synthetic_mag_z) {
		_mag_innov(2) = 0.0f;
	}

	// Perform an innovation consistency check and report the result
	bool all_innovation_checks_passed = true;

	for (uint8_t index = 0; index <= 2; index++) {
		_mag_test_ratio(index) = sq(_mag_innov(index)) / (sq(math::max(_params.mag_innov_gate, 1.0f)) * _mag_innov_var(index));

		if (_mag_test_ratio(index) > 1.0f) {
			all_innovation_checks_passed = false;
			_innov_check_fail_status.value |= (1 << (index + 3));

		} else {
			_innov_check_fail_status.value &= ~(1 << (index + 3));
		}
	}

	// we are no longer using heading fusion so set the reported test level to zero
	_yaw_test_ratio = 0.0f;

	// if any axis fails, abort the mag fusion
	if (!all_innovation_checks_passed) {
		return;
	}

	// For the first few seconds after in-flight alignment we allow the magnetic field state estimates to stabilise
	// before they are used to constrain heading drift
	const bool update_all_states = ((_imu_sample_delayed.time_us - _flt_mag_align_start_time) > (uint64_t)5e6);

	// Observation jacobian and Kalman gain vectors
	float Hfusion[24];
	Vector24f Kfusion;

	// update the states and covariance using sequential fusion of the magnetometer components
	for (uint8_t index = 0; index <= 2; index++) {

		// Calculate Kalman gains and observation jacobians
		if (index == 0) {
			// Calculate X axis observation jacobians
			memset(Hfusion, 0, sizeof(Hfusion));
			Hfusion[0] = HK4;
			Hfusion[1] = HK6;
			Hfusion[2] = 2*HK7 - 2*HK8 - 2*HK9;
			Hfusion[3] = HK14;
			Hfusion[16] = HK21;
			Hfusion[17] = 2*HK24;
			Hfusion[18] = 2*HK25 - 2*HK26;
			Hfusion[19] = 1;

			// Calculate X axis Kalman gains
			if (update_all_states) {
				Kfusion(0) = HK36*HK50;
				Kfusion(1) = HK48*HK50;
				Kfusion(2) = HK41*HK50;
				Kfusion(3) = HK45*HK50;

				for (unsigned row = 4; row <= 15; row++) {
					Kfusion(row) = HK50*(HK21*P(row,16) + HK27*P(row,17) - HK28*P(row,18) + HK29*P(1,row) + HK30*P(0,row) - HK32*P(2,row) + HK34*P(3,row) + P(row,19));
				}

				for (unsigned row = 22; row <= 23; row++) {
					Kfusion(row) = HK50*(HK21*P(row,16) + HK27*P(row,17) - HK28*P(row,18) + HK29*P(1,row) + HK30*P(0,row) - HK32*P(2,row) + HK34*P(3,row) + P(row,19));
				}
			} else {
				for (uint8_t i = 0; i < 16; i++) {
					Kfusion(i) = 0.0f;
				}

				Kfusion(22) = 0.0f;
				Kfusion(23) = 0.0f;
			}

			Kfusion(16) = HK37*HK50;
			Kfusion(17) = HK42*HK50;
			Kfusion(18) = HK38*HK50;
			Kfusion(19) = HK49*HK50;

			for (unsigned row = 20; row <= 21; row++) {
				Kfusion(row) = HK50*(HK21*P(16,row) + HK27*P(17,row) - HK28*P(18,row) + HK29*P(1,row) + HK30*P(0,row) - HK32*P(2,row) + HK34*P(3,row) + P(19,row));
			}
		} else if (index == 1) {
			// Calculate Y axis observation jacobians
			memset(Hfusion, 0, sizeof(Hfusion));
			Hfusion[0] = HK14;
			Hfusion[1] = HK51;
			Hfusion[2] = HK6;
			Hfusion[3] = -2*HK0 - 2*HK1 + 2*HK2;
			Hfusion[16] = -2*HK22 + 2*HK23;
			Hfusion[17] = HK53;
			Hfusion[18] = 2*HK56;
			Hfusion[20] = 1;

			// Calculate Y axis Kalman gains
			if (update_all_states) {
				Kfusion(0) = HK60*HK71;
				Kfusion(1) = HK67*HK71;
				Kfusion(2) = HK69*HK71;
				Kfusion(3) = HK64*HK71;

				for (unsigned row = 4; row <= 15; row++) {
					Kfusion(row) = HK71*(HK29*P(2,row) - HK30*P(3,row) + HK32*P(1,row) + HK34*P(0,row) + HK53*P(row,17) + HK57*P(row,18) - HK58*P(row,16) + P(row,20));
				}

				for (unsigned row = 22; row <= 23; row++) {
					Kfusion(row) = HK71*(HK29*P(2,row) - HK30*P(3,row) + HK32*P(1,row) + HK34*P(0,row) + HK53*P(row,17) + HK57*P(row,18) - HK58*P(row,16) + P(row,20));
				}
			} else {
				for (uint8_t i = 0; i < 16; i++) {
					Kfusion(i) = 0.0f;
				}

				Kfusion(22) = 0.0f;
				Kfusion(23) = 0.0f;
			}

			Kfusion(16) = HK62*HK71;
			Kfusion(17) = HK61*HK71;
			Kfusion(18) = HK65*HK71;
			Kfusion(19) = HK71*(HK29*P(2,19) - HK30*P(3,19) + HK32*P(1,19) + HK34*P(0,19) + HK53*P(17,19) + HK57*P(18,19) - HK58*P(16,19) + P(19,20));
			Kfusion(20) = HK70*HK71;
			Kfusion(21) = HK71*(HK29*P(2,21) - HK30*P(3,21) + HK32*P(1,21) + HK34*P(0,21) + HK53*P(17,21) + HK57*P(18,21) - HK58*P(16,21) + P(20,21));

		} else if (index == 2) {

			// we do not fuse synthesized magnetomter measurements when doing 3D fusion
			if (_control_status.flags.synthetic_mag_z) {
				continue;
			}

			// calculate Z axis observation jacobians
			memset(Hfusion, 0, sizeof(Hfusion));
			Hfusion[0] = HK51;
			Hfusion[1] = -2*HK10 - 2*HK11 + 2*HK12;
			Hfusion[2] = HK4;
			Hfusion[3] = HK6;
			Hfusion[16] = 2*HK72;
			Hfusion[17] = -2*HK54 + 2*HK55;
			Hfusion[18] = HK73;
			Hfusion[21] = 1;

			// Calculate Z axis Kalman gains
			if (update_all_states) {
				Kfusion(0) = HK76*HK84;
				Kfusion(1) = HK79*HK84;
				Kfusion(2) = HK82*HK84;
				Kfusion(3) = HK81*HK84;

				for (unsigned row = 4; row <= 15; row++) {
					Kfusion(row) = HK84*(HK29*P(3,row) + HK30*P(2,row) + HK32*P(0,row) - HK34*P(1,row) + HK73*P(row,18) + HK74*P(row,16) - HK75*P(row,17) + P(row,21));
				}

				for (unsigned row = 22; row <= 23; row++) {
					Kfusion(row) = HK84*(HK29*P(3,row) + HK30*P(2,row) + HK32*P(0,row) - HK34*P(1,row) + HK73*P(row,18) + HK74*P(row,16) - HK75*P(row,17) + P(row,21));
				}
			} else {
				for (uint8_t i = 0; i < 16; i++) {
					Kfusion(i) = 0.0f;
				}

				Kfusion(22) = 0.0f;
				Kfusion(23) = 0.0f;
			}

			Kfusion(16) = HK80*HK84;
			Kfusion(17) = HK78*HK84;
			Kfusion(18) = HK77*HK84;

			for (unsigned row = 19; row <= 20; row++) {
			Kfusion(row) = HK84*(HK29*P(3,row) + HK30*P(2,row) + HK32*P(0,row) - HK34*P(1,row) + HK73*P(18,row) + HK74*P(16,row) - HK75*P(17,row) + P(row,21));
			}

			Kfusion(21) = HK83*HK84;
		}

		// apply covariance correction via P_new = (I -K*H)*P
		// first calculate expression for KHP
		// then calculate P - KHP
		SquareMatrix24f KHP;
		float KH[10];

		for (unsigned row = 0; row < _k_num_states; row++) {

			KH[0] = Kfusion(row) * Hfusion[0];
			KH[1] = Kfusion(row) * Hfusion[1];
			KH[2] = Kfusion(row) * Hfusion[2];
			KH[3] = Kfusion(row) * Hfusion[3];
			KH[4] = Kfusion(row) * Hfusion[16];
			KH[5] = Kfusion(row) * Hfusion[17];
			KH[6] = Kfusion(row) * Hfusion[18];
			KH[7] = Kfusion(row) * Hfusion[19];
			KH[8] = Kfusion(row) * Hfusion[20];
			KH[9] = Kfusion(row) * Hfusion[21];

			for (unsigned column = 0; column < _k_num_states; column++) {
				float tmp = KH[0] * P(0,column);
				tmp += KH[1] * P(1,column);
				tmp += KH[2] * P(2,column);
				tmp += KH[3] * P(3,column);
				tmp += KH[4] * P(16,column);
				tmp += KH[5] * P(17,column);
				tmp += KH[6] * P(18,column);
				tmp += KH[7] * P(19,column);
				tmp += KH[8] * P(20,column);
				tmp += KH[9] * P(21,column);
				KHP(row,column) = tmp;
			}
		}

		const bool healthy = checkAndFixCovarianceUpdate(KHP);

		if (index == 0) {
			_fault_status.flags.bad_mag_x = !healthy;

		} else if (index == 1) {
			_fault_status.flags.bad_mag_y = !healthy;

		} else if (index == 2) {
			_fault_status.flags.bad_mag_z = !healthy;
		}

		if (healthy) {
			// apply the covariance corrections
			P -= KHP;

			fixCovarianceErrors(true);

			// apply the state corrections
			fuse(Kfusion, _mag_innov(index));

			// constrain the declination of the earth field states
			limitDeclination();
		}
	}
}

void Ekf::fuseYaw321(float yaw, float yaw_variance, bool zero_innovation)
{
	// assign intermediate state variables
	const float &q0 = _state.quat_nominal(0);
	const float &q1 = _state.quat_nominal(1);
	const float &q2 = _state.quat_nominal(2);
	const float &q3 = _state.quat_nominal(3);

	const float R_YAW = fmaxf(yaw_variance, 1.0e-4f);
	const float measurement = wrap_pi(yaw);
	float H_YAW[4];

    // calculate 321 yaw observation matrix
	// choose A or B computational paths to avoid singularity in derivation at +-90 degrees yaw
	bool canUseA = false;
	const float SA0 = 2*q3;
	const float SA1 = 2*q2;
	const float SA2 = SA0*q0 + SA1*q1;
	const float SA3 = sq(q0) + sq(q1) - sq(q2) - sq(q3);
	float SA4, SA5_inv;
	if (sq(SA3) > 1E-6f) {
		SA4 = 1.0F/sq(SA3);
		SA5_inv = sq(SA2)*SA4 + 1;
		canUseA = fabsf(SA5_inv) > 1E-6f;
	}

	bool canUseB = false;
	const float SB0 = 2*q0;
	const float SB1 = 2*q1;
	const float SB2 = SB0*q3 + SB1*q2;
	const float SB4 = sq(q0) + sq(q1) - sq(q2) - sq(q3);
	float SB3, SB5_inv;
	if (sq(SB2) > 1E-6f) {
		SB3 = 1.0F/sq(SB2);
		SB5_inv = SB3*sq(SB4) + 1;
		canUseB = fabsf(SB5_inv) > 1E-6f;
	}

	if (canUseA && (!canUseB || fabsf(SA5_inv) >= fabsf(SB5_inv))) {
		const float SA5 = 1.0F/SA5_inv;
		const float SA6 = 1.0F/SA3;
		const float SA7 = SA2*SA4;
		const float SA8 = 2*SA7;
		const float SA9 = 2*SA6;

		H_YAW[0] = SA5*(SA0*SA6 - SA8*q0);
		H_YAW[1] = SA5*(SA1*SA6 - SA8*q1);
		H_YAW[2] = SA5*(SA1*SA7 + SA9*q1);
		H_YAW[3] = SA5*(SA0*SA7 + SA9*q0);
	} else if (canUseB && (!canUseA || fabsf(SB5_inv) > fabsf(SA5_inv))) {
		const float SB5 = 1.0F/SB5_inv;
		const float SB6 = 1.0F/SB2;
		const float SB7 = SB3*SB4;
		const float SB8 = 2*SB7;
		const float SB9 = 2*SB6;

		H_YAW[0] = -SB5*(SB0*SB6 - SB8*q3);
		H_YAW[1] = -SB5*(SB1*SB6 - SB8*q2);
		H_YAW[2] = -SB5*(-SB1*SB7 - SB9*q2);
		H_YAW[3] = -SB5*(-SB0*SB7 - SB9*q3);
	} else {
		return;
	}

	// calculate the yaw innovation and wrap to the interval between +-pi
	float innovation;
	if (zero_innovation) {
		innovation = 0.0f;
	} else {
		innovation = wrap_pi(atan2f(_R_to_earth(1, 0), _R_to_earth(0, 0)) - measurement);
	}

	// define the innovation gate size
	float innov_gate = math::max(_params.heading_innov_gate, 1.0f);

	// Update the quaternion states and covariance matrix
	updateQuaternion(innovation, R_YAW, innov_gate, H_YAW);
}

void Ekf::fuseYaw312(float yaw, float yaw_variance, bool zero_innovation)
{
	// assign intermediate state variables
	const float q0 = _state.quat_nominal(0);
	const float q1 = _state.quat_nominal(1);
	const float q2 = _state.quat_nominal(2);
	const float q3 = _state.quat_nominal(3);

	const float R_YAW = fmaxf(yaw_variance, 1.0e-4f);
	const float measurement = wrap_pi(yaw);
	float H_YAW[4];

    // calculate 312 yaw observation matrix
	// choose A or B computational paths to avoid singularity in derivation at +-90 degrees yaw
	bool canUseA = false;
	const float SA0 = 2*q3;
	const float SA1 = 2*q2;
	const float SA2 = SA0*q0 - SA1*q1;
	const float SA3 = sq(q0) - sq(q1) + sq(q2) - sq(q3);
	float SA4, SA5_inv;
	if (sq(SA3) > 1E-6f) {
		SA4 = 1.0F/sq(SA3);
		SA5_inv = sq(SA2)*SA4 + 1;
		canUseA = fabsf(SA5_inv) > 1E-6f;
	}

	bool canUseB = false;
	const float SB0 = 2*q0;
	const float SB1 = 2*q1;
	const float SB2 = -SB0*q3 + SB1*q2;
	const float SB4 = -sq(q0) + sq(q1) - sq(q2) + sq(q3);
	float SB3, SB5_inv;
	if (sq(SB2) > 1E-6f) {
		SB3 = 1.0F/sq(SB2);
		SB5_inv = SB3*sq(SB4) + 1;
		canUseB = fabsf(SB5_inv) > 1E-6f;
	}

	if (canUseA && (!canUseB || fabsf(SA5_inv) >= fabsf(SB5_inv))) {
		const float SA5 = 1.0F/SA5_inv;
		const float SA6 = 1.0F/SA3;
		const float SA7 = SA2*SA4;
		const float SA8 = 2*SA7;
		const float SA9 = 2*SA6;

		H_YAW[0] = SA5*(SA0*SA6 - SA8*q0);
		H_YAW[1] = SA5*(-SA1*SA6 + SA8*q1);
		H_YAW[2] = SA5*(-SA1*SA7 - SA9*q1);
		H_YAW[3] = SA5*(SA0*SA7 + SA9*q0);
	} else if (canUseB && (!canUseA || fabsf(SB5_inv) > fabsf(SA5_inv))) {
		const float SB5 = 1.0F/SB5_inv;
		const float SB6 = 1.0F/SB2;
		const float SB7 = SB3*SB4;
		const float SB8 = 2*SB7;
		const float SB9 = 2*SB6;

		H_YAW[0] = -SB5*(-SB0*SB6 + SB8*q3);
		H_YAW[1] = -SB5*(SB1*SB6 - SB8*q2);
		H_YAW[2] = -SB5*(-SB1*SB7 - SB9*q2);
		H_YAW[3] = -SB5*(SB0*SB7 + SB9*q3);
	} else {
		return;
	}

	float innovation;
	if (zero_innovation) {
		innovation = 0.0f;
	} else {
		// calculate the the innovation and wrap to the interval between +-pi
		innovation = wrap_pi(atan2f(-_R_to_earth(0, 1), _R_to_earth(1, 1)) - measurement);
	}

	// define the innovation gate size
	float innov_gate = math::max(_params.heading_innov_gate, 1.0f);

	// Update the quaternion states and covariance matrix
	updateQuaternion(innovation, R_YAW, innov_gate, H_YAW);
}

// update quaternion states and covariances using the yaw innovation, yaw observation variance and yaw Jacobian
void Ekf::updateQuaternion(const float innovation, const float variance, const float gate_sigma, const float (&yaw_jacobian)[4])
{
	// Calculate innovation variance and Kalman gains, taking advantage of the fact that only the first 4 elements in H are non zero
	// calculate the innovation variance
	float PH[4];
	_heading_innov_var = variance;
	for (unsigned row = 0; row <= 3; row++) {
		PH[row] = 0.0f;

		for (uint8_t col = 0; col <= 3; col++) {
			PH[row] += P(row,col) * yaw_jacobian[col];
		}

		_heading_innov_var += yaw_jacobian[row] * PH[row];
	}

	float heading_innov_var_inv;

	// check if the innovation variance calculation is badly conditioned
	if (_heading_innov_var >= variance) {
		// the innovation variance contribution from the state covariances is not negative, no fault
		_fault_status.flags.bad_hdg = false;
		heading_innov_var_inv = 1.0f / _heading_innov_var;

	} else {
		// the innovation variance contribution from the state covariances is negative which means the covariance matrix is badly conditioned
		_fault_status.flags.bad_hdg = true;

		// we reinitialise the covariance matrix and abort this fusion step
		initialiseCovariance();
		ECL_ERR_TIMESTAMPED("mag yaw fusion numerical error - covariance reset");
		return;
	}

	// calculate the Kalman gains
	// only calculate gains for states we are using
	Vector24f Kfusion;

	for (uint8_t row = 0; row <= 15; row++) {
		for (uint8_t col = 0; col <= 3; col++) {
			Kfusion(row) += P(row,col) * yaw_jacobian[col];
		}

		Kfusion(row) *= heading_innov_var_inv;
	}

	if (_control_status.flags.wind) {
		for (uint8_t row = 22; row <= 23; row++) {
			for (uint8_t col = 0; col <= 3; col++) {
				Kfusion(row) += P(row,col) * yaw_jacobian[col];
			}

			Kfusion(row) *= heading_innov_var_inv;
		}
	}

	// innovation test ratio
	_yaw_test_ratio = sq(innovation) / (sq(gate_sigma) * _heading_innov_var);

	// we are no longer using 3-axis fusion so set the reported test levels to zero
	_mag_test_ratio.setZero();

	// set the magnetometer unhealthy if the test fails
	if (_yaw_test_ratio > 1.0f) {
		_innov_check_fail_status.flags.reject_yaw = true;

		// if we are in air we don't want to fuse the measurement
		// we allow to use it when on the ground because the large innovation could be caused
		// by interference or a large initial gyro bias
		if (_control_status.flags.in_air) {
			return;

		} else {
			// constrain the innovation to the maximum set by the gate
			float gate_limit = sqrtf((sq(gate_sigma) * _heading_innov_var));
			_heading_innov = math::constrain(innovation, -gate_limit, gate_limit);
		}

	} else {
		_innov_check_fail_status.flags.reject_yaw = false;
		_heading_innov = innovation;
	}

	// apply covariance correction via P_new = (I -K*H)*P
	// first calculate expression for KHP
	// then calculate P - KHP
	SquareMatrix24f KHP;
	float KH[4];

	for (unsigned row = 0; row < _k_num_states; row++) {

		KH[0] = Kfusion(row) * yaw_jacobian[0];
		KH[1] = Kfusion(row) * yaw_jacobian[1];
		KH[2] = Kfusion(row) * yaw_jacobian[2];
		KH[3] = Kfusion(row) * yaw_jacobian[3];

		for (unsigned column = 0; column < _k_num_states; column++) {
			float tmp = KH[0] * P(0,column);
			tmp += KH[1] * P(1,column);
			tmp += KH[2] * P(2,column);
			tmp += KH[3] * P(3,column);
			KHP(row,column) = tmp;
		}
	}

	const bool healthy = checkAndFixCovarianceUpdate(KHP);

	_fault_status.flags.bad_hdg = !healthy;

	if (healthy) {
		// apply the covariance corrections
		P -= KHP;

		fixCovarianceErrors(true);

		// apply the state corrections
		fuse(Kfusion, _heading_innov);

	}
}

void Ekf::fuseHeading()
{

	Vector3f mag_earth_pred;
	float measured_hdg;
	float predicted_hdg;

	// Calculate the observation variance
	float R_YAW;
	if (_control_status.flags.mag_hdg) {
		// using magnetic heading tuning parameter
		R_YAW = sq(_params.mag_heading_noise);

	} else if (_control_status.flags.ev_yaw) {
		// using error estimate from external vision data
		R_YAW = _ev_sample_delayed.angVar;

	} else {
		// default value
		R_YAW = 0.01f;
	}

	// update transformation matrix from body to world frame using the current state estimate
	_R_to_earth = Dcmf(_state.quat_nominal);

	// determine if a 321 or 312 Euler sequence is best
	if (fabsf(_R_to_earth(2, 0)) < fabsf(_R_to_earth(2, 1))) {
		// rolled more than pitched so use 321 rotation order to calculate the observed yaw angle
		Eulerf euler321(_state.quat_nominal);
		predicted_hdg = euler321(2);
		if (_control_status.flags.mag_hdg) {
			// Set the yaw angle to zero and rotate the measurements into earth frame using the zero yaw angle
			euler321(2) = 0.0f;
			const Dcmf R_to_earth(euler321);
			if (_control_status.flags.mag_3D) {
				// don't apply bias corrections if we are learning them
				mag_earth_pred = R_to_earth * _mag_sample_delayed.mag;

			} else {
				mag_earth_pred = R_to_earth * (_mag_sample_delayed.mag - _state.mag_B);
			}

			// the angle of the projection onto the horizontal gives the yaw angle
			measured_hdg = -atan2f(mag_earth_pred(1), mag_earth_pred(0)) + getMagDeclination();

		} else if (_control_status.flags.ev_yaw) {
			// calculate the yaw angle for a 321 sequence
			// Expressions obtained from yaw_input_321.c produced by https://github.com/PX4/ecl/blob/master/matlab/scripts/Inertial%20Nav%20EKF/quat2yaw321.m
			const float Tbn_1_0 = 2.0f*(_ev_sample_delayed.quat(0)*_ev_sample_delayed.quat(3)+_ev_sample_delayed.quat(1)*_ev_sample_delayed.quat(2));
			const float Tbn_0_0 = sq(_ev_sample_delayed.quat(0))+sq(_ev_sample_delayed.quat(1))-sq(_ev_sample_delayed.quat(2))-sq(_ev_sample_delayed.quat(3));
			measured_hdg = atan2f(Tbn_1_0,Tbn_0_0);

		} else {
			measured_hdg = predicted_hdg;
		}

		// handle special case where yaw measurement is unavailable
		bool fuse_zero_innov = false;
		if (_is_yaw_fusion_inhibited) {
			// The yaw measurement cannot be trusted but we need to fuse something to prevent a badly
			// conditioned covariance matrix developing over time.
			if (!_control_status.flags.vehicle_at_rest) {
				// Vehicle is not at rest so fuse a zero innovation if necessary to prevent
				// unconstrained quaternion variance growth and record the predicted heading
				// to use as an observation when movement ceases.
				// TODO a better way of determining when this is necessary
				const float sumQuatVar = P(0,0) + P(1,1) + P(2,2) + P(3,3);
				if (sumQuatVar > _params.quat_max_variance) {
					fuse_zero_innov = true;
					R_YAW = 0.25f;

				}
				_last_static_yaw = predicted_hdg;

			} else {
				// Vehicle is at rest so use the last moving prediction as an observation
				// to prevent the heading from drifting and to enable yaw gyro bias learning
				// before takeoff.
				measured_hdg = _last_static_yaw;

			}
		} else {
			_last_static_yaw = predicted_hdg;

		}

		fuseYaw321(measured_hdg, R_YAW, fuse_zero_innov);

	} else {

		// pitched more than rolled so use 312 rotation order to calculate the observed yaw angle
		predicted_hdg = atan2f(-_R_to_earth(0, 1), _R_to_earth(1, 1));
		if (_control_status.flags.mag_hdg) {

			// Calculate the body to earth frame rotation matrix from the euler angles using a 312 rotation sequence
			// with yaw angle set to to zero
			const Vector3f rotVec312(0.0f,  // yaw
						 asinf(_R_to_earth(2, 1)),  // roll
						 atan2f(-_R_to_earth(2, 0), _R_to_earth(2, 2)));  // pitch
			const Dcmf R_to_earth = taitBryan312ToRotMat(rotVec312);

			// rotate the magnetometer measurements into earth frame using a zero yaw angle
			if (_control_status.flags.mag_3D) {
				// don't apply bias corrections if we are learning them
				mag_earth_pred = R_to_earth * _mag_sample_delayed.mag;
			} else {
				mag_earth_pred = R_to_earth * (_mag_sample_delayed.mag - _state.mag_B);
			}

			// the angle of the projection onto the horizontal gives the yaw angle
			measured_hdg = -atan2f(mag_earth_pred(1), mag_earth_pred(0)) + getMagDeclination();

		} else if (_control_status.flags.ev_yaw) {
			// calculate the yaw angle for a 312 sequence
			// Values from yaw_input_312.c file produced by https://github.com/PX4/ecl/blob/master/matlab/scripts/Inertial%20Nav%20EKF/quat2yaw312.m
			float Tbn_0_1_neg = 2.0f*(_ev_sample_delayed.quat(0)*_ev_sample_delayed.quat(3)-_ev_sample_delayed.quat(1)*_ev_sample_delayed.quat(2));
			float Tbn_1_1 = sq(_ev_sample_delayed.quat(0))-sq(_ev_sample_delayed.quat(1))+sq(_ev_sample_delayed.quat(2))-sq(_ev_sample_delayed.quat(3));
			measured_hdg = atan2f(Tbn_0_1_neg,Tbn_1_1);

		} else {
			measured_hdg = predicted_hdg;
		}

		// handle special case where yaw measurement is unavailable
		bool fuse_zero_innov = false;
		if (_is_yaw_fusion_inhibited) {
			// The yaw measurement cannot be trusted but we need to fuse something to prevent a badly
			// conditioned covariance matrix developing over time.
			if (!_control_status.flags.vehicle_at_rest) {
				// Vehicle is not at rest so fuse a zero innovation if necessary to prevent
				// unconstrained quaterniion variance growth and record the predicted heading
				// to use as an observation when movement ceases.
				// TODO a better way of determining when this is necessary
				const float sumQuatVar = P(0,0) + P(1,1) + P(2,2) + P(3,3);
				if (sumQuatVar > _params.quat_max_variance) {
					fuse_zero_innov = true;
					R_YAW = 0.25f;

				}
				_last_static_yaw = predicted_hdg;

			} else {
				// Vehicle is at rest so use the last moving prediction as an observation
				// to prevent the heading from drifting and to enable yaw gyro bias learning
				// before takeoff.
				measured_hdg = _last_static_yaw;

			}
		} else {
			_last_static_yaw = predicted_hdg;

		}

		fuseYaw312(measured_hdg, R_YAW, fuse_zero_innov);

	}
}

void Ekf::fuseDeclination(float decl_sigma)
{
	// assign intermediate state variables
	const float &magN = _state.mag_I(0);
	const float &magE = _state.mag_I(1);

	// minimum North field strength before calculation becomes badly conditioned (T)
	const float N_field_min = 0.001f;

	// observation variance (rad**2)
	const float R_DECL = sq(decl_sigma);

	// Calculate intermediate variables
	const float magN_sq = sq(magN);
	if (magN_sq < sq(N_field_min)) {
		// calculation is badly conditioned close to +-90 deg declination
		return;
	}
	const float HK0 = 1.0F / magN_sq;
	const float HK1 = HK0*sq(magE) + 1.0F;
	const float HK2 = 1.0F/HK1;
	const float HK3 = 1.0F/magN;
	const float HK4 = HK2*HK3;
	const float HK5 = HK3*magE;
	const float HK6 = HK5*P(16,17) - P(17,17);
	const float HK7 = 1.0F/sq(HK1);
	const float HK8 = HK5*P(16,16) - P(16,17);
	const float innovation_variance = -HK0*HK6*HK7 + HK7*HK8*magE/(magN * magN_sq) + R_DECL;
	float HK9;
	if (innovation_variance > R_DECL) {
		HK9 = HK4/innovation_variance;
	} else {
		// variance calculation is badly conditioned
		return;
	}

	// Calculate the observation Jacobian
	// Note only 2 terms are non-zero which can be used in matrix operations for calculation of Kalman gains and covariance update to significantly reduce cost
	// Note Hfusion indices do not match state indices
	float Hfusion[2] = {};
	Hfusion[0] = -HK0*HK2*magE; 	// state index 16
	Hfusion[1] = HK4;				// state index 17

	// Calculate the Kalman gains
	Vector24f Kfusion;
	for (unsigned row = 0; row <= 15; row++) {
		Kfusion(row) = -HK9*(HK5*P(row,16) - P(row,17));
	}

	Kfusion(16) = -HK8*HK9;
	Kfusion(17) = -HK6*HK9;

	for (unsigned row = 18; row <= 23; row++) {
		Kfusion(row) = -HK9*(HK5*P(16,row) - P(17,row));
	}

	const float innovation = math::constrain(atan2f(magE, magN) - getMagDeclination(), -0.5f, 0.5f);

	// apply covariance correction via P_new = (I -K*H)*P
	// first calculate expression for KHP
	// then calculate P - KHP
	// take advantage of the empty columns in KH to reduce the number of operations
	SquareMatrix24f KHP;
	float KH[2];
	for (unsigned row = 0; row < _k_num_states; row++) {

		KH[0] = Kfusion(row) * Hfusion[0];
		KH[1] = Kfusion(row) * Hfusion[1];

		for (unsigned column = 0; column < _k_num_states; column++) {
			float tmp = KH[0] * P(16,column);
			tmp += KH[1] * P(17,column);
			KHP(row,column) = tmp;
		}
	}

	const bool healthy = checkAndFixCovarianceUpdate(KHP);

	_fault_status.flags.bad_mag_decl = !healthy;

	if (healthy) {
		// apply the covariance corrections
		P -= KHP;

		fixCovarianceErrors(true);

		// apply the state corrections
		fuse(Kfusion, innovation);

		// constrain the declination of the earth field states
		limitDeclination();
	}
}

void Ekf::limitDeclination()
{
	// get a reference value for the earth field declinaton and minimum plausible horizontal field strength
	// set to 50% of the horizontal strength from geo tables if location is known
	float decl_reference;
	float h_field_min = 0.001f;
	if (_params.mag_declination_source & MASK_USE_GEO_DECL) {
		// use parameter value until GPS is available, then use value returned by geo library
		if (_NED_origin_initialised) {
			decl_reference = _mag_declination_gps;
			h_field_min = fmaxf(h_field_min , 0.5f * _mag_strength_gps * cosf(_mag_inclination_gps));
		} else {
			decl_reference = math::radians(_params.mag_declination_deg);
		}
	} else {
		// always use the parameter value
		decl_reference = math::radians(_params.mag_declination_deg);
	}

	// do not allow the horizontal field length to collapse - this will make the declination fusion badly conditioned
	// and can result in a reversal of the NE field states which the filter cannot recover from
	// apply a circular limit
	float h_field = sqrtf(_state.mag_I(0)*_state.mag_I(0) + _state.mag_I(1)*_state.mag_I(1));
	if (h_field < h_field_min) {
		if (h_field > 0.001f * h_field_min) {
			float h_scaler = h_field_min / h_field;
			_state.mag_I(0) *= h_scaler;
			_state.mag_I(1) *= h_scaler;
		} else {
			// too small to scale radially so set to expected value
			float mag_declination = getMagDeclination();
			_state.mag_I(0) = 2.0f * h_field_min * cosf(mag_declination);
			_state.mag_I(1) = 2.0f * h_field_min * sinf(mag_declination);
		}
		h_field = h_field_min;
	}

	// do not allow the declination estimate to vary too much relative to the reference value
	constexpr float decl_tolerance = 0.5f;
	const float decl_max = decl_reference + decl_tolerance;
	const float decl_min = decl_reference - decl_tolerance;
	const float decl_estimate = atan2f(_state.mag_I(1) , _state.mag_I(0));
	if (decl_estimate > decl_max)  {
		_state.mag_I(0) = h_field * cosf(decl_max);
		_state.mag_I(1) = h_field * sinf(decl_max);
	} else if (decl_estimate < decl_min)  {
		_state.mag_I(0) = h_field * cosf(decl_min);
		_state.mag_I(1) = h_field * sinf(decl_min);
	}
}

float Ekf::calculate_synthetic_mag_z_measurement(const Vector3f& mag_meas, const Vector3f& mag_earth_predicted)
{
	// theoretical magnitude of the magnetometer Z component value given X and Y sensor measurement and our knowledge
	// of the earth magnetic field vector at the current location
	const float mag_z_abs = sqrtf(math::max(sq(mag_earth_predicted.length()) - sq(mag_meas(0)) - sq(mag_meas(1)), 0.0f));

	// calculate sign of synthetic magnetomter Z component based on the sign of the predicted magnetomer Z component
	const float mag_z_body_pred = mag_earth_predicted.dot(_R_to_earth.slice<3,1>(0,2));

	return mag_z_body_pred < 0 ? -mag_z_abs : mag_z_abs;
}
