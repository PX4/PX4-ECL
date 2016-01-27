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
 * @file ekf.h
 * Class for core functions for ekf attitude and position estimator.
 *
 * @author Roman Bast <bapstroman@gmail.com>
 *
 */

#include "estimator_base.h"

#define sq(_arg)	powf(_arg, 2.0f)

class Ekf : public EstimatorBase
{
public:

	Ekf();
	~Ekf();

	// should be called every time new data is pushed into the filter
	bool update();

	// gets the innovations of velocity and position measurements
	// 0-2 vel, 3-5 pos
	void get_vel_pos_innov(float vel_pos_innov[6]);

	// gets the innovations of the earth magnetic field measurements
	void get_mag_innov(float mag_innov[3]);

	// gets the innovations of the heading measurement
	void get_heading_innov(float *heading_innov);

	// gets the innovation variances of velocity and position measurements
	// 0-2 vel, 3-5 pos
	void get_vel_pos_innov_var(float vel_pos_innov_var[6]);

	// gets the innovation variances of the earth magnetic field measurements
	void get_mag_innov_var(float mag_innov_var[3]);

	// gets the innovation variance of the heading measurement
	void get_heading_innov_var(float *heading_innov_var);

	// get the state vector at the delayed time horizon
	void get_state_delayed(float *state);

	// get the diagonal elements of the covariance matrix
	void get_covariances(float *covariances);

private:

	static const uint8_t _k_num_states = 24;
	static constexpr float _k_earth_rate = 0.000072921f;

	bool _filter_initialised;
	bool _earth_rate_initialised;

	bool _fuse_height;	// baro height data should be fused
	bool _fuse_pos;		// gps position data should be fused
	bool _fuse_hor_vel;		// gps horizontal velocity measurement should be fused
	bool _fuse_vert_vel;	// gps vertical velocity measurement should be fused

	uint8_t _mag_fuse_index;	// counter for sequential mag axis fusion

	uint64_t _time_last_fake_gps;	// last timestamp in microseconds when we faked a gps measurement for static mode 

	Vector3f _earth_rate_NED;		// earth rotation vector (NED) in rad/s

	matrix::Dcm<float> _R_prev;		// transformation matrix from earth frame to body frame of previous ekf step

	float P[_k_num_states][_k_num_states];	// state covariance matrix

	float _vel_pos_innov[6];	// innovations: 0-2 vel,  3-5 pos
	float _mag_innov[3];		// earth magnetic field innovations
	float _heading_innov;		// heading measurement innovation

	float _vel_pos_innov_var[6]; // innovation variances: 0-2 vel, 3-5 pos
	float _mag_innov_var[3]; 	// earth magnetic field innovation variance
	float _heading_innov_var; // heading measurement innovation variance

	// complementary filter states
	Vector3f _delta_angle_corr;	// delta angle correction vector
	Vector3f _delta_vel_corr;	// delta velocity correction vector
	Vector3f _vel_corr;			// velocity correction vector

	// update the real time complementary filter states. This includes the prediction
	// and the correction step
	void calculateOutputStates();

	// initialise filter states of both the delayed ekf and the real time complementary filter
	bool initialiseFilter(void);

	// initialise ekf covariance matrix
	void initialiseCovariance();

	// predict ekf state
	void predictState();

	// predict ekf covariance
	void predictCovariance();

	// ekf sequential fusion of magnetometer measurements. index {0,1,2} defines the axis
	// which should be fused
	void fuseMag(uint8_t index);

	// fuse magnetometer heading measurement
	void fuseHeading();

	// fuse airspeed measurement
	void fuseAirspeed();

	// fuse range measurements
	void fuseRange();

	// fuse velocity and position measurements (also barometer height)
	void fuseVelPosHeight();

	// reset velocity states of the ekf
	void resetVelocity();

	// reset position state of the ekf
	void resetPosition();

	// limit the diagonal of the covariance matrix
	void limitCov();

	// make ekf covariance matrix symmetric
	void makeSymmetrical();

	// constrain the ekf states
	void constrainStates();

	// generic function which will perform a fusion step given a kalman gain
	// and a scalar innovation value
	void fuse(float *K, float innovation);

	// calculate the earth rotation vector from a given latitude
	void calcEarthRateNED(Vector3f &omega, double lat_rad) const;

};
