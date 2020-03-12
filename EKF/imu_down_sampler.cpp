#include "imu_down_sampler.hpp"

ImuDownSampler::ImuDownSampler(const int32_t& imu_samples) : _imu_samples{imu_samples} { reset(); }

// integrate imu samples until target dt reached
// assumes that dt of the gyroscope is close to the dt of the accelerometer
// returns true if target dt is reached
bool ImuDownSampler::update(const imuSample& imu_sample_new) {
	if (_do_reset) {
		reset();
	}

	_samples++;

	// accumulate time deltas
	_imu_down_sampled.delta_ang_dt += imu_sample_new.delta_ang_dt;
	_imu_down_sampled.delta_vel_dt += imu_sample_new.delta_vel_dt;
	_imu_down_sampled.time_us = imu_sample_new.time_us;

	// use a quaternion to accumulate delta angle data
	// this quaternion represents the rotation from the start to end of the accumulation period
	const Quatf delta_q(AxisAnglef(imu_sample_new.delta_ang));
	_delta_angle_accumulated = _delta_angle_accumulated * delta_q;
	_delta_angle_accumulated.normalize();

	// rotate the accumulated delta velocity data forward each time so it is always in the updated rotation frame
	const Dcmf delta_R(delta_q.inversed());
	_imu_down_sampled.delta_vel = delta_R * _imu_down_sampled.delta_vel;

	// accumulate the most recent delta velocity data at the updated rotation frame
	// assume effective sample time is halfway between the previous and current rotation frame
	_imu_down_sampled.delta_vel += (imu_sample_new.delta_vel + delta_R * imu_sample_new.delta_vel) * 0.5f;

	// check if we have enough IMU samples (and at least 2 ms of data) or if it's already more than 20 ms
	if (((_samples >= _imu_samples) && (_imu_down_sampled.delta_ang_dt >= 0.002f)) ||
	    (_imu_down_sampled.delta_ang_dt >= 0.020f)) {

		_imu_down_sampled.delta_ang = _delta_angle_accumulated.to_axis_angle();
		return true;
	} else {
		return false;
	}
}

void ImuDownSampler::reset() {
	_imu_down_sampled.delta_ang.setZero();
	_imu_down_sampled.delta_vel.setZero();
	_imu_down_sampled.delta_ang_dt = 0.0f;
	_imu_down_sampled.delta_vel_dt = 0.0f;
	_delta_angle_accumulated.setIdentity();
	_samples = 0;
	_do_reset = false;
}
