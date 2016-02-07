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
 * @author Paul Riseborough <p_riseborough@live.com.au>
 *
 */

#include "estimator_interface.h"

#define sq(_arg)    powf(_arg, 2.0f)

class Ekf : public EstimatorInterface
{
public:

    Ekf();
    ~Ekf();

    bool init(uint64_t timestamp);
    bool update();

    /**
     * \brief Gets the innovations of velocity and position measurements.
     * \arg vel_pos_innov Reference array to return the velocity and position innovation values: 0-2 velocity, 3-5 position.
     */
    void get_vel_pos_innov(float vel_pos_innov[6]);

    /**
     * \brief Gets the innovations of the earth magnetic field measurements.
     * \arg mag_innov Reference array to return the magnetometer innovation values.
     */
    void get_mag_innov(float mag_innov[3]);

    /**
     * \brief Gets the innovations of the heading measurement.
     * \arg heading_innov Reference array to return the heading inovation values.
     */
    void get_heading_innov(float *heading_innov);

    /**
     * \brief Gets the innovation variances of velocity and position measurements.
     * \arg vel_pos_innov_var Reference array to return the velocity and position innovation variances: 0-2 velocity, 3-5 position.
     */
    void get_vel_pos_innov_var(float vel_pos_innov_var[6]);

    /**
     * \brief Gets the innovation variances of the earth magnetic field measurements.
     * \arg mag_innov_var Reference array to return the magnetometer innovation variances.
     */
    void get_mag_innov_var(float mag_innov_var[3]);

    /**
     * \brief Gets the innovation variance of the heading measurement.
     * \arg heading_innov_var Pointer to the heading innovation variance.
     */
    void get_heading_innov_var(float *heading_innov_var);
    
    /**
     * \brief Gets the innovation variance of the heading measurement.
     * \arg state           - @TODO - Requires description.
     */
    void get_state_delayed(float *state);

    /**
     * \brief Gets the innovation variance of the heading measurement.
     * \arg covariances Point to the measurement and noise covariance values.   @TODO - Requires description.
     */
    void get_covariances(float *covariances);

    /**
     * \brief Get the ekf WGS-84 origin positoin and height and the system time it was last set.
     * \arg origin_time     - @TODO - Requires description.
     * \arg origin_pos      - @TODO - Requires description.
     * /arg origin_alt      - @TODO - Requires description.
     */
    void get_ekf_origin(uint64_t *origin_time, map_projection_reference_s *origin_pos, float *origin_alt);

    /**
     * \brief Asks estimator for sensor data collection decision and do any preprocessing if required, returns true if not defined.
     * \arg time_usec Collects the GPS sensor update data.
     * \arg gps Pointer to the GPS data member variable.
     * \return Returns true iff successful.
     */
    bool collect_gps(uint64_t time_usec, struct gps_message *gps);

    /**
     * \brief Collects the current IMU sensor data measurements.
     * \arg imu Current imu measurement sample data.
     * \return Returns true iff successful.
     */
    bool collect_imu(imuSample &imu);

    /**
     * \param       - @TODO - Requires description.
     */
    filter_control_status_u _control_status = {};

private:

    /** \param _k_num_states Length of the State Vector, (number of states to be estimated). */
    static const uint8_t _k_num_states = 24;

    /** \param _k_earth_rate Rotational rate of the earth (rad/s). */
    static constexpr float _k_earth_rate = 0.000072921f;

    /** \param _state       - @TODO - Requires description. */
    stateSample _state;

    /** \param _filter_initialised Boolean flag to indicate if the filter has been initialized. */
    bool _filter_initialised;

    /** \param _earth_rate_initialised Boolean flag to indicate if the earth rotational rate has been initialized. */
    bool _earth_rate_initialised;

    /** \param _fuse_height Boolean flag to indicate if barometric pressure altimeter measurement should be fused. */
    bool _fuse_height;
    
    /** \param _fuse_height Boolean flag to indicate if GPS position measurement should be fused. */
    bool _fuse_pos;
    
    /** \param _fuse_height Boolean flag to indicate if GPS horizontal velocity measurement should be fused. */
    bool _fuse_hor_vel;
    
    /** \param _fuse_height Boolean flag to indicate if GPS vertical velocity measurement should be fused. */
    bool _fuse_vert_vel;

    /** \param _time_last_fake_gps Clock time of the last faked GPS measurement (usec). */
    uint64_t _time_last_fake_gps;
    
    /** \param _time_last_pos_fuse Clock time the last fusion of horizotal position measurements was performed (usec). */
    uint64_t _time_last_pos_fuse;
    
    /** \param _time_last_vel_fuse Clock time the last fusion of velocity measurements was performed (usec). */
    uint64_t _time_last_vel_fuse;
    
    /** \param _time_last_hgt_fuse Clock time the last fusion of height measurements was performed (usec). */
    uint64_t _time_last_hgt_fuse;
    
    /** \param _time_last_of_fuse Clock time the last fusion of optical flow measurements were performed (usec). */
    uint64_t _time_last_of_fuse;
    
    /** \param _last_known_posNE Last known local NE position vector (m). */
    Vector2f _last_known_posNE;
    
    /** \param _last_disarmed_posD Vertical position recorded at arming (m). */
    float _last_disarmed_posD;

    /** \param _earth_rate_NED Earth rotational rate in the local NED coordinate frame. */
    Vector3f _earth_rate_NED;

    /** \param _R_prev              - @TODO - Requires description. */
    matrix::Dcm<float> _R_prev;

    /** \param P The NxN state covariance matrix. */
    float P[_k_num_states][_k_num_states];

    /** \param _vel_pos_innov Velocity and Position Innovations: 0-2 velocity, 3-5 position. */
    float _vel_pos_innov[6];
    
    /** \param _mag_innov Earth magnetic field innovations. */
    float _mag_innov[3];
    
    /** \param _heading_innov Heading measurement innovation. */
    float _heading_innov;

    /** \param _vel_pos_innov_var Velocity and Position Innovation variances: 0-2 velocity, 3-5 position. */
    float _vel_pos_innov_var[6];
    
    /** \param _mag_innov_var Earth magnetic field innovation variance. */
    float _mag_innov_var[3];
    
    /** \param _heading_innov_var Heading measurement innovation variance. */
    float _heading_innov_var;

    //----------------------------------------- Complementary Filter States -----------------------------------------//
    /** \param _delta_angle_corr        - @TODO - Requires description. */
    Vector3f _delta_angle_corr;
    
    /** \param _delta_vel_corr          - @TODO - Requires description. */
    Vector3f _delta_vel_corr;
    
    /** \param _vel_corr                - @TODO - Requires description. */
    Vector3f _vel_corr;
    
    /** \param _imu_down_sampled        - @TODO - Requires description. */
    imuSample _imu_down_sampled;
    
    /** \param _q_down_sampled          - @TODO - Requires description. */
    Quaternion _q_down_sampled;

    //--------------------------------- Variables used for the GPS quality checks -----------------------------------//
    /** \param _gpsDriftVelN GPS north position derivative (m/s). */
    float _gpsDriftVelN = 0.0f;
    
    /** \param _gpsDriftVelE GPS east position derivative (m/s). */
    float _gpsDriftVelE = 0.0f;
    
    /** \param _gps_drift_velD GPS down position derivative (m/s). */
    float _gps_drift_velD = 0.0f;
    
    /** \param _gps_velD_diff_filt GPS filtered Down velocity (m/s). */
    float _gps_velD_diff_filt = 0.0f;
    
    /** \param _gps_velN_filt GPS filtered North velocity (m/s). */
    float _gps_velN_filt = 0.0f;
    
    /** \param _gps_velE_filt GPS filtered East velocity (m/s). */
    float _gps_velE_filt = 0.0f;
    
    /** \param _last_gps_fail_us The last system time in usec that the GPS failed it's checks. */
    uint64_t _last_gps_fail_us = 0;


    //------------------ Variables used to publish the WGS-84 location of the EKF local NED origin ------------------//
    /** \param _last_gps_origin_time_us Clock time the origin was last set (uSec). */
    uint64_t _last_gps_origin_time_us = 0;
    
    /** \param _gps_alt_ref WGS-84 height (m). */
    float _gps_alt_ref = 0.0f;

    /** \param _gps_check_fail_status           - @TODO - Requires description. */
    gps_check_fail_status_u _gps_check_fail_status;

    /** \brief Calculates the Output States.*/
    void calculateOutputStates();

    /** \brief Initializes the EKF.*/
    bool initialiseFilter(void);

    /** \brief Initializes the Covariance matrices.*/
    void initialiseCovariance();

    /** \brief State Prediction step.*/
    void predictState();

    /** \brief Covariance Matrix Prediction step.*/
    void predictCovariance();

    /** \brief Fuse Magnetometer Measurements.*/
    void fuseMag();

    /** \brief Fuse Heading Measurements.*/
    void fuseHeading();

    /** \brief Fuse Airspeed measurement.*/
    void fuseAirspeed();

    /** \brief Fuse Distance Sensor/Range measurement.*/
    void fuseRange();

    /** \brief Fuse the Velocity, Position, and Height measurements.*/
    void fuseVelPosHeight();

    /** \brief Reset the velocity estimate.*/
    void resetVelocity();

    /** \brief Reset the Position estimate.*/
    void resetPosition();

    /** \brief Enforce Covariance Matrix symmetry.*/
    void makeCovSymetrical();

    /** \brief Bound value growth in the Covariance Matrix values.*/
    void limitCov();

    /**
     * \brief Output Covariance Matrix to file.
     * \arg filename The output filename.
     */
    void printCovToFile(char const *filename);

    /** \brief                  - @TODO - Requires description.*/
    void assertCovNiceness();

    /** \brief Enforce matrix symmetry.*/
    void makeSymmetrical();

    /** \brief                  - @TODO - Requires description.*/
    void constrainStates();

    /**
     * \brief               - @TODO - Requires description.
     * \arg K               - @TODO - Requires description.
     * \arg innovation      - @TODO - Requires description.
     */
    void fuse(float *K, float innovation);

    /** \brief Print State values to terminal.*/
    void printStates();

    /** \brief Print State values at high rate to terminal.*/
    void printStatesFast();

    /**
     * \brief Calculates the Earth Rotational Rate in the NED local coordinate frame.
     * \arg omega Rotational rate vector reference to return calculated value.
     * \arg lat_rad The Current Lattitude at which to calculate the local coordinate frame rotational rate.
     */
    void calcEarthRateNED(Vector3f &omega, double lat_rad) const;

    /**
     * \brief Return true id the GPS quality is good enough to set an origin and start aiding.
     * \arg gps Pointer to the GPS message object.
     */
    bool gps_is_good(struct gps_message *gps);

    /** \brief Control the filter fusion modes. */
    void controlFusionModes();

    /** \brief Determine if we are airborne or motors are armed. */
    void calculateVehicleStatus();
};
