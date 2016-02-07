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
 * @file estimator_base.h
 * Definition of base class for attitude estimators
 *
 * @author Roman Bast <bapstroman@gmail.com>
 * @author Siddharth Bharat Purohit <siddharthbharatpurohit@gmail.com>
 *
 */

namespace estimator
{
    struct gps_message {
        /** \param time_usec Clock time in mircoseconds. */
        uint64_t time_usec;

        /** \param lat Latitude in 1E-7 degrees. */
        int32_t lat;
        
        /** \param lon Longitude in 1E-7 degrees. */
        int32_t lon;
        
        /** \param alt Altitude in 1E-3 meters (millimeters) above MSL. */
        int32_t alt;
        
        /** \param fix_type Indicator Values: 0-1: no fix, 2: 2D fix, 3: 3D fix, 4: RTCM code differential, 5: Real-Time. */
        uint8_t fix_type;
        
        /** \param eph GPS horizontal position accuracy in m. */
        float eph;
        
        /** \param epv GPS vertical position accuracy in m.. */
        float epv;
        
        /** \param sacc GPS speed accuracy in m/s. */
        float sacc;
        
        /** \param time_usec_vel Timestamp for velocity informations. */
        uint64_t time_usec_vel;
        
        /** \param vel_m_s  GPS ground speed (m/s). */
        float vel_m_s;
        
        /** \param vel_ned GPS ground speed NED. */
        float vel_ned[3];
        
        /** \param vel_ned_valid GPS ground speed is valid. */
        bool vel_ned_valid;
        
        /** \param nsats Number of satellites used. */
        uint8_t nsats;
        
        /** \param gdop Geometric dilution of precision. */
        float gdop;
    };


    /** Type definitions specific to the Estimation Control library. */
    typedef matrix::Vector<float, 2> Vector2f;
    typedef matrix::Vector<float, 3> Vector3f;
    typedef matrix::Quaternion<float> Quaternion;
    typedef matrix::Matrix<float, 3, 3> Matrix3f;

    /**
     * \brief Data Structure for Position, Velocity, Attitude output from the Estimator in the local coordinate frame.
     */
    struct outputSample {
        /** \param quat_nominal Attitude quaternion relating the current body attitdue to the reference frame. */
        Quaternion  quat_nominal;

        /** \param vel Body velocity vector in the inertial reference frame. */
        Vector3f    vel;

        /** \param pos Position vector of the body in the inertial reference NED coordiante frame. */
        Vector3f    pos;

        /** \param time_us Timestamp associated with the estimated values (usec). */
        uint64_t    time_us;
    };

    /**
     * \brief Data Structure for IMU Measurement Sample data.
     */
    struct imuSample {

        /** \param delta_ang        - @TODO - Requires description. */
        Vector3f    delta_ang;
        
        /** \param delta_vel        - @TODO - Requires description.. */
        Vector3f    delta_vel;
        
        /** \param delta_ang_dt     - @TODO - Requires description.. */
        float       delta_ang_dt;
        
        /** \param delta_vel_dt     - @TODO - Requires description.. */
        float       delta_vel_dt;
        
        /** \param time_us Clock Time associated with the measurement sample (usec). */
        uint64_t    time_us;
    };

    /**
     * \brief Data Structure for GPS Measurement Sample data.
     */
    struct gpsSample {

        /** \param pos GPS Position estimate in Lat/Lon (deg)   - @TODO - Are assumptions of units here correct? */
        Vector2f    pos;

        /** \param hgt GPS Altitude measurement (m). */
        float       hgt;

        /** \param vel GPS Velocity measurement (m/s). */
        Vector3f    vel;

        /** \param time_us Clock Time associated with the measurement sample (usec). */
        uint64_t    time_us;
    };

    /**
     * \brief Data Structure for Magnetometer Sample data.
     */
    struct magSample {

        /** \param mag Vector of Magnetometer measurement sample (Gauss). */
        Vector3f    mag;

        /** \param time_us Clock Time associated with the measurement sample (usec). */
        uint64_t    time_us;
    };

    /**
     * \brief Data Structure for Barometric Pressure Sensor Sample data.
     */
    struct baroSample {

        /** \param hgt Barometric Altimeter absolute pressure reading converted to meters (m). */
        float       hgt;

        /** \param time_us Clock Time associated with the measurement sample (usec). */
        uint64_t    time_us;
    };

    /**
     * \brief Data Structure for Distance Sensor Sample data.
     */
    struct rangeSample {

        /** \param rng The Range sensor measurement sample converted to meters (m). */
        float       rng;

        /** \param time_us Clock Time associated with the measurement sample (usec). */
        uint64_t    time_us;
    };

    /**
     * \brief Data Structure for Dynamic Pressure Sensor Sample data.
     */
    struct airspeedSample {

        /** \param airspeed Dynamic Pressure Sensor measurement converted to speed (m/s). */
        float       airspeed;

        /** \param time_us Clock Time associated with the measurement sample (usec). */
        uint64_t    time_us;
    };

    /**
     * \brief Data Structure for Optical Flow sensor Sample data.
     */
    struct flowSample {

        /** \param flowRadXY        - @TODO - Requires description. */
        Vector2f    flowRadXY;

        /** \param flowRadXYcomp    - @TODO - Requires description. */
        Vector2f    flowRadXYcomp;

        /** \param time_us Clock Time associated with the measurement sample (usec). */
        uint64_t    time_us;
    };

    /**
     * \brief State Estimation Filter Parameter Data Structure.
     */
    struct parameters {

        /** \param Magnetometer measurement delay relative to the IMU. */
        float mag_delay_ms = 0.0f;

        /** \param baro_delay_ms The barometer height measurement delay relative to the IMU. */
        float baro_delay_ms = 0.0f;

        /** \param gps_delay_ms GPS measurement delay relative to the IMU. */
        float gps_delay_ms = 200.0f;

        /** \param airspeed_delay_ms Airspeed measurement delay relative to the IMU. */
        float airspeed_delay_ms = 200.0f;

        //------------------------------------------- Input Noise Parameters --------------------------------------------//
        /** \param gyro_noise IMU angular rate noise used for covariance prediction. */
        float gyro_noise = 1.0e-3f;

        /** \param accel_noise IMU acceleration noise use for covariance prediction. */
        float accel_noise = 2.5e-1f;

        //------------------------------------------ Process Noise Parameters -------------------------------------------//
        /** \param gyro_bias_p_noise Process noise for IMU delta angle bias prediction. */
        float gyro_bias_p_noise = 7.0e-5f;
        
        /** \param accel_bias_p_noise Process noise for IMU delta velocity bias prediction. */
        float accel_bias_p_noise = 1.0e-4f;
        
        /** \param  gyro_scale_p_noise Process noise for gyro scale factor prediction. */
        float gyro_scale_p_noise = 3.0e-3f;
        
        /** \param mag_p_noise Process noise for magnetic field prediction. */
        float mag_p_noise = 2.5e-2f;
        
        /** \param wind_vel_p_noise Process noise for wind velocity prediction. */
        float wind_vel_p_noise = 1.0e-1f;

        /** \param gps_pos_noise Observation noise for gps velocity fusion. */
        float gps_vel_noise = 5.0e-1f;
        
        /** \param gps_pos_noise Observation noise for gps position fusion. */
        float gps_pos_noise = 1.0f;
        
        /** \param pos_noaid_noise Observation noise for non-aiding position fusion. */
        float pos_noaid_noise = 10.0f;
        
        /** \param baro_noise Observation noise for barometric height fusion. */
        float baro_noise = 3.0f;

        /** \param baro_innov_gate Barometric height innovation consistency gate size in standard deviations. */
        float baro_innov_gate = 3.0f;

        /** \param posNE_innov_gate GPS horizontal position innovation consistency gate size in standard deviations. */
        float posNE_innov_gate = 3.0f;

        /** \param vel_innov_gate GPS velocity innovation consistency gate size in standard deviations. */
        float vel_innov_gate = 3.0f;

        /** \param mag_heading_noise Measurement noise used for simple heading fusion. */
        float mag_heading_noise = 1.7e-1f;

        /** \param mag_noise Measurement noise used for 3-axis magnetoemeter fusion. */
        float mag_noise = 5.0e-2f;
        
        /** \param mag_declination_deg Magnetic declination in degrees. */
        float mag_declination_deg = 0.0f; 
        
        /** \param heading_innov_gate Heading fusion innovation consistency gate size in standard deviations. */
        float heading_innov_gate = 3.0f;
        
        /** \param mag_innov_gate Magnetometer fusion innovation consistency gate size in standard deviations. */
        float mag_innov_gate = 3.0f;

        // ------------------------------------------ GPS Quality Checks ------------------------------------------------//
        /**
         * \detail These parameters control the strictness of GPS quality checks used to determine uf the GPS is
         *         good enough to set a local origin and commence aiding.
         */
        
        /** \param gps_check_mask Bitmask used to control which GPS quality checks are used. */
        int gps_check_mask = 21;
        
        /** \param req_hacc Maximum acceptable horizontal position error. */
        float req_hacc = 5.0f;
        
        /** \param req_vacc Maximum acceptable vertical position error. */
        float req_vacc = 8.0f;
        
        /** \param req_sacc Maximum acceptable speed error. */
        float req_sacc = 1.0f;
        
        /** \param req_nsats Minimum acceptable satellite count. */
        int req_nsats = 6;
        
        /** \param req_gdop Maximum acceptable geometric dilution of precision*/
        float req_gdop = 2.0f;
        
        /** \param req_hdrift Maximum acceptable horizontal drift speed. */
        float req_hdrift = 0.3f;
        
        /** \param req_vdrift Maximum acceptable vertical drift speed. */
        float req_vdrift = 0.5f;
    };

    /**
     * \brief Data Structure for Sampled Sensor States.
     */
    struct stateSample {
        
        /** \param ang_error        - @TODO - Requires description. */
        Vector3f    ang_error;
        
        /** \param vel              - @TODO - Requires description. */
        Vector3f    vel;
        
        /** \param pos              - @TODO - Requires description. */
        Vector3f    pos;
        
        /** \param gyro_bias        - @TODO - Requires description. */
        Vector3f    gyro_bias;
        
        /** \param gyro_scale       - @TODO - Requires description. */
        Vector3f    gyro_scale;
        
        /** \param accel_z_bias     - @TODO - Requires description. */
        float       accel_z_bias;
        
        /** \param mag_I            - @TODO - Requires description. */
        Vector3f    mag_I;
        
        /** \param mag_B            - @TODO - Requires description. */
        Vector3f    mag_B;
        
        /** \param wind_vel         - @TODO - Requires description. */
        Vector2f    wind_vel;
        
        /** \param quat_nominal     - @TODO - Requires description. */
        Quaternion  quat_nominal;
    };

    /**
     * \brief Data Structure for Fault Status indicators.
     */
    struct fault_status_t {
        
        /** \param bad_mag_x Status Indicator of bad Magnetometer X reading. */
        bool bad_mag_x: 1;
        
        /** \param bad_mag_y Status Indicator of bad Magnetometer Y reading. */
        bool bad_mag_y: 1;
        
        /** \param bad_mag_z Status Indicator of bad Magnetometer Z reading. */
        bool bad_mag_z: 1;
        
        /** \param bad_airspeed Status Indicator of bad Dynamic Pressure Sensor reading. */
        bool bad_airspeed: 1;
        
        /** \param bad_sideslip Status Indicator of Sideslip reading (?). */
        bool bad_sideslip: 1;
    };

    /**
     * \brief Union Data Structure for publishing the status of GPS quality checks
     */
    union gps_check_fail_status_u {
        struct {
            
            /** \param fix Value: 0 - true if the fix type is insufficient, (no 3D solution). */
            uint16_t fix    : 1; //
            
            /** \param nsats Value: 1 - true if number of satellites used is insufficient. */
            uint16_t nsats  : 1;
            
            /** \param gdop Value: 2 - true if geometric dilution of precision is insufficient. */
            uint16_t gdop   : 1;
            
            /** \param hacc Value: 3 - true if reported horizontal accuracy is insufficient. */
            uint16_t hacc   : 1;
            
            /** \param vacc Value: 4 - true if reported vertical accuracy is insufficient. */
            uint16_t vacc   : 1;
            
            /** \param sacc Value: 5 - true if reported speed accuracy is insufficient. */
            uint16_t sacc   : 1;
            
            /** \param hdrift Value: 6 - true if horizontal drift is excessive (can only be used when stationary on ground). */
            uint16_t hdrift : 1; 
            
            /** \param vdrift Value: 7 - true if vertical drift is excessive (can only be used when stationary on ground). */
            uint16_t vdrift : 1;
            
            /** \param hspeed Value: 8 - true if horizontal speed is excessive (can only be used when stationary on ground). */
            uint16_t hspeed : 1;
            
            /** \param vspeed Value: 9 - true if vertical speed error is excessive. */
            uint16_t vspeed : 1;
        } flags;
        uint16_t value;
    };

    /**
     * \brief Union Data Structure Bitmask containing filter control status.
     */
    union filter_control_status_u {
        struct {
            
            /** \param angle_align Value: 0 - true if the filter angular alignment is complete. */
            uint8_t angle_align : 1;
            
            /** \param gps Value: 1 - true if GPS measurements are being fused. */
            uint8_t gps         : 1;
            
            /** \param opt_flow Value: 2 - true if optical flow measurements are being fused. */
            uint8_t opt_flow    : 1;
            
            /** \param mag_hdg Value: 3 - true if a simple magnetic heading is being fused. */
            uint8_t mag_hdg     : 1;
            
            /** \param mag_3D Value: 4 - true if 3-axis magnetometer measurement are being fused. */
            uint8_t mag_3D      : 1;
            
            /** \param mag_dec Value: 5 - true if synthetic magnetic declination measurements are being fused. */
            uint8_t mag_dec     : 1;
            
            /** \param in_air Value: 6 - true when the vehicle is airborne. */
            uint8_t in_air      : 1;
            
            /** \param armed Value: 7 - true when the vehicle motors are armed. */
            uint8_t armed       : 1;
        } flags;
        uint16_t value;
    };
}
