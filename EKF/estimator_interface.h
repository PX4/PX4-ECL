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
 *
 */

#include <stdint.h>
#include <matrix/matrix/math.hpp>
#include <lib/geo/geo.h>
#include "RingBuffer.h"

#include "common.h"

using namespace estimator;
class EstimatorInterface
{

public:
    EstimatorInterface();
    ~EstimatorInterface();

    virtual bool init(uint64_t timestamp) = 0;
    virtual bool update() = 0;

    /**
     * \brief Gets the innovations of velocity and position measurements.
     * \arg vel_pos_innov Array of velocity and position innovation values: 0-2 velocity, 3-5 position.
     */
    virtual void get_vel_pos_innov(float vel_pos_innov[6]) = 0;

    /**
     * \brief Gets the innovations of the earth magnetic field measurements.
     * \arg mag_innov Array to store the magnetometer innovation values.
     */
    virtual void get_mag_innov(float mag_innov[3]) = 0;

    /**
     * \brief Gets the innovations of the heading measurement.
     * \arg heading_innov Array to store the heading inovation values.
     */
    virtual void get_heading_innov(float *heading_innov) = 0;

    /**
     * \brief Gets the innovation variances of velocity and position measurements.
     * \arg vel_pos_innov_var Array to store the velocity and position innovation variances: 0-2 velocity, 3-5 position.
     */
    virtual void get_vel_pos_innov_var(float vel_pos_innov_var[6]) = 0;

    /**
     * \brief Gets the innovation variances of the earth magnetic field measurements.
     * \arg mag_innov_var Array to store the magnetometer innovation variances.
     */
    virtual void get_mag_innov_var(float mag_innov_var[3]) = 0;

    /**
     * \brief Gets the innovation variance of the heading measurement.
     * \arg heading_innov_var Pointer to the heading innovation variance.
     */
    virtual void get_heading_innov_var(float *heading_innov_var) = 0;
    
    /**
     * \brief Gets the innovation variance of the heading measurement.
     * \arg state   @TODO - Requires description.
     */
    virtual void get_state_delayed(float *state) = 0;

    /**
     * \brief Gets the innovation variance of the heading measurement.
     * \arg covariances Point to the measurement and noise covariance values.   @TODO - Requires description.
     */
    virtual void get_covariances(float *covariances) = 0;

    /**
     * \brief Get the ekf WGS-84 origin positoin and height and the system time it was last set.
     * \arg origin_time     @TODO - Requires description.
     * \arg origin_pos      @TODO - Requires description.
     * /arg origin_alt      @TODO - Requires description.
     */
    virtual void get_ekf_origin(uint64_t *origin_time, map_projection_reference_s *origin_pos, float *origin_alt) = 0;

    /**
     * \brief Asks estimator for sensor data collection decision and do any preprocessing if required, returns true if not defined.
     * \arg time_usec Collects the GPS sensor update data.
     * \arg gps Pointer to the GPS data member variable.
     * \return Returns true iff successful.
     */
    virtual bool collect_gps(uint64_t time_usec, struct gps_message *gps) { return true; }

    /**
     * \brief Collects the current IMU sensor data measurements.
     * \arg imu Current imu measurement sample data.
     * \return Returns true iff successful.
     */
    virtual bool collect_imu(imuSample &imu) { return true; }

    /**
     * \brief Collects the current Magnetometer data measurements at \var time_usec.
     * \arg time_usec Current time of the sensor update.
     * \arg data Pointer to the magnetometer member variable.
     * \return Returns true iff successful.
     */
    virtual bool collect_mag(uint64_t time_usec, float *data) { return true; }

    /**
     * \brief Collects the current Barometric pressure sensor data measurements at \var time_usec.
     * \arg time_usec Current time of the sensor update.
     * \arg data Pointer to the barometric pressure sensor member variable.
     * \return Returns true iff successful.
     */
    virtual bool collect_baro(uint64_t time_usec, float *data) { return true; }

    /**
     * \brief Collects the current Dynamic pressure (Airspeed) sensor data measurements at \var time_usec.
     * \arg time_usec Current time of the sensor update.
     * \arg data Pointer to the dynamic pressure (airspeed) sensor member variable.
     * \return Returns true iff successful.
     */
    virtual bool collect_airspeed(uint64_t time_usec, float *data) { return true; }

    /**
     * \brief Collects the current Distance/Range sensor data measurements at \var time_usec.
     * \arg time_usec Current time of the sensor update.
     * \arg data Pointer to the distance/range sensor member variable.
     * \return Returns true iff successful.
     */
    virtual bool collect_range(uint64_t time_usec, float *data) { return true; }

    /**
     * \brief Collects the current Optical Flow sensor data measurements at \var time_usec.
     * \arg time_usec Current time of the sensor update.
     * \arg data Pointer to the optical flow sensor member variable.
     * \return Returns true iff successful.
     */
    virtual bool collect_opticalflow(uint64_t time_usec, float *data) { return true; }

    /**
     * \brief Sets delta angle imu data member variable at \var time_usec.
     * \arg time_usec Current time of the sensor update.
     * \arg delta_ang_dt                                                @TODO - Requires description. 
     * \arg delta_vel_dt                                                @TODO - Requires description. 
     * \arg delta_ang Pointer to the \var delta_ang member variable.    @TODO - Requires description review.
     * \arg delta_vel Pointer to the \var delta_vel member variable.    @TODO - Requires description review. 
     */
    void setIMUData(uint64_t time_usec, uint64_t delta_ang_dt, uint64_t delta_vel_dt, float *delta_ang, float *delta_vel);

    /**
     * \brief Sets magnetometer data member variable at \var time_usec.
     * \arg time_usec Current time of the sensor update.
     * \arg data Pointer to the magnetometer sensor data member variable.
     */
    void setMagData(uint64_t time_usec, float *data);

    /**
     * \brief Sets gps data member variable at \var time_usec.
     * \arg time_usec Current time of the sensor update.
     * \arg data data Pointer to the GPS data member variable.
     */
    void setGpsData(uint64_t time_usec, struct gps_message *gps);

    /**
     * \brief Sets baro data member variable at \var time_usec.
     * \arg time_usec Current time of the sensor update.
     * \arg data Pointer to the barometric pressure sensor data member variable.
     */
    void setBaroData(uint64_t time_usec, float *data);

    /**
     * \brief Sets airspeed data member variable at \var time_usec.
     * \arg time_usec Current time of the sensor update.
     * \arg data Pointer to the dynamic pressure airspeed sensor data member variable.
     */
    void setAirspeedData(uint64_t time_usec, float *data);

    /**
     * \brief Sets range data at \var time_usec.
     * \arg time_usec Current time of the sensor update.
     * \arg data Pointer to the range/distance sensor data member variable.
     */
    void setRangeData(uint64_t time_usec, float *data);

    /**
     * \brief Sets optical flow data at \var time_usec.
     * \arg time_usec Current time of the sensor update.
     * \arg data Pointer to the optical flow sensor data member variable.
     */
    void setOpticalFlowData(uint64_t time_usec, float *data);

    /**
     * \brief Returns the handle to the parameters struct in order to give access to the application.   @TODO - Requires description review. 
     */
    parameters *getParamHandle() {return &_params;}

    /**
     * \brief Sets vehicle arm status data.
     * \arg data The armed status.
     */
    void set_arm_status(bool data) { _vehicle_armed = data; }

    /**
     * \brief Prints the current IMU measurement sample data to terminal.
     * \arg data The IMU Sample data to printout.
     */
    void printIMU(struct imuSample *data);

    /**
     * \brief  Prints the Stored IMU measurement sample data to terminal.
     */
    void printStoredIMU();

    /**
     * \brief Prints the Quaternion values to terminal.
     * \arg q The Quaternion to printout.
     */
    void printQuaternion(Quaternion &q);

    /**
     * \brief Prints the average IMU sample data timing interval to terminal.
     */
    void print_imu_avg_time();

    /**
     * \brief Prints the Magnetometer measurement sample data to terminal.
     * \arg data The Magnetometer data to prinout.
     */
    void printMag(struct magSample *data);

    /**
     * \brief Prints the Stored Magnetometer data to terminal.
     */
    void printStoredMag();

    /**
     * \brief Prints the Barometric Pressure Sensor measurement sample data to terminal.
     * \arg data The barometric pressure sensor data to printout.
     */
    void printBaro(struct baroSample *data);

    /**
     * \brief Prints the Stored Barometric Pressure Sensor measurement data to terminal.
     */
    void printStoredBaro();

    /**
     * \brief Prints the GPS sample data to terminal.
     * \arg data The GPS sample data to prinout.
     */
    void printGps(struct gpsSample *data);

    /**
     * \brief Prints the Stored GPS data to terminal.
     */
    void printStoredGps();

    /**
     * \brief Returns true if the current position estimate is valid.
     * \return Returns true if the current position estimate is valid.
     */
    bool position_is_valid();


    /**
     * \brief Copies the quaternion estimation to the \var quat member variable.
     * \arg quat Pointer to the attitude quaternion estimate member variable.
     */
    void copy_quaternion(float *quat)
    {
        for (unsigned i = 0; i < 4; i++) {
            quat[i] = _output_new.quat_nominal(i);
        }
    }

    /**
     * \brief Copies the velocity estimation to the \var vel member variable.
     * \arg vel Pointer to the velocity estimate member variable.
     */
    void copy_velocity(float *vel)
    {
        for (unsigned i = 0; i < 3; i++) {
            vel[i] = _output_new.vel(i);
        }
    }

    /**
     * \brief Copies the position estimate to the \var pos member variable.
     * \arg Pointer to the position estimate member variable.
     */
    void copy_position(float *pos)
    {
        for (unsigned i = 0; i < 3; i++) {
            pos[i] = _output_new.pos(i);
        }
    }

    /**
     * \brief Copies the current timestamp value.
     * \arg time_us Pointer to the time stamp member variable.
     */
    void copy_timestamp(uint64_t *time_us)
    {
        *time_us = _imu_time_last;
    }

protected:

    /** \var _params Filter parameters. */
    parameters _params;

    /** \var OBS_BUFFER_LENGTH Buffer length. */
    static const uint8_t OBS_BUFFER_LENGTH = 10;

    /** \var IMU_BUFFER_LENGTH Buffer length. */
    static const uint8_t IMU_BUFFER_LENGTH = 30;

    /** \var FILTER_UPDATE_PERIOD_MS Update period in msec. */
    static const unsigned FILTER_UPDATE_PERIOD_MS = 10;

    /** \var _dt_imu_avg Average IMU update period in msec. */
    float _dt_imu_avg;

    /** \var _imu_time_last Time value for the previous IMU update. */
    uint64_t _imu_time_last;

    /** \var _imu_sample_delayed    @TODO - Requires description. */
    imuSample _imu_sample_delayed;

    /** \var _mag_sample_delayed    @TODO - Requires description. */
    magSample _mag_sample_delayed;

    /** \var _baro_sample_delayed   @TODO - Requires description. */
    baroSample _baro_sample_delayed;
    
    /** \var _gps_sample_delayed    @TODO - Requires description. */
    gpsSample _gps_sample_delayed;
    
    /** \var _range_sample_delayed  @TODO - Requires description. */
    rangeSample _range_sample_delayed;
    
    /** \var _airspeed_sample_delayed   @TODO - Requires description. */
    airspeedSample _airspeed_sample_delayed;
    
    /** \var _flow_sample_delayed   @TODO - Requires description. */
    flowSample _flow_sample_delayed;

    /** \var _output_sample_delayed @TODO - Requires description. */
    outputSample _output_sample_delayed;
    
    /** \var _output_new    @TODO - Requires description. */
    outputSample _output_new;
    
    /** \var _imu_sample_new    @TODO - Requires description. */
    imuSample _imu_sample_new;

    /** \var _imu_ticks @TODO - Requires description. */
    uint64_t _imu_ticks;

    /** \var _imu_updated Boolean flag to indicate if the IMU data has been updated. */
    bool _imu_updated = false;
    
    /** \var _initialised Boolean flag to indicate if the filter has been initialized.*/
    bool _initialised = false;
    
    /** \var _vehicle_armed Boolean flag to indicate if the vehicle has been armed. Used to disable functionality while on the ground. */
    bool _vehicle_armed = false;

    /** \var _NED_origin_initialised Boolean flag to indicate if the NED Origin has been initialized.*/
    bool _NED_origin_initialised = false;
    
    /** \var _gps_speed_valid Boolean flag to indicate if the GPS ground speed data has been updated.*/
    bool _gps_speed_valid = false;
    
    /** \var _gps_speed_accuracy GPS receiver reported speed accuracy (m/s). */
    float _gps_speed_accuracy = 0.0f; 
    
    /** \var _pos_ref Contains WGS-84 position latitude and longitude (radians). */
    struct map_projection_reference_s _pos_ref = {};    // 

    /** \var _mag_healthy Boolean Flag to indicate if the magnetometer estimate health. Computed by mag innovation test. */
    bool _mag_healthy = false;
    
    /** \var _yaw_test_ratio Yaw innovation consistency check ratio. */
    float _yaw_test_ratio;
    
    /** \var _mag_test_ratio Magnetometer XYZ innovation consistency check ratios. */
    float _mag_test_ratio[3];

    /** \var _vel_pos_test_ratio Velocity and Position innovation consistency check ratios*/
    float _vel_pos_test_ratio[6];

    /** \var _imu_buffer IMU Data buffer. */
    RingBuffer<imuSample> _imu_buffer;
    
    /** \var _gps_buffer GPS Data buffer. */
    RingBuffer<gpsSample> _gps_buffer;
    
    /** \var _mag_buffer Magentometer Data buffer. */
    RingBuffer<magSample> _mag_buffer;
    
    /** \var _baro_buffer Barometric Pressure Sensor Data buffer. */
    RingBuffer<baroSample> _baro_buffer;
    
    /** \var _range_buffer Range/Distance Sensor Data buffer. */
    RingBuffer<rangeSample> _range_buffer;
    
    /** \var _airspeed_buffer Dynamic Pressure/Airspeed Sensor Data buffer. */
    RingBuffer<airspeedSample> _airspeed_buffer;
    
    /** \var _flow_buffer Optical Flow Sensor Data buffer. */
    RingBuffer<flowSample>  _flow_buffer;
    
    /** \var _output_buffer Output Data buffer. @TODO - Requires description review. */
    RingBuffer<outputSample> _output_buffer;

    /** \var _time_last_imu The clock value for the most recent IMU update. */
    uint64_t _time_last_imu;
    
    /** \var _time_last_gps The clock value for the most recent GPS update. */
    uint64_t _time_last_gps;
    
    /** \var _time_last_mag The clock value for the most recent Magenetometer update. */
    uint64_t _time_last_mag;
    
    /** \var _time_last_baro The clock value for the most recent Barometric Presure sensor update. */
    uint64_t _time_last_baro;
    
    /** \var _time_last_range The clock value for the most recent range/distance sensor update. */
    uint64_t _time_last_range;
    
    /** \var _time_last_airspeed The clock value for the most recent dynamic presure (airspeed) sensor update. */
    uint64_t _time_last_airspeed;

    /** \var _fault_status Indicates a fault. @TODO - Requires description review. */
    fault_status_t _fault_status;
    
    /**
     * \brief Initializes The ECL Interface member variables. @TODO - Requires description review.
     * \arg timestamp The current time to demarcate initializion instance.
     */
    bool initialise_interface(uint64_t timestamp);
    
    /** 
     * \brief Cleanup/Deallocate all sensor buffer memory. @TODO - Requires description review.
     */
    void unallocate_buffers();
};
