/****************************************************************************
 *
 *   Copyright (c) 2013 Estimation and Control Library (ECL). All rights reserved.
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
 * @file ecl_l1_pos_controller.h
 * Implementation of L1 position control.
 * Authors and acknowledgements in header.
 *
 */

#include <float.h>

#include "ecl_l1_pos_controller.h"

using matrix::Vector2f;
using matrix::wrap_pi;


void ECL_L1_Pos_Controller::update_roll_setpoint()
{
	float roll_new = atanf(_lateral_accel * 1.0f / CONSTANTS_ONE_G);
	roll_new = math::constrain(roll_new, -_roll_lim_rad, _roll_lim_rad);

	// no slew rate limiting active
	if (_dt <= 0.0f || _roll_slew_rate <= 0.0f) {
		_roll_setpoint = roll_new;
		return;
	}

	// slew rate limiting active
	roll_new = math::constrain(roll_new, _roll_setpoint - _roll_slew_rate * _dt, _roll_setpoint + _roll_slew_rate * _dt);

	_roll_setpoint = roll_new;

}

float ECL_L1_Pos_Controller::switch_distance(float wp_radius)
{
	/* following [2], switching on L1 distance */
	return math::min(wp_radius, _L1_distance);
}

float ECL_L1_Pos_Controller::airspeed_incr(float airspeed, const float airsp_ref, const float airsp_max,
		const float wind_speed, float min_ground_speed)
{
	/* following [3] */
	// NOTE: uses current "lambda" angle, this function should only be called *after any other navigation function

	const float airsp_incr_max = math::max(airsp_max - airsp_ref, 0.0f);

	airspeed = math::max(airspeed, 0.1f);
	min_ground_speed = math::max(min_ground_speed, 0.0f);

	const float wind_ratio = (wind_speed + min_ground_speed) / airspeed;
	const float feas = get_bearing_feasibility(_lambda, wind_ratio, 0.5f / airspeed);

	/* increase airspeed reference to match wind speed (or minimum ground speed) when nav bearing is infeasible, as airspeed max allows */
	return math::constrain(wind_speed + min_ground_speed - airsp_ref, 0.0f, airsp_incr_max) * (1.0f - feas);
}

void ECL_L1_Pos_Controller::navigate_waypoints(const Vector2f &vector_A, const Vector2f &vector_B,
		const Vector2f &vector_curr_position, const Vector2f &ground_speed_vector, const Vector2f &wind_speed_vector)
{
	/* this follows logic presented in [1], [2], and [3] */

	/* reset _L1_ratio */
	set_l1_period(_L1_period);

	/* get the direction between the last (visited) and next waypoint */
	_target_bearing = get_bearing_to_next_waypoint((double)vector_curr_position(0), (double)vector_curr_position(1),
			  (double)vector_B(0), (double)vector_B(1));

	/* enforce a minimum ground speed/airspeed of 0.1 m/s to avoid singularities */
	const float ground_speed = math::max(ground_speed_vector.length(), 0.1f);
	Vector2f airspeed_vector = ground_speed_vector - wind_speed_vector;
	const float airspeed = math::max(airspeed_vector.length(), 0.1f);

	/* calculate the L1 length required for the desired period */
	_L1_distance = _L1_ratio * ground_speed;

	/* calculate vector from A to B */
	Vector2f vector_AB = get_local_planar_vector(vector_A, vector_B);

	/*
	 * check if waypoints are on top of each other. If yes,
	 * skip A and directly continue to B
	 */
	if (vector_AB.length() < 1.0e-6f) {
		vector_AB = get_local_planar_vector(vector_curr_position, vector_B);
	}

	vector_AB.normalize();

	/* calculate the vector from waypoint A to the aircraft */
	Vector2f vector_A_to_airplane = get_local_planar_vector(vector_A, vector_curr_position);

	/* calculate crosstrack error (output only) */
	_crosstrack_error = vector_AB % vector_A_to_airplane;

	/*
	 * If the current position is in a +-135 degree angle behind waypoint A
	 * and further away from A than the L1 distance, then A becomes the L1 point.
	 * If the aircraft is already between A and B normal L1 logic is applied.
	 */
	float distance_A_to_airplane = vector_A_to_airplane.length();
	float alongTrackDist = vector_A_to_airplane * vector_AB;

	/* estimate airplane position WRT to B */
	Vector2f vector_B_to_P_unit = get_local_planar_vector(vector_B, vector_curr_position).normalized();

	/* calculate angle of airplane position vector relative to line) */

	// XXX this could probably also be based solely on the dot product
	float AB_to_BP_bearing = atan2f(vector_B_to_P_unit % vector_AB, vector_B_to_P_unit * vector_AB);

	/* extension from [2], fly directly to A */
	if (distance_A_to_airplane > _L1_distance && alongTrackDist / math::max(distance_A_to_airplane, 1.0f) < -0.7071f) {

		/* calculate eta to fly to waypoint A */

		/* bearing from current position to L1 point */
		_nav_bearing = atan2f(-vector_A_to_airplane(1), -vector_A_to_airplane(0));

		/*
		 * If the AB vector and the vector from B to airplane point in the same
		 * direction, we have missed the waypoint. At +- 90 degrees we are just passing it.
		 */

	} else if (fabsf(AB_to_BP_bearing) < math::radians(100.0f)) {
		/*
		 * Extension, fly back to waypoint.
		 *
		 * This corner case is possible if the system was following
		 * the AB line from waypoint A to waypoint B, then is
		 * switched to manual mode (or otherwise misses the waypoint)
		 * and behind the waypoint continues to follow the AB line.
		 */

		/* calculate eta to fly to waypoint B */

		/* bearing from current position to L1 point */
		_nav_bearing = atan2f(-vector_B_to_P_unit(1), -vector_B_to_P_unit(0));

	} else {

		/* calculate eta to fly along the line between A and B */

		/* calculate eta1 (angle to L1 point) */
		float xtrackErr = vector_A_to_airplane % vector_AB;
		float sine_eta1 = xtrackErr / math::max(_L1_distance, 0.1f);

		/* limit output to 45 degrees */
		sine_eta1 = math::constrain(sine_eta1, -0.7071f, 0.7071f); //sin(pi/4) = 0.7071

		/* bearing from current position to L1 point */
		_nav_bearing = atan2f(vector_AB(1), vector_AB(0)) + asinf(sine_eta1);
	}

	/* L1 (nav) vector */
	Vector2f nav_vector {cosf(_nav_bearing), sinf(_nav_bearing)};

	/* angle "lambda" between wind and look-ahead vectors, in [-pi,pi] */
	_lambda = atan2f(wind_speed_vector % nav_vector, wind_speed_vector * nav_vector);

	/* wind speed to airspeed ratio */
	const float wind_ratio = wind_speed_vector.length() / airspeed;

	/* nav bearing feasibility
	 * 0.5f / airspeed = wind ratio buffer (avoids numerical issues at feasibility boundary) */
	const float feas = get_bearing_feasibility(_lambda, wind_ratio, 0.5f / airspeed);

	/* transition nav vector from ground speed to airspeed vector in over-wind scenarios to
	 * mitigate run-away and maintain continuity in acceleration commands (no jumps!) */
	Vector2f nav_speed {feas * ground_speed_vector(0) + (1.0f - feas) *airspeed_vector(0), feas * ground_speed_vector(1) + (1.0f - feas) *airspeed_vector(1)};

	/* calculate error angle (limit to 90 degrees) */
	const float nav_speed_course = atan2f(nav_speed(1), nav_speed(0));
	float eta = wrap_pi(_nav_bearing - nav_speed_course);
	eta = math::constrain(eta, -M_PI_2_F, M_PI_2_F);

	/* demanded lateral acceleration */
	_lateral_accel = _K_L1 * nav_speed.length() / _L1_ratio * sinf(eta);

	/* flying to waypoints, not circling them */
	_circle_mode = false;

	/* the bearing angle, in NED frame */
	_bearing_error = eta;

	update_roll_setpoint();
}

void
ECL_L1_Pos_Controller::navigate_loiter(const Vector2f &vector_A, const Vector2f &vector_curr_position, float radius,
				       int8_t loiter_direction, const Vector2f &ground_speed_vector, const Vector2f &wind_speed_vector)
{
	/* the guidance logic in this section was proposed by [3] */

	_circle_mode = false;

	/* reset _L1_ratio */
	set_l1_period(_L1_period);

	/* update bearing to next waypoint */
	_target_bearing = get_bearing_to_next_waypoint((double)vector_curr_position(0), (double)vector_curr_position(1),
			  (double)vector_A(0), (double)vector_A(1));

	/* enforce a minimum ground speed/airspeed of 0.1 m/s to avoid singularities */
	const float ground_speed = math::max(ground_speed_vector.length(), 0.1f);
	Vector2f airspeed_vector = ground_speed_vector - wind_speed_vector;
	const float airspeed = math::max(airspeed_vector.length(), 0.1f);

	/* calculate the L1 length required for the desired period */
	_L1_distance = _L1_ratio * ground_speed;

	/* calculate the vector from waypoint A to current position */
	Vector2f vector_A_to_airplane = get_local_planar_vector(vector_A, vector_curr_position);

	/* avoid erroneous bearing calculations when aircraft is in center of loiter
	 * -- set arbitrary direction until out of center */
	if (vector_A_to_airplane.length() < 0.1f) {
		vector_A_to_airplane(0) = 0.1;
		vector_A_to_airplane(1) = 0.0;
	}

	const float norm_vector_A_to_airplane = vector_A_to_airplane.length();

	/* radial distance from the loiter circle (not center) */
	const float xtrack_err_circle = norm_vector_A_to_airplane - radius;

	/* check circle tracking feasibility */
	if (_L1_distance > radius && fabsf(xtrack_err_circle) <= _L1_distance) {
		// re-calculate L1 ratio & distance (effectively reduce period to allow tracking small radii)
		// NOTE: when far from the loiter radius, operator defined period and damping are unaffected
		_L1_ratio = math::max(fabsf(xtrack_err_circle), radius) / ground_speed;
		_L1_distance = _L1_ratio * ground_speed;
	}

	/* calculate L1 (nav) bearing
	 * use law of cosines to calculate the angle "gamma" between the vector from the aircraft to the
	 * circle center and the L1 vector (from aircraft to point on the circle) */
	float cos_gamma = (_L1_distance * _L1_distance + norm_vector_A_to_airplane * norm_vector_A_to_airplane - radius *
			   radius) * 0.5f / _L1_distance / norm_vector_A_to_airplane;
	// constrain arccosine input to safe domain -- also produces perpendicular approach to path when outside L1 distance
	cos_gamma = math::constrain(cos_gamma, -1.0f, 1.0f);
	// rotate bearing from aircraft to A by gamma in loiter direction to calculate nav bearing
	const float bearing_airplane_to_A = atan2f(-vector_A_to_airplane(1), -vector_A_to_airplane(0));
	_nav_bearing = wrap_pi(bearing_airplane_to_A - (float)loiter_direction * acosf(cos_gamma));

	/* L1 (nav) vector */
	Vector2f nav_vector {cosf(_nav_bearing), sinf(_nav_bearing)};

	/* angle "lambda" between wind and look-ahead vectors, in [-pi,pi] */
	_lambda = atan2f(wind_speed_vector % nav_vector, wind_speed_vector * nav_vector);

	/* wind speed to airspeed ratio */
	const float wind_ratio = wind_speed_vector.length() / airspeed;

	/* nav bearing feasibility
	 * 0.5f / airspeed = wind ratio buffer (avoids numerical issues at feasibility boundary) */
	const float feas = get_bearing_feasibility(_lambda, wind_ratio, 0.5f / airspeed);

	/* transition nav vector from ground speed to airspeed vector in over-wind scenarios to
	 * mitigate run-away and maintain continuity in acceleration commands (no jumps!) */
	Vector2f nav_speed {feas * ground_speed_vector(0) + (1.0f - feas) *airspeed_vector(0), feas * ground_speed_vector(1) + (1.0f - feas) *airspeed_vector(1)};

	/* calculate error angle (limit to 90 degrees) */
	const float nav_speed_course = atan2f(nav_speed(1), nav_speed(0));
	float eta = wrap_pi(_nav_bearing - nav_speed_course);
	eta = math::constrain(eta, -M_PI_2_F, M_PI_2_F);

	/* feedback */
	_crosstrack_error = xtrack_err_circle;
	_lateral_accel =  _K_L1 * nav_speed.length() / _L1_ratio * sinf(eta);
	_bearing_error = eta;

	if (fabsf(cos_gamma) < 1.0f) { _circle_mode = true; }

	update_roll_setpoint();
}

void ECL_L1_Pos_Controller::navigate_heading(float navigation_heading, float current_heading,
		const Vector2f &ground_speed_vector, const Vector2f &wind_speed_vector)
{
	/* the complete guidance logic in this section was proposed by [2] */

	/*
	 * As the commanded heading is the only reference
	 * (and no crosstrack correction occurs),
	 * target and navigation bearing become the same
	 */
	_target_bearing = _nav_bearing = wrap_pi(navigation_heading);

	/* L1 (nav) vector */
	Vector2f nav_vector {cosf(_nav_bearing), sinf(_nav_bearing)};

	/* angle "lambda" between wind and look-ahead vectors, in [-pi,pi] */
	_lambda = atan2f(wind_speed_vector % nav_vector, wind_speed_vector * nav_vector);

	float eta = wrap_pi(_target_bearing - wrap_pi(current_heading));

	/* consequently the bearing error is exactly eta: */
	_bearing_error = eta;

	/* ground speed is the length of the ground speed vector */
	float ground_speed = ground_speed_vector.length();

	/* adjust L1 distance to keep constant frequency */
	_L1_distance = ground_speed / _heading_omega;
	float omega_vel = ground_speed * _heading_omega;

	/* not circling a waypoint */
	_circle_mode = false;

	/* navigating heading means by definition no crosstrack error */
	_crosstrack_error = 0;

	/* limit eta to 90 degrees */
	eta = math::constrain(eta, -M_PI_2_F, M_PI_2_F);
	_lateral_accel = 2.0f * sinf(eta) * omega_vel;

	update_roll_setpoint();
}

void ECL_L1_Pos_Controller::navigate_level_flight(float current_heading)
{
	/* the logic in this section is trivial, but originally proposed by [2] */

	/* reset all heading / error measures resulting in zero roll */
	_target_bearing = current_heading;
	_nav_bearing = current_heading;
	_bearing_error = 0;
	_crosstrack_error = 0;
	_lateral_accel = 0;

	/* not circling a waypoint when flying level */
	_circle_mode = false;

	update_roll_setpoint();
}

Vector2f ECL_L1_Pos_Controller::get_local_planar_vector(const Vector2f &origin, const Vector2f &target) const
{
	/* this is an approximation for small angles, proposed by [2] */
	Vector2f out(math::radians((target(0) - origin(0))),
		     math::radians((target(1) - origin(1))*cosf(math::radians(origin(0)))));

	return out * static_cast<float>(CONSTANTS_RADIUS_OF_EARTH);
}

void ECL_L1_Pos_Controller::set_l1_period(float period)
{
	_L1_period = period;

	/* calculate the ratio introduced in [2] */
	_L1_ratio = 1.0f / M_PI_F * _L1_damping * _L1_period;

	/* calculate normalized frequency for heading tracking */
	_heading_omega = sqrtf(2.0f) * M_PI_F / _L1_period;
}

void ECL_L1_Pos_Controller::set_l1_damping(float damping)
{
	_L1_damping = damping;

	/* calculate the ratio introduced in [2] */
	_L1_ratio = 1.0f / M_PI_F * _L1_damping * _L1_period;

	/* calculate the L1 gain (following [2]) */
	_K_L1 = 4.0f * _L1_damping * _L1_damping;
}

float ECL_L1_Pos_Controller::get_bearing_feasibility(float lambda, const float wind_ratio, const float wind_ratio_buf)
{
	/* following [3] */

	/* bound lambda -- angle between wind and bearing */
	lambda = math::constrain(fabsf(lambda), 0.0f, M_PI_2_F);

	/* upper and lower feasibility barriers */
	float wind_ratio_ub;
	float wind_ratio_lb;

	if (lambda < LAMBDA_CO) {
		/* linear finite cut-off */
		const float mx = M_CO * (LAMBDA_CO - lambda);
		const float wind_ratio_ub_co = ONE_OVER_S_LAMBDA_CO;
		wind_ratio_ub = wind_ratio_ub_co + mx;
		const float wind_ratio_lb_co = (ONE_OVER_S_LAMBDA_CO - 2.0f) * wind_ratio_buf + 1.0f;
		wind_ratio_lb = wind_ratio_lb_co + wind_ratio_buf * mx;

	} else {
		const float one_over_s_lambda = 1.0f / sinf(lambda);
		wind_ratio_ub = one_over_s_lambda;
		wind_ratio_lb = (one_over_s_lambda - 2.0f) * wind_ratio_buf + 1.0f;
	}

	/* calculate bearing feasibility */
	float feas;

	if (wind_ratio > wind_ratio_ub) {
		// infeasible
		feas = 0.0f;

	} else if (wind_ratio > wind_ratio_lb) {
		// partially feasible
		// smoothly transition from fully feasible to fully infeasible
		feas = 0.5f + 0.5f * cosf(2.0f * M_PI_2_F * math::constrain((wind_ratio - wind_ratio_lb) /
					  (wind_ratio_ub - wind_ratio_lb), 0.0f, 1.0f)); // cos^2(x) = 0.5+0.5cos(2x)

	} else {
		// feasible
		feas = 1.0f;
	}

	return feas;
}
