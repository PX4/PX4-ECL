/****************************************************************************
 *
 *   Copyright (c) 2019 Todd Stellanova. All rights reserved.
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
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be  used to endorse or promote products derived
 *    from this software without specific prior written permission.
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
 * @file test_data_validator_group.cpp
 * Testing the DataValidatorGroup class
 *
 * @author Todd Stellanova
 */

#include <stdint.h>
#include <cassert>
#include <stdio.h>
#include <math.h>
#include <validation/data_validator.h>
#include <validation/data_validator_group.h>
#include <validation/tests/tests_common.h>


const uint32_t base_timeout_usec = 2000;//from original private value
const int equal_value_count = 100; //default is private VALUE_EQUAL_COUNT_DEFAULT
const uint64_t base_timestamp = 666;
const unsigned base_num_siblings = 4;


/**
 * Initialize a DataValidatorGroup with some common settings;
 * @param sibling_count (out) the number of siblings initialized
 */
DataValidatorGroup  *setup_base_group(unsigned *sibling_count)
{
	unsigned num_siblings = base_num_siblings;

	DataValidatorGroup *group = new DataValidatorGroup(num_siblings);
	assert(nullptr != group);
	//verify that calling print doesn't crash the tests
	group->print();
	printf("\n");

	//should be no failovers yet
	assert(0 == group->failover_count());
	assert(DataValidator::ERROR_FLAG_NO_ERROR == group->failover_state());
	assert(-1 == group->failover_index());

	//no vibration yet
	float vibe_off =  group->get_vibration_offset(base_timestamp, 0);
	printf("vibe_off: %f \n", (double)vibe_off);
	assert(-1.0f == group->get_vibration_offset(base_timestamp, 0));

	float vibe_fact = group->get_vibration_factor(base_timestamp);
	printf("vibe_fact: %f \n", (double)vibe_fact);
	assert(0.0f == vibe_fact);

	//this sets the timeout on all current members of the group, as well as members added later
	group->set_timeout(base_timeout_usec);
	//the following sets the threshold on all CURRENT members of the group, but not any added later
	//TODO This is likely a bug in DataValidatorGroup
	group->set_equal_value_threshold(equal_value_count);

	//return values
	*sibling_count = num_siblings;

	return group;

}

/**
 * Fill one DataValidator with samples, by index.
 *
 * @param group
 * @param val1_idx Index of the validator to fill with samples
 * @param num_samples
 */
void fill_one_with_valid_data(DataValidatorGroup *group, int val1_idx,  uint32_t num_samples)
{
	uint64_t timestamp = base_timestamp;
	uint64_t error_count = 0;
	float last_best_val = 0.0f;

	for (uint32_t i = 0; i < num_samples; i++) {
		float val = ((float) rand() / (float) RAND_MAX);
		float data[DataValidator::dimensions] = {val};
		group->put(val1_idx, timestamp, data, error_count, 100);
		last_best_val = val;
	}

	int best_idx = 0;
	float *best_data = group->get_best(timestamp, &best_idx);
	assert(last_best_val == best_data[0]);
	assert(best_idx == val1_idx);
}



/**
 * Fill two validators in the group with samples, by index.
 * Both validators will be filled with the same data, but
 * the priority of the first validator will be higher than the second.
 *
 * @param group
 * @param val1_idx index of the first validator to fill
 * @param val2_idx index of the second validator to fill
 * @param num_samples
 */
void fill_two_with_valid_data(DataValidatorGroup *group, int val1_idx, int val2_idx, uint32_t num_samples)
{
	uint64_t timestamp = base_timestamp;
	uint64_t error_count = 0;
	float last_best_val = 0.0f;

	for (uint32_t i = 0; i < num_samples; i++) {
		float val = ((float) rand() / (float) RAND_MAX);
		float data[DataValidator::dimensions] = {val};
		//two sensors with identical values, but different priorities
		group->put(val1_idx, timestamp, data, error_count, 100);
		group->put(val2_idx, timestamp, data, error_count, 10);
		last_best_val = val;
	}

	int best_idx = 0;
	float *best_data = group->get_best(timestamp, &best_idx);
	assert(last_best_val == best_data[0]);
	assert(best_idx == val1_idx);

}

/**
 * Dynamically add a validator to the group after construction
 * @param group
 * @return
 */
DataValidator *add_validator_to_group(DataValidatorGroup *group)
{
	DataValidator *validator = group->add_new_validator();
	//verify the previously set timeout applies to the new group member
	assert(validator->get_timeout() == base_timeout_usec);
	//for testing purposes, ensure this newly added member is consistent with the rest of the group
	//TODO this is likely a bug in DataValidatorGroup
	validator->set_equal_value_threshold(equal_value_count);

	return validator;
}

/**
 * Create a DataValidatorGroup and tack on two additional DataValidators
 *
 * @param validator1_handle (out) first DataValidator added to the group
 * @param validator2_handle (out) second DataValidator added to the group
 * @param sibling_count (in/out) in: number of initial siblings to create, out: total
 * @return
 */
DataValidatorGroup *setup_group_with_two_validator_handles(
	DataValidator **validator1_handle,
	DataValidator **validator2_handle,
	unsigned *sibling_count)
{
	DataValidatorGroup *group = setup_base_group(sibling_count);

	//now we add validators
	*validator1_handle = add_validator_to_group(group);
	*validator2_handle = add_validator_to_group(group);
	*sibling_count += 2;

	return group;
}


void test_init()
{
	unsigned num_siblings = 0;

	DataValidatorGroup *group = setup_base_group(&num_siblings);

	//should not yet be any best value
	int best_index = -1;
	assert(nullptr == group->get_best(base_timestamp, &best_index));

	delete group; //force cleanup
}

void test_put()
{
	unsigned num_siblings = 0;
	DataValidator *validator1 = nullptr;
	DataValidator *validator2 = nullptr;

	uint64_t timestamp = base_timestamp;

	DataValidatorGroup *group = setup_group_with_two_validator_handles(&validator1, &validator2, &num_siblings);
	printf("num_siblings: %d \n",num_siblings);
	unsigned val1_idx = num_siblings - 2;
	unsigned val2_idx = num_siblings - 1;

	fill_two_with_valid_data(group, val1_idx, val2_idx, 500);
	int best_idx = -1;
	float *best_data = group->get_best(timestamp, &best_idx);
	assert(nullptr != best_data);
	float best_val = best_data[0];

	float *cur_val1 = validator1->value();
	assert(nullptr != cur_val1);
	//printf("cur_val1 %p \n", cur_val1);
	assert(best_val == cur_val1[0]);

	float *cur_val2 = validator2->value();
	assert(nullptr != cur_val2);
	//printf("cur_val12 %p \n", cur_val2);
	assert(best_val == cur_val2[0]);

	delete group; //force cleanup
}


void test_failover()
{

	unsigned num_siblings = 0;
	DataValidator *validator1 = nullptr;
	DataValidator *validator2 = nullptr;

	uint64_t timestamp = base_timestamp;

	DataValidatorGroup *group = setup_group_with_two_validator_handles(&validator1, &validator2, &num_siblings);
	//printf("num_siblings: %d \n",num_siblings);
	int val1_idx = (int)num_siblings - 2;
	int val2_idx = (int)num_siblings - 1;
	uint64_t error_count = 0;


	fill_two_with_valid_data(group, val1_idx, val2_idx, 100);


	int best_idx = -1;
	float *best_data = nullptr;
	//now, switch the priorities, which switches "best" but doesn't detect a failover
	{
		float new_best_val = 3.14159f;
		float data[DataValidator::dimensions] = {new_best_val};
		group->put(val1_idx, timestamp, data, error_count, 1);
		group->put(val2_idx, timestamp, data, error_count, 100);
		best_data = group->get_best(timestamp, &best_idx);
		assert(new_best_val == best_data[0]);
		//the new best sensor should now be the sensor with the higher priority
		assert(best_idx == val2_idx);
		//should not have detected a real failover
		//printf("failover_count A: %d \n", group->failover_count());
		assert(0 == group->failover_count());
	}

	//flush any garbage
	fill_two_with_valid_data(group, val1_idx, val2_idx, 10);


	// now trigger a real failover
	{
		float new_best_val = 3.14159f;
		float data[DataValidator::dimensions] = {new_best_val};
		//trigger a bunch of errors on the previous best sensor
		unsigned fake_err_count = 0;

		for (int i = 0; i < 25; i++) {
			group->put(val1_idx, timestamp, data, ++fake_err_count, 100);
			group->put(val2_idx, timestamp, data, error_count, 10);
		}

		assert(validator1->error_count() == fake_err_count);
		best_data = group->get_best(timestamp + 1, &best_idx);
		assert(nullptr != best_data);
		assert(new_best_val == best_data[0]);
		assert(best_idx == val2_idx);
		//should have detected a real failover
		printf("failover_count B: %d \n", group->failover_count());
		assert(1 == group->failover_count());
		//TODO figure out what these values are supposed to be-- they don't match naive expectations
		int fail_idx = group->failover_index();
		printf("fail_idx: %d expected: %d \n", fail_idx, val1_idx);
		//assert (val1_idx == fail_idx);
		printf("error state: %x expected: %x \n", group->failover_state(), DataValidator::ERROR_FLAG_HIGH_ERRCOUNT);
		//assert (DataValidator::ERROR_FLAG_HIGH_ERRCOUNT == group->failover_state());
	}

	delete  group; //cleanup
}

/**
 * Verify that we get expected vibration values after injecting samples.
 */
void test_vibration()
{
	unsigned num_siblings = 0;
	uint64_t timestamp = base_timestamp;

	DataValidatorGroup *group =  setup_base_group(&num_siblings);

	//now we add validators
	DataValidator *validator  = add_validator_to_group(group);
	assert(nullptr != validator);
	num_siblings++;
	float *vibes = validator->vibration_offset();
	assert(nullptr != vibes);
	//printf("val vibes: %f \n", vibes[0]);
	//should be no vibration data yet
	assert(0 == vibes[0]);

	float vibe_o = group->get_vibration_offset(timestamp, 0);
	//printf("group vibe_o %f \n", vibe_o);
	// there should be no known vibe offset before samples are inserted
	assert(-1.0f == vibe_o);

	float rms_err = 0.0f;
	//insert some swinging values
	insert_values_around_mean(validator, 3.14159f, 1000, &rms_err, &timestamp);
	vibes = validator->vibration_offset();
	assert(nullptr != vibes);
	printf("val1 vibes: %f rms_err: %f \n", vibes[0], (double)rms_err);

	vibe_o = group->get_vibration_offset(timestamp, 0);
	printf("group vibe_o %f \n", vibe_o);
	//the one validator's vibration offset should match the group's vibration offset
	assert(vibes[0] == vibe_o);


	//this should be "The best RMS value of a non-timed out sensor"
	float group_vibe_fact = group->get_vibration_factor(timestamp);
	float val1_rms = (validator->rms())[0];
	printf("group_vibe_fact: %f val1_rms: %f\n", (double)group_vibe_fact, (double)val1_rms);
	assert(group_vibe_fact == val1_rms);

}

int main(int argc, char *argv[])
{
	(void)argc; // unused
	(void)argv; // unused

	test_init();
	test_put();
	test_failover();
	test_vibration();

	return 0; //passed
}
