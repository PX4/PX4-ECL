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
#include "../data_validator.h"
#include "../data_validator_group.h"


const uint32_t timeout_usec = 2000;//from original private value
const int equal_value_count = 100; //default is private VALUE_EQUAL_COUNT_DEFAULT
const uint64_t base_timestamp = 666;
const unsigned base_num_siblings = 4;

/**
 * Initialize a DataValidatorGroup with some common settings;
 * @param group_handle
 * @param last_validator_handle
 * @param sibling_count
 */
void setup_group(DataValidatorGroup **group_handle,
        DataValidator **last_validator_handle,
        unsigned *sibling_count)
{
    unsigned num_siblings = base_num_siblings;

    DataValidatorGroup *group = new DataValidatorGroup(num_siblings);
    assert(nullptr != group);
    //verify that calling print doesn't crash the tests
    group->print();

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
    group->set_timeout(timeout_usec);
    //the following sets the threshold on all CURRENT members of the group, but not any added later //TODO BUG?
    group->set_equal_value_threshold(equal_value_count);

    //dynamically add a validator to the group after constructor
    DataValidator *validator = group->add_new_validator();
    //verify the previously set timeout applies to the new group member
    assert(validator->get_timeout() == timeout_usec);
    //for testing purposes, ensure this newly added member is consistent with the rest of the group
    validator->set_equal_value_threshold(equal_value_count);
    num_siblings++;

    //return values
    *group_handle = group;
    *last_validator_handle = validator;
    *sibling_count = num_siblings;

}

void test_init()
{
    DataValidatorGroup *group = nullptr;
    DataValidator *validator = nullptr;
    unsigned num_siblings = 0;

    setup_group(&group,
                &validator,
                &num_siblings);

    //should not yet be any best value
    int best_index = -1;
    assert(nullptr == group->get_best(base_timestamp, &best_index));

    delete group; //force cleanup
}

void test_put()
{
    DataValidatorGroup *group = nullptr;
    DataValidator *validator = nullptr;
    unsigned num_siblings = 0;

    setup_group(&group,
                &validator,
                &num_siblings);


    delete group; //force cleanup
}

int main(int argc, char *argv[])
{
    (void)argc; // unused
    (void)argv; // unused

    test_init();

    return 0; //passed
}
