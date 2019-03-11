//
// Created by Todd Stellanova on 2019-03-11.
//

#include <stdio.h>
#include "tests_common.h"

/**
 * Insert a series of samples around a mean value
 * @param validator The validator under test
 * @param mean The mean value
 * @param count The number of samples to insert in the validator
 * @param rms_err (out) calculated rms error of the inserted samples
 */
void insert_values_around_mean(DataValidator *validator, const float mean, uint32_t count, float *rms_err)
{
    uint64_t timestamp = 500;
    uint64_t timestamp_incr = 5;
    const uint64_t error_count = 0;
    const int priority = 50;
    const float swing = 1E-2f;
    double sum_dev_squares = 0.0f;

    //insert a series of values that swing around the mean
    for (uint32_t i = 0; i < count;  i++) {
        float iter_swing = (0 == (i % 2)) ? swing : -swing;
        float iter_val = mean + iter_swing;
        float iter_dev = iter_val - mean;
        sum_dev_squares += (iter_dev * iter_dev);
        timestamp += timestamp_incr;
        validator->put(timestamp, iter_val, error_count, priority);
    }

    double rms = sqrt(sum_dev_squares / (double)count);
    //note: this should be approximately equal to "swing"
    *rms_err = (float)rms;
}

void dump_validator_state(DataValidator* validator)
{
    uint32_t state = validator->state();
    printf("state: 0x%x no_data: %d stale: %d timeout:%d\n",
           validator->state(),
           DataValidator::ERROR_FLAG_NO_DATA & state,
           DataValidator::ERROR_FLAG_STALE_DATA & state,
           DataValidator::ERROR_FLAG_TIMEOUT & state
    );
    validator->print();
}