//
// Created by Todd Stellanova on 2019-03-11.
//

#ifndef ECL_TESTS_COMMON_H
#define ECL_TESTS_COMMON_H

#include "../data_validator.h"

void insert_values_around_mean(DataValidator *validator, const float mean, uint32_t count, float *rms_err,
                               uint64_t *last_timestamp);

#endif //ECL_TESTS_COMMON_H
