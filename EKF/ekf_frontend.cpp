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
 * @author Siddharth B Purohit, 3DRobotics Inc. <siddharthbharatpurohit@gmail.com>
 *
 */

#include "ekf_frontend.h"

Ekf::Ekf():
_instances(0)
{}

Ekf::~Ekf()
{
	Ekf_list *temp;
	for(uint8_t i=0; i < _instances; i++) {
		temp = ekf_list_start->next_instance;
		if(temp == NULL) {
			return;
		}
		delete ekf_list_start;
		ekf_list_start = temp;
	}
	delete ekf_list_start;
}

//creates a new instance of Ekf returns -1 if failed
int8_t Ekf::create()
{
	int8_t ret = -1;
	if(_instances == 0) {
		ekf_list_start = new Ekf_list;
		if( ekf_list_start == NULL) {
			goto end;
		}
		ekf_list_end = ekf_list_start;
		_instances++;
	} else {
		ekf_list_end->next_instance = new Ekf_list;
		if(ekf_list_end->next_instance == NULL) {
			goto end;
		}
		ekf_list_end = ekf_list_end->next_instance;
		_instances++;
	}
	ret = _instances;
end:
	return ret;
}

Ekf_core* Ekf::instance(uint8_t instance)
{
	Ekf_list* ekf;
	ekf = ekf_list_start;
	for(uint8_t i=0; i < instance; i++) {
		ekf = ekf->next_instance;
		if(ekf == NULL) {
			return NULL;
		}
	}
	return &ekf->instance;
}
