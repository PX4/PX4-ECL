#include <math.h>
#include <stdio.h>
#include <cstdlib>
#include "../../../../matrix/matrix/math.hpp"

typedef matrix::Vector<float, 24> Vector24f;
typedef matrix::SquareMatrix<float, 24> SquareMatrix24f;

float sq(float in) {
    return in * in;
}

int main()
{
    // Compare calculation of observation Jacobians and Kalman gains for sympy and matlab generated equations

	float Hfusion[24];
    Vector24f H_TAS;
    Vector24f Kfusion;
    float airspeed_innov_var;

	const float R_TAS = sq(1.5f);

    const bool update_wind_only = false;

	// get latest velocity in earth frame
	const float vn = 9.0f;
	const float ve = 12.0f;
	const float vd = -1.5f;

	// get latest wind velocity in earth frame
	const float vwn = -4.0f;
	const float vwe = 3.0f;

    // create a symmetrical positive dfinite matrix with off diagonals between -1 and 1 and diagonals between 0 and 1
    SquareMatrix24f P;
    for (int col=0; col<=23; col++) {
        for (int row=0; row<=col; row++) {
            if (row == col) {
                P(row,col) = (float)rand();
            } else {
                P(col,row) = P(row,col) = 2.0f * ((float)rand() - 0.5f);
            }
        }
    }

    // First calculate observationjacobians and Kalman gains using sympy generated equations

    // Intermediate variables
    const float HK0 = vn - vwn;
    const float HK1 = ve - vwe;
    const float HK2 = powf(HK0, 2) + powf(HK1, 2) + powf(vd, 2);
    const float HK3 = powf(HK2, -1.0F/2.0F);
    const float HK4 = HK0*HK3;
    const float HK5 = HK1*HK3;
    const float HK6 = 1.0F/HK2;
    const float HK7 = HK0*P(4,6) - HK0*P(6,22) + HK1*P(5,6) - HK1*P(6,23) + P(6,6)*vd;
    const float HK8 = HK1*P(5,23);
    const float HK9 = HK0*P(4,5) - HK0*P(5,22) + HK1*P(5,5) - HK8 + P(5,6)*vd;
    const float HK10 = HK1*HK6;
    const float HK11 = HK0*P(4,22);
    const float HK12 = HK0*P(4,4) - HK1*P(4,23) + HK1*P(4,5) - HK11 + P(4,6)*vd;
    const float HK13 = HK0*HK6;
    const float HK14 = -HK0*P(22,23) + HK0*P(4,23) - HK1*P(23,23) + HK8 + P(6,23)*vd;
    const float HK15 = -HK0*P(22,22) - HK1*P(22,23) + HK1*P(5,22) + HK11 + P(6,22)*vd;
    const float HK16 = HK3/(-HK10*HK14 + HK10*HK9 + HK12*HK13 - HK13*HK15 + HK6*HK7*vd + R_TAS);

	// Observation Jacobians
	memset(&Hfusion, 0, sizeof(Hfusion));
	Hfusion[4] = HK4;    // corresponds to state index 4
	Hfusion[5] = HK5;    // corresponds to state index 5
	Hfusion[6] = HK3*vd; // corresponds to state index 6
	Hfusion[22] = -HK4;   // corresponds to state index 22
	Hfusion[23] = -HK5;   // corresponds to state index 23

    Kfusion(0) = HK16*(-HK0*P(0,22) + HK0*P(0,4) - HK1*P(0,23) + HK1*P(0,5) + P(0,6)*vd);
    Kfusion(1) = HK16*(-HK0*P(1,22) + HK0*P(1,4) - HK1*P(1,23) + HK1*P(1,5) + P(1,6)*vd);
    Kfusion(2) = HK16*(-HK0*P(2,22) + HK0*P(2,4) - HK1*P(2,23) + HK1*P(2,5) + P(2,6)*vd);
    Kfusion(3) = HK16*(-HK0*P(3,22) + HK0*P(3,4) - HK1*P(3,23) + HK1*P(3,5) + P(3,6)*vd);
    Kfusion(4) = HK12*HK16;
    Kfusion(5) = HK16*HK9;
    Kfusion(6) = HK16*HK7;
    Kfusion(7) = HK16*(HK0*P(4,7) - HK0*P(7,22) + HK1*P(5,7) - HK1*P(7,23) + P(6,7)*vd);
    Kfusion(8) = HK16*(HK0*P(4,8) - HK0*P(8,22) + HK1*P(5,8) - HK1*P(8,23) + P(6,8)*vd);
    Kfusion(9) = HK16*(HK0*P(4,9) - HK0*P(9,22) + HK1*P(5,9) - HK1*P(9,23) + P(6,9)*vd);
    Kfusion(10) = HK16*(-HK0*P(10,22) + HK0*P(4,10) - HK1*P(10,23) + HK1*P(5,10) + P(6,10)*vd);
    Kfusion(11) = HK16*(-HK0*P(11,22) + HK0*P(4,11) - HK1*P(11,23) + HK1*P(5,11) + P(6,11)*vd);
    Kfusion(12) = HK16*(-HK0*P(12,22) + HK0*P(4,12) - HK1*P(12,23) + HK1*P(5,12) + P(6,12)*vd);
    Kfusion(13) = HK16*(-HK0*P(13,22) + HK0*P(4,13) - HK1*P(13,23) + HK1*P(5,13) + P(6,13)*vd);
    Kfusion(14) = HK16*(-HK0*P(14,22) + HK0*P(4,14) - HK1*P(14,23) + HK1*P(5,14) + P(6,14)*vd);
    Kfusion(15) = HK16*(-HK0*P(15,22) + HK0*P(4,15) - HK1*P(15,23) + HK1*P(5,15) + P(6,15)*vd);
    Kfusion(16) = HK16*(-HK0*P(16,22) + HK0*P(4,16) - HK1*P(16,23) + HK1*P(5,16) + P(6,16)*vd);
    Kfusion(17) = HK16*(-HK0*P(17,22) + HK0*P(4,17) - HK1*P(17,23) + HK1*P(5,17) + P(6,17)*vd);
    Kfusion(18) = HK16*(-HK0*P(18,22) + HK0*P(4,18) - HK1*P(18,23) + HK1*P(5,18) + P(6,18)*vd);
    Kfusion(19) = HK16*(-HK0*P(19,22) + HK0*P(4,19) - HK1*P(19,23) + HK1*P(5,19) + P(6,19)*vd);
    Kfusion(20) = HK16*(-HK0*P(20,22) + HK0*P(4,20) - HK1*P(20,23) + HK1*P(5,20) + P(6,20)*vd);
    Kfusion(21) = HK16*(-HK0*P(21,22) + HK0*P(4,21) - HK1*P(21,23) + HK1*P(5,21) + P(6,21)*vd);
    Kfusion(22) = HK15*HK16;
    Kfusion(23) = HK14*HK16;

    // save output and repeat calculation using legacy matlab generated code
    float Hfusion_sympy[24];
    Vector24f Kfusion_sympy;
    for (int row=0; row<24; row++) {
        Hfusion_sympy[row] = Hfusion[row];
        Kfusion_sympy(row) = Kfusion(row);
    }

    // repeat calculation using matlab generated equations

	const float v_tas_pred = sqrtf((ve - vwe) * (ve - vwe) + (vn - vwn) * (vn - vwn) + vd * vd);

    // intermediate variable from algebraic optimisation
    float SH_TAS[3];
    SH_TAS[0] = 1.0f/v_tas_pred;
    SH_TAS[1] = (SH_TAS[0]*(2.0f*ve - 2.0f*vwe))*0.5f;
    SH_TAS[2] = (SH_TAS[0]*(2.0f*vn - 2.0f*vwn))*0.5f;

    // Observation Jacobian
    memset(&H_TAS, 0, sizeof(H_TAS));
    H_TAS(4) = SH_TAS[2];
    H_TAS(5) = SH_TAS[1];
    H_TAS(6) = vd*SH_TAS[0];
    H_TAS(22) = -SH_TAS[2];
    H_TAS(23) = -SH_TAS[1];

    airspeed_innov_var = (R_TAS + SH_TAS[2]*(P(4,4)*SH_TAS[2] + P(5,4)*SH_TAS[1] - P(22,4)*SH_TAS[2] - P(23,4)*SH_TAS[1] + P(6,4)*vd*SH_TAS[0]) + SH_TAS[1]*(P(4,5)*SH_TAS[2] + P(5,5)*SH_TAS[1] - P(22,5)*SH_TAS[2] - P(23,5)*SH_TAS[1] + P(6,5)*vd*SH_TAS[0]) - SH_TAS[2]*(P(4,22)*SH_TAS[2] + P(5,22)*SH_TAS[1] - P(22,22)*SH_TAS[2] - P(23,22)*SH_TAS[1] + P(6,22)*vd*SH_TAS[0]) - SH_TAS[1]*(P(4,23)*SH_TAS[2] + P(5,23)*SH_TAS[1] - P(22,23)*SH_TAS[2] - P(23,23)*SH_TAS[1] + P(6,23)*vd*SH_TAS[0]) + vd*SH_TAS[0]*(P(4,6)*SH_TAS[2] + P(5,6)*SH_TAS[1] - P(22,6)*SH_TAS[2] - P(23,6)*SH_TAS[1] + P(6,6)*vd*SH_TAS[0]));

    float SK_TAS[2];
    SK_TAS[0] = 1.0f / airspeed_innov_var;
    SK_TAS[1] = SH_TAS[1];

    // Kalman gain
    Kfusion(0) = SK_TAS[0]*(P(0,4)*SH_TAS[2] - P(0,22)*SH_TAS[2] + P(0,5)*SK_TAS[1] - P(0,23)*SK_TAS[1] + P(0,6)*vd*SH_TAS[0]);
    Kfusion(1) = SK_TAS[0]*(P(1,4)*SH_TAS[2] - P(1,22)*SH_TAS[2] + P(1,5)*SK_TAS[1] - P(1,23)*SK_TAS[1] + P(1,6)*vd*SH_TAS[0]);
    Kfusion(2) = SK_TAS[0]*(P(2,4)*SH_TAS[2] - P(2,22)*SH_TAS[2] + P(2,5)*SK_TAS[1] - P(2,23)*SK_TAS[1] + P(2,6)*vd*SH_TAS[0]);
    Kfusion(3) = SK_TAS[0]*(P(3,4)*SH_TAS[2] - P(3,22)*SH_TAS[2] + P(3,5)*SK_TAS[1] - P(3,23)*SK_TAS[1] + P(3,6)*vd*SH_TAS[0]);
    Kfusion(4) = SK_TAS[0]*(P(4,4)*SH_TAS[2] - P(4,22)*SH_TAS[2] + P(4,5)*SK_TAS[1] - P(4,23)*SK_TAS[1] + P(4,6)*vd*SH_TAS[0]);
    Kfusion(5) = SK_TAS[0]*(P(5,4)*SH_TAS[2] - P(5,22)*SH_TAS[2] + P(5,5)*SK_TAS[1] - P(5,23)*SK_TAS[1] + P(5,6)*vd*SH_TAS[0]);
    Kfusion(6) = SK_TAS[0]*(P(6,4)*SH_TAS[2] - P(6,22)*SH_TAS[2] + P(6,5)*SK_TAS[1] - P(6,23)*SK_TAS[1] + P(6,6)*vd*SH_TAS[0]);
    Kfusion(7) = SK_TAS[0]*(P(7,4)*SH_TAS[2] - P(7,22)*SH_TAS[2] + P(7,5)*SK_TAS[1] - P(7,23)*SK_TAS[1] + P(7,6)*vd*SH_TAS[0]);
    Kfusion(8) = SK_TAS[0]*(P(8,4)*SH_TAS[2] - P(8,22)*SH_TAS[2] + P(8,5)*SK_TAS[1] - P(8,23)*SK_TAS[1] + P(8,6)*vd*SH_TAS[0]);
    Kfusion(9) = SK_TAS[0]*(P(9,4)*SH_TAS[2] - P(9,22)*SH_TAS[2] + P(9,5)*SK_TAS[1] - P(9,23)*SK_TAS[1] + P(9,6)*vd*SH_TAS[0]);
    Kfusion(10) = SK_TAS[0]*(P(10,4)*SH_TAS[2] - P(10,22)*SH_TAS[2] + P(10,5)*SK_TAS[1] - P(10,23)*SK_TAS[1] + P(10,6)*vd*SH_TAS[0]);
    Kfusion(11) = SK_TAS[0]*(P(11,4)*SH_TAS[2] - P(11,22)*SH_TAS[2] + P(11,5)*SK_TAS[1] - P(11,23)*SK_TAS[1] + P(11,6)*vd*SH_TAS[0]);
    Kfusion(12) = SK_TAS[0]*(P(12,4)*SH_TAS[2] - P(12,22)*SH_TAS[2] + P(12,5)*SK_TAS[1] - P(12,23)*SK_TAS[1] + P(12,6)*vd*SH_TAS[0]);
    Kfusion(13) = SK_TAS[0]*(P(13,4)*SH_TAS[2] - P(13,22)*SH_TAS[2] + P(13,5)*SK_TAS[1] - P(13,23)*SK_TAS[1] + P(13,6)*vd*SH_TAS[0]);
    Kfusion(14) = SK_TAS[0]*(P(14,4)*SH_TAS[2] - P(14,22)*SH_TAS[2] + P(14,5)*SK_TAS[1] - P(14,23)*SK_TAS[1] + P(14,6)*vd*SH_TAS[0]);
    Kfusion(15) = SK_TAS[0]*(P(15,4)*SH_TAS[2] - P(15,22)*SH_TAS[2] + P(15,5)*SK_TAS[1] - P(15,23)*SK_TAS[1] + P(15,6)*vd*SH_TAS[0]);
    Kfusion(16) = SK_TAS[0]*(P(16,4)*SH_TAS[2] - P(16,22)*SH_TAS[2] + P(16,5)*SK_TAS[1] - P(16,23)*SK_TAS[1] + P(16,6)*vd*SH_TAS[0]);
    Kfusion(17) = SK_TAS[0]*(P(17,4)*SH_TAS[2] - P(17,22)*SH_TAS[2] + P(17,5)*SK_TAS[1] - P(17,23)*SK_TAS[1] + P(17,6)*vd*SH_TAS[0]);
    Kfusion(18) = SK_TAS[0]*(P(18,4)*SH_TAS[2] - P(18,22)*SH_TAS[2] + P(18,5)*SK_TAS[1] - P(18,23)*SK_TAS[1] + P(18,6)*vd*SH_TAS[0]);
    Kfusion(19) = SK_TAS[0]*(P(19,4)*SH_TAS[2] - P(19,22)*SH_TAS[2] + P(19,5)*SK_TAS[1] - P(19,23)*SK_TAS[1] + P(19,6)*vd*SH_TAS[0]);
    Kfusion(20) = SK_TAS[0]*(P(20,4)*SH_TAS[2] - P(20,22)*SH_TAS[2] + P(20,5)*SK_TAS[1] - P(20,23)*SK_TAS[1] + P(20,6)*vd*SH_TAS[0]);
    Kfusion(21) = SK_TAS[0]*(P(21,4)*SH_TAS[2] - P(21,22)*SH_TAS[2] + P(21,5)*SK_TAS[1] - P(21,23)*SK_TAS[1] + P(21,6)*vd*SH_TAS[0]);
    Kfusion(22) = SK_TAS[0]*(P(22,4)*SH_TAS[2] - P(22,22)*SH_TAS[2] + P(22,5)*SK_TAS[1] - P(22,23)*SK_TAS[1] + P(22,6)*vd*SH_TAS[0]);
    Kfusion(23) = SK_TAS[0]*(P(23,4)*SH_TAS[2] - P(23,22)*SH_TAS[2] + P(23,5)*SK_TAS[1] - P(23,23)*SK_TAS[1] + P(23,6)*vd*SH_TAS[0]);

    // find largest observation variance difference as a fraction of the matlab value
    float max_diff_fraction = 0.0f;
    int max_row;
    float max_old, max_new;
    for (int row=0; row<24; row++) {
        float diff_fraction;
        if (H_TAS(row) != 0.0f) {
            diff_fraction = fabsf(Hfusion_sympy[row] - H_TAS(row)) / fabsf(H_TAS(row));
        } else if (Hfusion_sympy[row] != 0.0f) {
            diff_fraction = fabsf(Hfusion_sympy[row] - H_TAS(row)) / fabsf(Hfusion_sympy[row]);
        } else {
            diff_fraction = 0.0f;
        }
        if (Hfusion_sympy[row] - H_TAS(row) != 0.0f) {
            printf("new,old Hfusion(%i) = %e,%e\n",row,Hfusion_sympy[row],H_TAS(row));
        }
        if (diff_fraction > max_diff_fraction) {
            max_diff_fraction = diff_fraction;
            max_row = row;
            max_old = H_TAS(row);
            max_new = Hfusion_sympy[row];
        }
    }

	if (max_diff_fraction > 1E-5f) {
		printf("Fail: Airspeed Hfusion max diff fraction = %e , old = %e , new = %e , location index = %i\n",max_diff_fraction, max_old, max_new, max_row);
	} else {
		printf("Pass: Airspeed Hfusion max diff fraction = %e\n",max_diff_fraction);
	}

    // find largest Kalman gain difference as a fraction of the matlab value
    max_diff_fraction = 0.0f;
    for (int row=0; row<4; row++) {
        float diff_fraction;
        if (Kfusion(row) != 0.0f) {
            diff_fraction = fabsf(Kfusion_sympy(row) - Kfusion(row)) / fabsf(Kfusion(row));
        } else if (Kfusion_sympy(row) != 0.0f) {
            diff_fraction = fabsf(Kfusion_sympy(row) - Kfusion(row)) / fabsf(Kfusion_sympy(row));
        } else {
            diff_fraction = 0.0f;
        }
        // if (Kfusion_sympy(row) - Kfusion(row) != 0.0f) {
        //     printf("new,old Kfusion(%i) = %e,%e\n",row,Kfusion_sympy(row),Kfusion(row));
        // }
        if (diff_fraction > max_diff_fraction) {
            max_diff_fraction = diff_fraction;
            max_row = row;
            max_old = Kfusion(row);
            max_new = Kfusion_sympy(row);
        }
    }

	if (max_diff_fraction > 1E-5f) {
		printf("Fail: Airspeed Kfusion max diff fraction = %e , old = %e , new = %e , location index = %i\n",max_diff_fraction, max_old, max_new, max_row);
	} else {
		printf("Pass: Airspeed Kfusion max diff fraction = %e\n",max_diff_fraction);
	}

    return 0;
}
