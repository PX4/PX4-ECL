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
    Vector24f H_MAG;
    Vector24f Kfusion;
    float mag_innov_var;

    // quaternion inputs must be normalised
    float q0 = 2.0f * ((float)rand() - 0.5f);
    float q1 = 2.0f * ((float)rand() - 0.5f);
    float q2 = 2.0f * ((float)rand() - 0.5f);
    float q3 = 2.0f * ((float)rand() - 0.5f);
    const float length = sqrtf(sq(q0) + sq(q1) + sq(q2) + sq(q3));
    q0 /= length;
    q1 /= length;
    q2 /= length;
    q3 /= length;

	const float magN = 2.0f * ((float)rand() - 0.5f);
	const float magE = 2.0f * ((float)rand() - 0.5f);
	const float magD = 2.0f * ((float)rand() - 0.5f);

	const float R_MAG = sq(0.05f);
    const bool update_all_states = true;

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

    // common expressions used by sympy generated equations
	const float HK0 = magE*q3;
	const float HK1 = magN*q0;
	const float HK2 = magD*q2;
	const float HK3 = HK0 + HK1 - HK2;
	const float HK4 = 2*HK3;
	const float HK5 = magD*q3 + magE*q2 + magN*q1;
	const float HK6 = 2*HK5;
	const float HK7 = magE*q1;
	const float HK8 = magD*q0;
	const float HK9 = magN*q2;
	const float HK10 = magD*q1;
	const float HK11 = magE*q0;
	const float HK12 = magN*q3;
	const float HK13 = HK10 + HK11 - HK12;
	const float HK14 = 2*HK13;
	const float HK15 = powf(q1, 2);
	const float HK16 = powf(q2, 2);
	const float HK17 = -HK16;
	const float HK18 = powf(q0, 2);
	const float HK19 = powf(q3, 2);
	const float HK20 = HK18 - HK19;
	const float HK21 = HK15 + HK17 + HK20;
	const float HK22 = q0*q3;
	const float HK23 = q1*q2;
	const float HK24 = HK22 + HK23;
	const float HK25 = q1*q3;
	const float HK26 = q0*q2;
	const float HK27 = 2*HK24;
	const float HK28 = -2*HK25 + 2*HK26;
	const float HK29 = 2*HK5;
	const float HK30 = 2*HK3;
	const float HK31 = -HK7 + HK8 + HK9;
	const float HK32 = 2*HK31;
	const float HK33 = HK32*P(0,2);
	const float HK34 = 2*HK13;
	const float HK35 = HK34*P(0,3);
	const float HK36 = HK21*P(0,16) + HK27*P(0,17) - HK28*P(0,18) + HK29*P(0,1) + HK30*P(0,0) - HK33 + HK35 + P(0,19);
	const float HK37 = HK21*P(16,16) + HK27*P(16,17) - HK28*P(16,18) + HK29*P(1,16) + HK30*P(0,16) - HK32*P(2,16) + HK34*P(3,16) + P(16,19);
	const float HK38 = HK21*P(16,18) + HK27*P(17,18) - HK28*P(18,18) + HK29*P(1,18) + HK30*P(0,18) - HK32*P(2,18) + HK34*P(3,18) + P(18,19);
	const float HK39 = HK29*P(1,2);
	const float HK40 = HK30*P(0,2);
	const float HK41 = HK21*P(2,16) + HK27*P(2,17) - HK28*P(2,18) - HK32*P(2,2) + HK34*P(2,3) + HK39 + HK40 + P(2,19);
	const float HK42 = HK21*P(16,17) + HK27*P(17,17) - HK28*P(17,18) + HK29*P(1,17) + HK30*P(0,17) - HK32*P(2,17) + HK34*P(3,17) + P(17,19);
	const float HK43 = HK29*P(1,3);
	const float HK44 = HK30*P(0,3);
	const float HK45 = HK21*P(3,16) + HK27*P(3,17) - HK28*P(3,18) - HK32*P(2,3) + HK34*P(3,3) + HK43 + HK44 + P(3,19);
	const float HK46 = HK32*P(1,2);
	const float HK47 = HK34*P(1,3);
	const float HK48 = HK21*P(1,16) + HK27*P(1,17) - HK28*P(1,18) + HK29*P(1,1) + HK30*P(0,1) - HK46 + HK47 + P(1,19);
	const float HK49 = HK21*P(16,19) + HK27*P(17,19) - HK28*P(18,19) + HK29*P(1,19) + HK30*P(0,19) - HK32*P(2,19) + HK34*P(3,19) + P(19,19);

    const float HK51 = 2*HK31;
    const float HK52 = -HK15;
    const float HK53 = HK16 + HK20 + HK52;
    const float HK54 = q0*q1;
    const float HK55 = q2*q3;
    const float HK56 = HK54 + HK55;
    const float HK57 = 2*HK56;
    const float HK58 = 2*HK22 - 2*HK23;
    const float HK59 = HK32*P(0,1);
    const float HK60 = HK29*P(0,2) + HK34*P(0,0) - HK44 + HK53*P(0,17) + HK57*P(0,18) - HK58*P(0,16) + HK59 + P(0,20);
    const float HK61 = HK29*P(2,17) - HK30*P(3,17) + HK32*P(1,17) + HK34*P(0,17) + HK53*P(17,17) + HK57*P(17,18) - HK58*P(16,17) + P(17,20);
    const float HK62 = HK29*P(2,16) - HK30*P(3,16) + HK32*P(1,16) + HK34*P(0,16) + HK53*P(16,17) + HK57*P(16,18) - HK58*P(16,16) + P(16,20);
    const float HK63 = HK29*P(2,3);
    const float HK64 = -HK30*P(3,3) + HK32*P(1,3) + HK35 + HK53*P(3,17) + HK57*P(3,18) - HK58*P(3,16) + HK63 + P(3,20);
    const float HK65 = HK29*P(2,18) - HK30*P(3,18) + HK32*P(1,18) + HK34*P(0,18) + HK53*P(17,18) + HK57*P(18,18) - HK58*P(16,18) + P(18,20);
    const float HK66 = HK34*P(0,1);
    const float HK67 = -HK30*P(1,3) + HK32*P(1,1) + HK39 + HK53*P(1,17) + HK57*P(1,18) - HK58*P(1,16) + HK66 + P(1,20);
    const float HK68 = HK30*P(2,3);
    const float HK69 = HK29*P(2,2) + HK34*P(0,2) + HK46 + HK53*P(2,17) + HK57*P(2,18) - HK58*P(2,16) - HK68 + P(2,20);
    const float HK70 = HK29*P(2,20) - HK30*P(3,20) + HK32*P(1,20) + HK34*P(0,20) + HK53*P(17,20) + HK57*P(18,20) - HK58*P(16,20) + P(20,20);

    const float HK72 = HK25 + HK26;
    const float HK73 = HK17 + HK18 + HK19 + HK52;
    const float HK74 = 2*HK72;
    const float HK75 = 2*HK54 - 2*HK55;
    const float HK76 = HK29*P(0,3) + HK32*P(0,0) + HK40 - HK66 + HK73*P(0,18) + HK74*P(0,16) - HK75*P(0,17) + P(0,21);
    const float HK77 = HK29*P(3,18) + HK30*P(2,18) + HK32*P(0,18) - HK34*P(1,18) + HK73*P(18,18) + HK74*P(16,18) - HK75*P(17,18) + P(18,21);
    const float HK78 = HK29*P(3,17) + HK30*P(2,17) + HK32*P(0,17) - HK34*P(1,17) + HK73*P(17,18) + HK74*P(16,17) - HK75*P(17,17) + P(17,21);
    const float HK79 = HK30*P(1,2) - HK34*P(1,1) + HK43 + HK59 + HK73*P(1,18) + HK74*P(1,16) - HK75*P(1,17) + P(1,21);
    const float HK80 = HK29*P(3,16) + HK30*P(2,16) + HK32*P(0,16) - HK34*P(1,16) + HK73*P(16,18) + HK74*P(16,16) - HK75*P(16,17) + P(16,21);
    const float HK81 = HK29*P(3,3) + HK32*P(0,3) - HK47 + HK68 + HK73*P(3,18) + HK74*P(3,16) - HK75*P(3,17) + P(3,21);
    const float HK82 = HK30*P(2,2) + HK33 - HK34*P(1,2) + HK63 + HK73*P(2,18) + HK74*P(2,16) - HK75*P(2,17) + P(2,21);
    const float HK83 = HK29*P(3,21) + HK30*P(2,21) + HK32*P(0,21) - HK34*P(1,21) + HK73*P(18,21) + HK74*P(16,21) - HK75*P(17,21) + P(21,21);

    // common expressions used by matlab generated equations
	float SH_MAG[9];
	SH_MAG[0] = 2.0f*magD*q3 + 2.0f*magE*q2 + 2.0f*magN*q1;
	SH_MAG[1] = 2.0f*magD*q0 - 2.0f*magE*q1 + 2.0f*magN*q2;
	SH_MAG[2] = 2.0f*magD*q1 + 2.0f*magE*q0 - 2.0f*magN*q3;
	SH_MAG[3] = sq(q3);
	SH_MAG[4] = sq(q2);
	SH_MAG[5] = sq(q1);
	SH_MAG[6] = sq(q0);
	SH_MAG[7] = 2.0f*magN*q0;
	SH_MAG[8] = 2.0f*magE*q3;

    // Compare X axis equations
    {
        mag_innov_var = (HK21*HK37 + HK27*HK42 - HK28*HK38 + HK29*HK48 + HK30*HK36 - HK32*HK41 + HK34*HK45 + HK49 + R_MAG);
        float HK50 = 1.0F/mag_innov_var;

        memset(Hfusion, 0, sizeof(Hfusion));
        Hfusion[0] = HK4;
        Hfusion[1] = HK6;
        Hfusion[2] = 2*HK7 - 2*HK8 - 2*HK9;
        Hfusion[3] = HK14;
        Hfusion[16] = HK21;
        Hfusion[17] = 2*HK24;
        Hfusion[18] = 2*HK25 - 2*HK26;
        Hfusion[19] = 1;

        if (update_all_states) {
            Kfusion(0) = HK36*HK50;
            Kfusion(1) = HK48*HK50;
            Kfusion(2) = HK41*HK50;
            Kfusion(3) = HK45*HK50;

            for (unsigned row = 4; row <= 15; row++) {
                Kfusion(row) = HK50*(HK21*P(row,16) + HK27*P(row,17) - HK28*P(row,18) + HK29*P(1,row) + HK30*P(0,row) - HK32*P(2,row) + HK34*P(3,row) + P(row,19));
            }

            for (unsigned row = 22; row <= 23; row++) {
                Kfusion(row) = HK50*(HK21*P(row,16) + HK27*P(row,17) - HK28*P(row,18) + HK29*P(1,row) + HK30*P(0,row) - HK32*P(2,row) + HK34*P(3,row) + P(row,19));
            }
        } else {
            for (uint8_t i = 0; i < 16; i++) {
                Kfusion(i) = 0.0f;
            }

            Kfusion(22) = 0.0f;
            Kfusion(23) = 0.0f;
        }

        Kfusion(16) = HK37*HK50;
        Kfusion(17) = HK42*HK50;
        Kfusion(18) = HK38*HK50;
        Kfusion(19) = HK49*HK50;

        for (unsigned row = 20; row <= 21; row++) {
            Kfusion(row) = HK50*(HK21*P(16,row) + HK27*P(17,row) - HK28*P(18,row) + HK29*P(1,row) + HK30*P(0,row) - HK32*P(2,row) + HK34*P(3,row) + P(19,row));
        }

        // save output and repeat calculation using legacy matlab generated code
        float Hfusion_sympy[24];
        Vector24f Kfusion_sympy;
        for (int row=0; row<24; row++) {
            Hfusion_sympy[row] = Hfusion[row];
            Kfusion_sympy(row) = Kfusion(row);
        }

        // repeat calculation using matlab generated equations
     	// X axis innovation variance
	    mag_innov_var = (P(19,19) + R_MAG + P(1,19)*SH_MAG[0] - P(2,19)*SH_MAG[1] + P(3,19)*SH_MAG[2] - P(16,19)*(SH_MAG[3] + SH_MAG[4] - SH_MAG[5] - SH_MAG[6]) + (2.0f*q0*q3 + 2.0f*q1*q2)*(P(19,17) + P(1,17)*SH_MAG[0] - P(2,17)*SH_MAG[1] + P(3,17)*SH_MAG[2] - P(16,17)*(SH_MAG[3] + SH_MAG[4] - SH_MAG[5] - SH_MAG[6]) + P(17,17)*(2.0f*q0*q3 + 2.0f*q1*q2) - P(18,17)*(2.0f*q0*q2 - 2.0f*q1*q3) + P(0,17)*(SH_MAG[7] + SH_MAG[8] - 2.0f*magD*q2)) - (2.0f*q0*q2 - 2.0f*q1*q3)*(P(19,18) + P(1,18)*SH_MAG[0] - P(2,18)*SH_MAG[1] + P(3,18)*SH_MAG[2] - P(16,18)*(SH_MAG[3] + SH_MAG[4] - SH_MAG[5] - SH_MAG[6]) + P(17,18)*(2.0f*q0*q3 + 2.0f*q1*q2) - P(18,18)*(2.0f*q0*q2 - 2.0f*q1*q3) + P(0,18)*(SH_MAG[7] + SH_MAG[8] - 2.0f*magD*q2)) + (SH_MAG[7] + SH_MAG[8] - 2.0f*magD*q2)*(P(19,0) + P(1,0)*SH_MAG[0] - P(2,0)*SH_MAG[1] + P(3,0)*SH_MAG[2] - P(16,0)*(SH_MAG[3] + SH_MAG[4] - SH_MAG[5] - SH_MAG[6]) + P(17,0)*(2.0f*q0*q3 + 2.0f*q1*q2) - P(18,0)*(2.0f*q0*q2 - 2.0f*q1*q3) + P(0,0)*(SH_MAG[7] + SH_MAG[8] - 2.0f*magD*q2)) + P(17,19)*(2.0f*q0*q3 + 2.0f*q1*q2) - P(18,19)*(2.0f*q0*q2 - 2.0f*q1*q3) + SH_MAG[0]*(P(19,1) + P(1,1)*SH_MAG[0] - P(2,1)*SH_MAG[1] + P(3,1)*SH_MAG[2] - P(16,1)*(SH_MAG[3] + SH_MAG[4] - SH_MAG[5] - SH_MAG[6]) + P(17,1)*(2.0f*q0*q3 + 2.0f*q1*q2) - P(18,1)*(2.0f*q0*q2 - 2.0f*q1*q3) + P(0,1)*(SH_MAG[7] + SH_MAG[8] - 2.0f*magD*q2)) - SH_MAG[1]*(P(19,2) + P(1,2)*SH_MAG[0] - P(2,2)*SH_MAG[1] + P(3,2)*SH_MAG[2] - P(16,2)*(SH_MAG[3] + SH_MAG[4] - SH_MAG[5] - SH_MAG[6]) + P(17,2)*(2.0f*q0*q3 + 2.0f*q1*q2) - P(18,2)*(2.0f*q0*q2 - 2.0f*q1*q3) + P(0,2)*(SH_MAG[7] + SH_MAG[8] - 2.0f*magD*q2)) + SH_MAG[2]*(P(19,3) + P(1,3)*SH_MAG[0] - P(2,3)*SH_MAG[1] + P(3,3)*SH_MAG[2] - P(16,3)*(SH_MAG[3] + SH_MAG[4] - SH_MAG[5] - SH_MAG[6]) + P(17,3)*(2.0f*q0*q3 + 2.0f*q1*q2) - P(18,3)*(2.0f*q0*q2 - 2.0f*q1*q3) + P(0,3)*(SH_MAG[7] + SH_MAG[8] - 2.0f*magD*q2)) - (SH_MAG[3] + SH_MAG[4] - SH_MAG[5] - SH_MAG[6])*(P(19,16) + P(1,16)*SH_MAG[0] - P(2,16)*SH_MAG[1] + P(3,16)*SH_MAG[2] - P(16,16)*(SH_MAG[3] + SH_MAG[4] - SH_MAG[5] - SH_MAG[6]) + P(17,16)*(2.0f*q0*q3 + 2.0f*q1*q2) - P(18,16)*(2.0f*q0*q2 - 2.0f*q1*q3) + P(0,16)*(SH_MAG[7] + SH_MAG[8] - 2.0f*magD*q2)) + P(0,19)*(SH_MAG[7] + SH_MAG[8] - 2.0f*magD*q2));

			// Calculate X axis observation jacobians
			H_MAG.setZero();
			H_MAG(0) = SH_MAG[7] + SH_MAG[8] - 2.0f*magD*q2;
			H_MAG(1) = SH_MAG[0];
			H_MAG(2) = -SH_MAG[1];
			H_MAG(3) = SH_MAG[2];
			H_MAG(16) = SH_MAG[5] - SH_MAG[4] - SH_MAG[3] + SH_MAG[6];
			H_MAG(17) = 2.0f*q0*q3 + 2.0f*q1*q2;
			H_MAG(18) = 2.0f*q1*q3 - 2.0f*q0*q2;
			H_MAG(19) = 1.0f;

			// Calculate X axis Kalman gains
			float SK_MX[5];
			SK_MX[0] = 1.0f / mag_innov_var;
			SK_MX[1] = SH_MAG[3] + SH_MAG[4] - SH_MAG[5] - SH_MAG[6];
			SK_MX[2] = SH_MAG[7] + SH_MAG[8] - 2.0f*magD*q2;
			SK_MX[3] = 2.0f*q0*q2 - 2.0f*q1*q3;
			SK_MX[4] = 2.0f*q0*q3 + 2.0f*q1*q2;

			if (update_all_states) {
				Kfusion(0) = SK_MX[0]*(P(0,19) + P(0,1)*SH_MAG[0] - P(0,2)*SH_MAG[1] + P(0,3)*SH_MAG[2] + P(0,0)*SK_MX[2] - P(0,16)*SK_MX[1] + P(0,17)*SK_MX[4] - P(0,18)*SK_MX[3]);
				Kfusion(1) = SK_MX[0]*(P(1,19) + P(1,1)*SH_MAG[0] - P(1,2)*SH_MAG[1] + P(1,3)*SH_MAG[2] + P(1,0)*SK_MX[2] - P(1,16)*SK_MX[1] + P(1,17)*SK_MX[4] - P(1,18)*SK_MX[3]);
				Kfusion(2) = SK_MX[0]*(P(2,19) + P(2,1)*SH_MAG[0] - P(2,2)*SH_MAG[1] + P(2,3)*SH_MAG[2] + P(2,0)*SK_MX[2] - P(2,16)*SK_MX[1] + P(2,17)*SK_MX[4] - P(2,18)*SK_MX[3]);
				Kfusion(3) = SK_MX[0]*(P(3,19) + P(3,1)*SH_MAG[0] - P(3,2)*SH_MAG[1] + P(3,3)*SH_MAG[2] + P(3,0)*SK_MX[2] - P(3,16)*SK_MX[1] + P(3,17)*SK_MX[4] - P(3,18)*SK_MX[3]);
				Kfusion(4) = SK_MX[0]*(P(4,19) + P(4,1)*SH_MAG[0] - P(4,2)*SH_MAG[1] + P(4,3)*SH_MAG[2] + P(4,0)*SK_MX[2] - P(4,16)*SK_MX[1] + P(4,17)*SK_MX[4] - P(4,18)*SK_MX[3]);
				Kfusion(5) = SK_MX[0]*(P(5,19) + P(5,1)*SH_MAG[0] - P(5,2)*SH_MAG[1] + P(5,3)*SH_MAG[2] + P(5,0)*SK_MX[2] - P(5,16)*SK_MX[1] + P(5,17)*SK_MX[4] - P(5,18)*SK_MX[3]);
				Kfusion(6) = SK_MX[0]*(P(6,19) + P(6,1)*SH_MAG[0] - P(6,2)*SH_MAG[1] + P(6,3)*SH_MAG[2] + P(6,0)*SK_MX[2] - P(6,16)*SK_MX[1] + P(6,17)*SK_MX[4] - P(6,18)*SK_MX[3]);
				Kfusion(7) = SK_MX[0]*(P(7,19) + P(7,1)*SH_MAG[0] - P(7,2)*SH_MAG[1] + P(7,3)*SH_MAG[2] + P(7,0)*SK_MX[2] - P(7,16)*SK_MX[1] + P(7,17)*SK_MX[4] - P(7,18)*SK_MX[3]);
				Kfusion(8) = SK_MX[0]*(P(8,19) + P(8,1)*SH_MAG[0] - P(8,2)*SH_MAG[1] + P(8,3)*SH_MAG[2] + P(8,0)*SK_MX[2] - P(8,16)*SK_MX[1] + P(8,17)*SK_MX[4] - P(8,18)*SK_MX[3]);
				Kfusion(9) = SK_MX[0]*(P(9,19) + P(9,1)*SH_MAG[0] - P(9,2)*SH_MAG[1] + P(9,3)*SH_MAG[2] + P(9,0)*SK_MX[2] - P(9,16)*SK_MX[1] + P(9,17)*SK_MX[4] - P(9,18)*SK_MX[3]);
				Kfusion(10) = SK_MX[0]*(P(10,19) + P(10,1)*SH_MAG[0] - P(10,2)*SH_MAG[1] + P(10,3)*SH_MAG[2] + P(10,0)*SK_MX[2] - P(10,16)*SK_MX[1] + P(10,17)*SK_MX[4] - P(10,18)*SK_MX[3]);
				Kfusion(11) = SK_MX[0]*(P(11,19) + P(11,1)*SH_MAG[0] - P(11,2)*SH_MAG[1] + P(11,3)*SH_MAG[2] + P(11,0)*SK_MX[2] - P(11,16)*SK_MX[1] + P(11,17)*SK_MX[4] - P(11,18)*SK_MX[3]);
				Kfusion(12) = SK_MX[0]*(P(12,19) + P(12,1)*SH_MAG[0] - P(12,2)*SH_MAG[1] + P(12,3)*SH_MAG[2] + P(12,0)*SK_MX[2] - P(12,16)*SK_MX[1] + P(12,17)*SK_MX[4] - P(12,18)*SK_MX[3]);
				Kfusion(13) = SK_MX[0]*(P(13,19) + P(13,1)*SH_MAG[0] - P(13,2)*SH_MAG[1] + P(13,3)*SH_MAG[2] + P(13,0)*SK_MX[2] - P(13,16)*SK_MX[1] + P(13,17)*SK_MX[4] - P(13,18)*SK_MX[3]);
				Kfusion(14) = SK_MX[0]*(P(14,19) + P(14,1)*SH_MAG[0] - P(14,2)*SH_MAG[1] + P(14,3)*SH_MAG[2] + P(14,0)*SK_MX[2] - P(14,16)*SK_MX[1] + P(14,17)*SK_MX[4] - P(14,18)*SK_MX[3]);
				Kfusion(15) = SK_MX[0]*(P(15,19) + P(15,1)*SH_MAG[0] - P(15,2)*SH_MAG[1] + P(15,3)*SH_MAG[2] + P(15,0)*SK_MX[2] - P(15,16)*SK_MX[1] + P(15,17)*SK_MX[4] - P(15,18)*SK_MX[3]);
				Kfusion(22) = SK_MX[0]*(P(22,19) + P(22,1)*SH_MAG[0] - P(22,2)*SH_MAG[1] + P(22,3)*SH_MAG[2] + P(22,0)*SK_MX[2] - P(22,16)*SK_MX[1] + P(22,17)*SK_MX[4] - P(22,18)*SK_MX[3]);
				Kfusion(23) = SK_MX[0]*(P(23,19) + P(23,1)*SH_MAG[0] - P(23,2)*SH_MAG[1] + P(23,3)*SH_MAG[2] + P(23,0)*SK_MX[2] - P(23,16)*SK_MX[1] + P(23,17)*SK_MX[4] - P(23,18)*SK_MX[3]);
			}

			Kfusion(16) = SK_MX[0]*(P(16,19) + P(16,1)*SH_MAG[0] - P(16,2)*SH_MAG[1] + P(16,3)*SH_MAG[2] + P(16,0)*SK_MX[2] - P(16,16)*SK_MX[1] + P(16,17)*SK_MX[4] - P(16,18)*SK_MX[3]);
			Kfusion(17) = SK_MX[0]*(P(17,19) + P(17,1)*SH_MAG[0] - P(17,2)*SH_MAG[1] + P(17,3)*SH_MAG[2] + P(17,0)*SK_MX[2] - P(17,16)*SK_MX[1] + P(17,17)*SK_MX[4] - P(17,18)*SK_MX[3]);
			Kfusion(18) = SK_MX[0]*(P(18,19) + P(18,1)*SH_MAG[0] - P(18,2)*SH_MAG[1] + P(18,3)*SH_MAG[2] + P(18,0)*SK_MX[2] - P(18,16)*SK_MX[1] + P(18,17)*SK_MX[4] - P(18,18)*SK_MX[3]);
			Kfusion(19) = SK_MX[0]*(P(19,19) + P(19,1)*SH_MAG[0] - P(19,2)*SH_MAG[1] + P(19,3)*SH_MAG[2] + P(19,0)*SK_MX[2] - P(19,16)*SK_MX[1] + P(19,17)*SK_MX[4] - P(19,18)*SK_MX[3]);
			Kfusion(20) = SK_MX[0]*(P(20,19) + P(20,1)*SH_MAG[0] - P(20,2)*SH_MAG[1] + P(20,3)*SH_MAG[2] + P(20,0)*SK_MX[2] - P(20,16)*SK_MX[1] + P(20,17)*SK_MX[4] - P(20,18)*SK_MX[3]);
			Kfusion(21) = SK_MX[0]*(P(21,19) + P(21,1)*SH_MAG[0] - P(21,2)*SH_MAG[1] + P(21,3)*SH_MAG[2] + P(21,0)*SK_MX[2] - P(21,16)*SK_MX[1] + P(21,17)*SK_MX[4] - P(21,18)*SK_MX[3]);

        // find largest observation variance difference as a fraction of the matlab value
        float max_diff_fraction = 0.0f;
        int max_row;
        float max_old, max_new;
        for (int row=0; row<24; row++) {
            float diff_fraction;
            if (H_MAG(row) != 0.0f) {
                diff_fraction = fabsf(Hfusion_sympy[row] - H_MAG(row)) / fabsf(H_MAG(row));
            } else if (Hfusion_sympy[row] != 0.0f) {
                diff_fraction = fabsf(Hfusion_sympy[row] - H_MAG(row)) / fabsf(Hfusion_sympy[row]);
            } else {
                diff_fraction = 0.0f;
            }
            if (diff_fraction > max_diff_fraction) {
                max_diff_fraction = diff_fraction;
                max_row = row;
                max_old = H_MAG(row);
                max_new = Hfusion_sympy[row];
            }
        }

        printf("X axis Hfusion max_diff_fraction = %e , old = %e , new = %e , location index = %i\n",max_diff_fraction, max_old, max_new, max_row);

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
            if (diff_fraction > max_diff_fraction) {
                max_diff_fraction = diff_fraction;
                max_row = row;
                max_old = Kfusion(row);
                max_new = Kfusion_sympy(row);
            }
        }

        printf("X axis Kfusion max_diff_fraction = %e , old = %e , new = %e , location index = %i\n",max_diff_fraction, max_old, max_new, max_row);

    }

    // Compare Y axis equations
    {
        mag_innov_var = (HK29*HK69 - HK30*HK64 + HK32*HK67 + HK34*HK60 + HK53*HK61 + HK57*HK65 - HK58*HK62 + HK70 + R_MAG);
        float HK71 = 1.0F/mag_innov_var;

        // Calculate Y axis observation jacobians
        memset(Hfusion, 0, sizeof(Hfusion));
        Hfusion[0] = HK14;
        Hfusion[1] = HK51;
        Hfusion[2] = HK6;
        Hfusion[3] = -2*HK0 - 2*HK1 + 2*HK2;
        Hfusion[16] = -2*HK22 + 2*HK23;
        Hfusion[17] = HK53;
        Hfusion[18] = 2*HK56;
        Hfusion[20] = 1;

        // Calculate Y axis Kalman gains
        if (update_all_states) {
            Kfusion(0) = HK60*HK71;
            Kfusion(1) = HK67*HK71;
            Kfusion(2) = HK69*HK71;
            Kfusion(3) = HK64*HK71;

            for (unsigned row = 4; row <= 15; row++) {
                Kfusion(row) = HK71*(HK29*P(2,row) - HK30*P(3,row) + HK32*P(1,row) + HK34*P(0,row) + HK53*P(row,17) + HK57*P(row,18) - HK58*P(row,16) + P(row,20));
            }

            for (unsigned row = 22; row <= 23; row++) {
                Kfusion(row) = HK71*(HK29*P(2,row) - HK30*P(3,row) + HK32*P(1,row) + HK34*P(0,row) + HK53*P(row,17) + HK57*P(row,18) - HK58*P(row,16) + P(row,20));
            }
        } else {
            for (uint8_t i = 0; i < 16; i++) {
                Kfusion(i) = 0.0f;
            }

            Kfusion(22) = 0.0f;
            Kfusion(23) = 0.0f;
        }

        Kfusion(16) = HK62*HK71;
        Kfusion(17) = HK61*HK71;
        Kfusion(18) = HK65*HK71;
        Kfusion(19) = HK71*(HK29*P(2,19) - HK30*P(3,19) + HK32*P(1,19) + HK34*P(0,19) + HK53*P(17,19) + HK57*P(18,19) - HK58*P(16,19) + P(19,20));
        Kfusion(20) = HK70*HK71;
        Kfusion(21) = HK71*(HK29*P(2,21) - HK30*P(3,21) + HK32*P(1,21) + HK34*P(0,21) + HK53*P(17,21) + HK57*P(18,21) - HK58*P(16,21) + P(20,21));

        // save output and repeat calculation using legacy matlab generated code
        float Hfusion_sympy[24];
        Vector24f Kfusion_sympy;
        for (int row=0; row<24; row++) {
            Hfusion_sympy[row] = Hfusion[row];
            Kfusion_sympy(row) = Kfusion(row);
        }

        // repeat calculation using matlab generated equations

        // Y axis innovation variance
        mag_innov_var = (P(20,20) + R_MAG + P(0,20)*SH_MAG[2] + P(1,20)*SH_MAG[1] + P(2,20)*SH_MAG[0] - P(17,20)*(SH_MAG[3] - SH_MAG[4] + SH_MAG[5] - SH_MAG[6]) - (2.0f*q0*q3 - 2.0f*q1*q2)*(P(20,16) + P(0,16)*SH_MAG[2] + P(1,16)*SH_MAG[1] + P(2,16)*SH_MAG[0] - P(17,16)*(SH_MAG[3] - SH_MAG[4] + SH_MAG[5] - SH_MAG[6]) - P(16,16)*(2.0f*q0*q3 - 2.0f*q1*q2) + P(18,16)*(2.0f*q0*q1 + 2.0f*q2*q3) - P(3,16)*(SH_MAG[7] + SH_MAG[8] - 2.0f*magD*q2)) + (2.0f*q0*q1 + 2.0f*q2*q3)*(P(20,18) + P(0,18)*SH_MAG[2] + P(1,18)*SH_MAG[1] + P(2,18)*SH_MAG[0] - P(17,18)*(SH_MAG[3] - SH_MAG[4] + SH_MAG[5] - SH_MAG[6]) - P(16,18)*(2.0f*q0*q3 - 2.0f*q1*q2) + P(18,18)*(2.0f*q0*q1 + 2.0f*q2*q3) - P(3,18)*(SH_MAG[7] + SH_MAG[8] - 2.0f*magD*q2)) - (SH_MAG[7] + SH_MAG[8] - 2.0f*magD*q2)*(P(20,3) + P(0,3)*SH_MAG[2] + P(1,3)*SH_MAG[1] + P(2,3)*SH_MAG[0] - P(17,3)*(SH_MAG[3] - SH_MAG[4] + SH_MAG[5] - SH_MAG[6]) - P(16,3)*(2.0f*q0*q3 - 2.0f*q1*q2) + P(18,3)*(2.0f*q0*q1 + 2.0f*q2*q3) - P(3,3)*(SH_MAG[7] + SH_MAG[8] - 2.0f*magD*q2)) - P(16,20)*(2.0f*q0*q3 - 2.0f*q1*q2) + P(18,20)*(2.0f*q0*q1 + 2.0f*q2*q3) + SH_MAG[2]*(P(20,0) + P(0,0)*SH_MAG[2] + P(1,0)*SH_MAG[1] + P(2,0)*SH_MAG[0] - P(17,0)*(SH_MAG[3] - SH_MAG[4] + SH_MAG[5] - SH_MAG[6]) - P(16,0)*(2.0f*q0*q3 - 2.0f*q1*q2) + P(18,0)*(2.0f*q0*q1 + 2.0f*q2*q3) - P(3,0)*(SH_MAG[7] + SH_MAG[8] - 2.0f*magD*q2)) + SH_MAG[1]*(P(20,1) + P(0,1)*SH_MAG[2] + P(1,1)*SH_MAG[1] + P(2,1)*SH_MAG[0] - P(17,1)*(SH_MAG[3] - SH_MAG[4] + SH_MAG[5] - SH_MAG[6]) - P(16,1)*(2.0f*q0*q3 - 2.0f*q1*q2) + P(18,1)*(2.0f*q0*q1 + 2.0f*q2*q3) - P(3,1)*(SH_MAG[7] + SH_MAG[8] - 2.0f*magD*q2)) + SH_MAG[0]*(P(20,2) + P(0,2)*SH_MAG[2] + P(1,2)*SH_MAG[1] + P(2,2)*SH_MAG[0] - P(17,2)*(SH_MAG[3] - SH_MAG[4] + SH_MAG[5] - SH_MAG[6]) - P(16,2)*(2.0f*q0*q3 - 2.0f*q1*q2) + P(18,2)*(2.0f*q0*q1 + 2.0f*q2*q3) - P(3,2)*(SH_MAG[7] + SH_MAG[8] - 2.0f*magD*q2)) - (SH_MAG[3] - SH_MAG[4] + SH_MAG[5] - SH_MAG[6])*(P(20,17) + P(0,17)*SH_MAG[2] + P(1,17)*SH_MAG[1] + P(2,17)*SH_MAG[0] - P(17,17)*(SH_MAG[3] - SH_MAG[4] + SH_MAG[5] - SH_MAG[6]) - P(16,17)*(2.0f*q0*q3 - 2.0f*q1*q2) + P(18,17)*(2.0f*q0*q1 + 2.0f*q2*q3) - P(3,17)*(SH_MAG[7] + SH_MAG[8] - 2.0f*magD*q2)) - P(3,20)*(SH_MAG[7] + SH_MAG[8] - 2.0f*magD*q2));

        // Calculate Y axis observation jacobians
        H_MAG.setZero();
        H_MAG(0) = SH_MAG[2];
        H_MAG(1) = SH_MAG[1];
        H_MAG(2) = SH_MAG[0];
        H_MAG(3) = 2.0f*magD*q2 - SH_MAG[8] - SH_MAG[7];
        H_MAG(16) = 2.0f*q1*q2 - 2.0f*q0*q3;
        H_MAG(17) = SH_MAG[4] - SH_MAG[3] - SH_MAG[5] + SH_MAG[6];
        H_MAG(18) = 2.0f*q0*q1 + 2.0f*q2*q3;
        H_MAG(20) = 1.0f;

        // Calculate Y axis Kalman gains
        float SK_MY[5];
        SK_MY[0] = 1.0f / mag_innov_var;
        SK_MY[1] = SH_MAG[3] - SH_MAG[4] + SH_MAG[5] - SH_MAG[6];
        SK_MY[2] = SH_MAG[7] + SH_MAG[8] - 2.0f*magD*q2;
        SK_MY[3] = 2.0f*q0*q3 - 2.0f*q1*q2;
        SK_MY[4] = 2.0f*q0*q1 + 2.0f*q2*q3;

        if (update_all_states) {
            Kfusion(0) = SK_MY[0]*(P(0,20) + P(0,0)*SH_MAG[2] + P(0,1)*SH_MAG[1] + P(0,2)*SH_MAG[0] - P(0,3)*SK_MY[2] - P(0,17)*SK_MY[1] - P(0,16)*SK_MY[3] + P(0,18)*SK_MY[4]);
            Kfusion(1) = SK_MY[0]*(P(1,20) + P(1,0)*SH_MAG[2] + P(1,1)*SH_MAG[1] + P(1,2)*SH_MAG[0] - P(1,3)*SK_MY[2] - P(1,17)*SK_MY[1] - P(1,16)*SK_MY[3] + P(1,18)*SK_MY[4]);
            Kfusion(2) = SK_MY[0]*(P(2,20) + P(2,0)*SH_MAG[2] + P(2,1)*SH_MAG[1] + P(2,2)*SH_MAG[0] - P(2,3)*SK_MY[2] - P(2,17)*SK_MY[1] - P(2,16)*SK_MY[3] + P(2,18)*SK_MY[4]);
            Kfusion(3) = SK_MY[0]*(P(3,20) + P(3,0)*SH_MAG[2] + P(3,1)*SH_MAG[1] + P(3,2)*SH_MAG[0] - P(3,3)*SK_MY[2] - P(3,17)*SK_MY[1] - P(3,16)*SK_MY[3] + P(3,18)*SK_MY[4]);
            Kfusion(4) = SK_MY[0]*(P(4,20) + P(4,0)*SH_MAG[2] + P(4,1)*SH_MAG[1] + P(4,2)*SH_MAG[0] - P(4,3)*SK_MY[2] - P(4,17)*SK_MY[1] - P(4,16)*SK_MY[3] + P(4,18)*SK_MY[4]);
            Kfusion(5) = SK_MY[0]*(P(5,20) + P(5,0)*SH_MAG[2] + P(5,1)*SH_MAG[1] + P(5,2)*SH_MAG[0] - P(5,3)*SK_MY[2] - P(5,17)*SK_MY[1] - P(5,16)*SK_MY[3] + P(5,18)*SK_MY[4]);
            Kfusion(6) = SK_MY[0]*(P(6,20) + P(6,0)*SH_MAG[2] + P(6,1)*SH_MAG[1] + P(6,2)*SH_MAG[0] - P(6,3)*SK_MY[2] - P(6,17)*SK_MY[1] - P(6,16)*SK_MY[3] + P(6,18)*SK_MY[4]);
            Kfusion(7) = SK_MY[0]*(P(7,20) + P(7,0)*SH_MAG[2] + P(7,1)*SH_MAG[1] + P(7,2)*SH_MAG[0] - P(7,3)*SK_MY[2] - P(7,17)*SK_MY[1] - P(7,16)*SK_MY[3] + P(7,18)*SK_MY[4]);
            Kfusion(8) = SK_MY[0]*(P(8,20) + P(8,0)*SH_MAG[2] + P(8,1)*SH_MAG[1] + P(8,2)*SH_MAG[0] - P(8,3)*SK_MY[2] - P(8,17)*SK_MY[1] - P(8,16)*SK_MY[3] + P(8,18)*SK_MY[4]);
            Kfusion(9) = SK_MY[0]*(P(9,20) + P(9,0)*SH_MAG[2] + P(9,1)*SH_MAG[1] + P(9,2)*SH_MAG[0] - P(9,3)*SK_MY[2] - P(9,17)*SK_MY[1] - P(9,16)*SK_MY[3] + P(9,18)*SK_MY[4]);
            Kfusion(10) = SK_MY[0]*(P(10,20) + P(10,0)*SH_MAG[2] + P(10,1)*SH_MAG[1] + P(10,2)*SH_MAG[0] - P(10,3)*SK_MY[2] - P(10,17)*SK_MY[1] - P(10,16)*SK_MY[3] + P(10,18)*SK_MY[4]);
            Kfusion(11) = SK_MY[0]*(P(11,20) + P(11,0)*SH_MAG[2] + P(11,1)*SH_MAG[1] + P(11,2)*SH_MAG[0] - P(11,3)*SK_MY[2] - P(11,17)*SK_MY[1] - P(11,16)*SK_MY[3] + P(11,18)*SK_MY[4]);
            Kfusion(12) = SK_MY[0]*(P(12,20) + P(12,0)*SH_MAG[2] + P(12,1)*SH_MAG[1] + P(12,2)*SH_MAG[0] - P(12,3)*SK_MY[2] - P(12,17)*SK_MY[1] - P(12,16)*SK_MY[3] + P(12,18)*SK_MY[4]);
            Kfusion(13) = SK_MY[0]*(P(13,20) + P(13,0)*SH_MAG[2] + P(13,1)*SH_MAG[1] + P(13,2)*SH_MAG[0] - P(13,3)*SK_MY[2] - P(13,17)*SK_MY[1] - P(13,16)*SK_MY[3] + P(13,18)*SK_MY[4]);
            Kfusion(14) = SK_MY[0]*(P(14,20) + P(14,0)*SH_MAG[2] + P(14,1)*SH_MAG[1] + P(14,2)*SH_MAG[0] - P(14,3)*SK_MY[2] - P(14,17)*SK_MY[1] - P(14,16)*SK_MY[3] + P(14,18)*SK_MY[4]);
            Kfusion(15) = SK_MY[0]*(P(15,20) + P(15,0)*SH_MAG[2] + P(15,1)*SH_MAG[1] + P(15,2)*SH_MAG[0] - P(15,3)*SK_MY[2] - P(15,17)*SK_MY[1] - P(15,16)*SK_MY[3] + P(15,18)*SK_MY[4]);
            Kfusion(22) = SK_MY[0]*(P(22,20) + P(22,0)*SH_MAG[2] + P(22,1)*SH_MAG[1] + P(22,2)*SH_MAG[0] - P(22,3)*SK_MY[2] - P(22,17)*SK_MY[1] - P(22,16)*SK_MY[3] + P(22,18)*SK_MY[4]);
            Kfusion(23) = SK_MY[0]*(P(23,20) + P(23,0)*SH_MAG[2] + P(23,1)*SH_MAG[1] + P(23,2)*SH_MAG[0] - P(23,3)*SK_MY[2] - P(23,17)*SK_MY[1] - P(23,16)*SK_MY[3] + P(23,18)*SK_MY[4]);
        }

        Kfusion(16) = SK_MY[0]*(P(16,20) + P(16,0)*SH_MAG[2] + P(16,1)*SH_MAG[1] + P(16,2)*SH_MAG[0] - P(16,3)*SK_MY[2] - P(16,17)*SK_MY[1] - P(16,16)*SK_MY[3] + P(16,18)*SK_MY[4]);
        Kfusion(17) = SK_MY[0]*(P(17,20) + P(17,0)*SH_MAG[2] + P(17,1)*SH_MAG[1] + P(17,2)*SH_MAG[0] - P(17,3)*SK_MY[2] - P(17,17)*SK_MY[1] - P(17,16)*SK_MY[3] + P(17,18)*SK_MY[4]);
        Kfusion(18) = SK_MY[0]*(P(18,20) + P(18,0)*SH_MAG[2] + P(18,1)*SH_MAG[1] + P(18,2)*SH_MAG[0] - P(18,3)*SK_MY[2] - P(18,17)*SK_MY[1] - P(18,16)*SK_MY[3] + P(18,18)*SK_MY[4]);
        Kfusion(19) = SK_MY[0]*(P(19,20) + P(19,0)*SH_MAG[2] + P(19,1)*SH_MAG[1] + P(19,2)*SH_MAG[0] - P(19,3)*SK_MY[2] - P(19,17)*SK_MY[1] - P(19,16)*SK_MY[3] + P(19,18)*SK_MY[4]);
        Kfusion(20) = SK_MY[0]*(P(20,20) + P(20,0)*SH_MAG[2] + P(20,1)*SH_MAG[1] + P(20,2)*SH_MAG[0] - P(20,3)*SK_MY[2] - P(20,17)*SK_MY[1] - P(20,16)*SK_MY[3] + P(20,18)*SK_MY[4]);
        Kfusion(21) = SK_MY[0]*(P(21,20) + P(21,0)*SH_MAG[2] + P(21,1)*SH_MAG[1] + P(21,2)*SH_MAG[0] - P(21,3)*SK_MY[2] - P(21,17)*SK_MY[1] - P(21,16)*SK_MY[3] + P(21,18)*SK_MY[4]);

        // find largest observation variance difference as a fraction of the matlab value
        float max_diff_fraction = 0.0f;
        int max_row;
        float max_old, max_new;
        for (int row=0; row<24; row++) {
            float diff_fraction;
            if (H_MAG(row) != 0.0f) {
                diff_fraction = fabsf(Hfusion_sympy[row] - H_MAG(row)) / fabsf(H_MAG(row));
            } else if (Hfusion_sympy[row] != 0.0f) {
                diff_fraction = fabsf(Hfusion_sympy[row] - H_MAG(row)) / fabsf(Hfusion_sympy[row]);
            } else {
                diff_fraction = 0.0f;
            }
            if (diff_fraction > max_diff_fraction) {
                max_diff_fraction = diff_fraction;
                max_row = row;
                max_old = H_MAG(row);
                max_new = Hfusion_sympy[row];
            }
        }

        printf("Y axis Hfusion max_diff_fraction = %e , old = %e , new = %e , location index = %i\n",max_diff_fraction, max_old, max_new, max_row);

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
            if (diff_fraction > max_diff_fraction) {
                max_diff_fraction = diff_fraction;
                max_row = row;
                max_old = Kfusion(row);
                max_new = Kfusion_sympy(row);
            }
        }

        printf("Y axis Kfusion max_diff_fraction = %e , old = %e , new = %e , location index = %i\n",max_diff_fraction, max_old, max_new, max_row);

    }

    // Compare Z axis equations
    {
        mag_innov_var = (HK29*HK81 + HK30*HK82 + HK32*HK76 - HK34*HK79 + HK73*HK77 + HK74*HK80 - HK75*HK78 + HK83 + R_MAG);
        float HK84 = 1.0F/mag_innov_var;

        // calculate Z axis observation jacobians
        memset(Hfusion, 0, sizeof(Hfusion));
        Hfusion[0] = HK51;
        Hfusion[1] = -2*HK10 - 2*HK11 + 2*HK12;
        Hfusion[2] = HK4;
        Hfusion[3] = HK6;
        Hfusion[16] = 2*HK72;
        Hfusion[17] = -2*HK54 + 2*HK55;
        Hfusion[18] = HK73;
        Hfusion[21] = 1;

        // Calculate Z axis Kalman gains
        if (update_all_states) {
            Kfusion(0) = HK76*HK84;
            Kfusion(1) = HK79*HK84;
            Kfusion(2) = HK82*HK84;
            Kfusion(3) = HK81*HK84;

            for (unsigned row = 4; row <= 15; row++) {
                Kfusion(row) = HK84*(HK29*P(3,row) + HK30*P(2,row) + HK32*P(0,row) - HK34*P(1,row) + HK73*P(row,18) + HK74*P(row,16) - HK75*P(row,17) + P(row,21));
            }

            for (unsigned row = 22; row <= 23; row++) {
                Kfusion(row) = HK84*(HK29*P(3,row) + HK30*P(2,row) + HK32*P(0,row) - HK34*P(1,row) + HK73*P(row,18) + HK74*P(row,16) - HK75*P(row,17) + P(row,21));
            }
        } else {
            for (uint8_t i = 0; i < 16; i++) {
                Kfusion(i) = 0.0f;
            }

            Kfusion(22) = 0.0f;
            Kfusion(23) = 0.0f;
        }

        Kfusion(16) = HK80*HK84;
        Kfusion(17) = HK78*HK84;
        Kfusion(18) = HK77*HK84;

        for (unsigned row = 19; row <= 20; row++) {
        Kfusion(row) = HK84*(HK29*P(3,row) + HK30*P(2,row) + HK32*P(0,row) - HK34*P(1,row) + HK73*P(18,row) + HK74*P(16,row) - HK75*P(17,row) + P(row,21));
        }

        Kfusion(21) = HK83*HK84;

        // save output and repeat calculation using legacy matlab generated code
        float Hfusion_sympy[24];
        Vector24f Kfusion_sympy;
        for (int row=0; row<24; row++) {
            Hfusion_sympy[row] = Hfusion[row];
            Kfusion_sympy(row) = Kfusion(row);
        }

        // repeat calculation using matlab generated equations

        // Z axis innovation variance
        mag_innov_var = (P(21,21) + R_MAG + P(0,21)*SH_MAG[1] - P(1,21)*SH_MAG[2] + P(3,21)*SH_MAG[0] + P(18,21)*(SH_MAG[3] - SH_MAG[4] - SH_MAG[5] + SH_MAG[6]) + (2.0f*q0*q2 + 2.0f*q1*q3)*(P(21,16) + P(0,16)*SH_MAG[1] - P(1,16)*SH_MAG[2] + P(3,16)*SH_MAG[0] + P(18,16)*(SH_MAG[3] - SH_MAG[4] - SH_MAG[5] + SH_MAG[6]) + P(16,16)*(2.0f*q0*q2 + 2.0f*q1*q3) - P(17,16)*(2.0f*q0*q1 - 2.0f*q2*q3) + P(2,16)*(SH_MAG[7] + SH_MAG[8] - 2.0f*magD*q2)) - (2.0f*q0*q1 - 2.0f*q2*q3)*(P(21,17) + P(0,17)*SH_MAG[1] - P(1,17)*SH_MAG[2] + P(3,17)*SH_MAG[0] + P(18,17)*(SH_MAG[3] - SH_MAG[4] - SH_MAG[5] + SH_MAG[6]) + P(16,17)*(2.0f*q0*q2 + 2.0f*q1*q3) - P(17,17)*(2.0f*q0*q1 - 2.0f*q2*q3) + P(2,17)*(SH_MAG[7] + SH_MAG[8] - 2.0f*magD*q2)) + (SH_MAG[7] + SH_MAG[8] - 2.0f*magD*q2)*(P(21,2) + P(0,2)*SH_MAG[1] - P(1,2)*SH_MAG[2] + P(3,2)*SH_MAG[0] + P(18,2)*(SH_MAG[3] - SH_MAG[4] - SH_MAG[5] + SH_MAG[6]) + P(16,2)*(2.0f*q0*q2 + 2.0f*q1*q3) - P(17,2)*(2.0f*q0*q1 - 2.0f*q2*q3) + P(2,2)*(SH_MAG[7] + SH_MAG[8] - 2.0f*magD*q2)) + P(16,21)*(2.0f*q0*q2 + 2.0f*q1*q3) - P(17,21)*(2.0f*q0*q1 - 2.0f*q2*q3) + SH_MAG[1]*(P(21,0) + P(0,0)*SH_MAG[1] - P(1,0)*SH_MAG[2] + P(3,0)*SH_MAG[0] + P(18,0)*(SH_MAG[3] - SH_MAG[4] - SH_MAG[5] + SH_MAG[6]) + P(16,0)*(2.0f*q0*q2 + 2.0f*q1*q3) - P(17,0)*(2.0f*q0*q1 - 2.0f*q2*q3) + P(2,0)*(SH_MAG[7] + SH_MAG[8] - 2.0f*magD*q2)) - SH_MAG[2]*(P(21,1) + P(0,1)*SH_MAG[1] - P(1,1)*SH_MAG[2] + P(3,1)*SH_MAG[0] + P(18,1)*(SH_MAG[3] - SH_MAG[4] - SH_MAG[5] + SH_MAG[6]) + P(16,1)*(2.0f*q0*q2 + 2.0f*q1*q3) - P(17,1)*(2.0f*q0*q1 - 2.0f*q2*q3) + P(2,1)*(SH_MAG[7] + SH_MAG[8] - 2.0f*magD*q2)) + SH_MAG[0]*(P(21,3) + P(0,3)*SH_MAG[1] - P(1,3)*SH_MAG[2] + P(3,3)*SH_MAG[0] + P(18,3)*(SH_MAG[3] - SH_MAG[4] - SH_MAG[5] + SH_MAG[6]) + P(16,3)*(2.0f*q0*q2 + 2.0f*q1*q3) - P(17,3)*(2.0f*q0*q1 - 2.0f*q2*q3) + P(2,3)*(SH_MAG[7] + SH_MAG[8] - 2.0f*magD*q2)) + (SH_MAG[3] - SH_MAG[4] - SH_MAG[5] + SH_MAG[6])*(P(21,18) + P(0,18)*SH_MAG[1] - P(1,18)*SH_MAG[2] + P(3,18)*SH_MAG[0] + P(18,18)*(SH_MAG[3] - SH_MAG[4] - SH_MAG[5] + SH_MAG[6]) + P(16,18)*(2.0f*q0*q2 + 2.0f*q1*q3) - P(17,18)*(2.0f*q0*q1 - 2.0f*q2*q3) + P(2,18)*(SH_MAG[7] + SH_MAG[8] - 2.0f*magD*q2)) + P(2,21)*(SH_MAG[7] + SH_MAG[8] - 2.0f*magD*q2));

        // calculate Z axis observation jacobians
        H_MAG.setZero();
        H_MAG(0) = SH_MAG[1];
        H_MAG(1) = -SH_MAG[2];
        H_MAG(2) = SH_MAG[7] + SH_MAG[8] - 2.0f*magD*q2;
        H_MAG(3) = SH_MAG[0];
        H_MAG(16) = 2.0f*q0*q2 + 2.0f*q1*q3;
        H_MAG(17) = 2.0f*q2*q3 - 2.0f*q0*q1;
        H_MAG(18) = SH_MAG[3] - SH_MAG[4] - SH_MAG[5] + SH_MAG[6];
        H_MAG(21) = 1.0f;

        // Calculate Z axis Kalman gains
        float SK_MZ[5];
        SK_MZ[0] = 1.0f / mag_innov_var;
        SK_MZ[1] = SH_MAG[3] - SH_MAG[4] - SH_MAG[5] + SH_MAG[6];
        SK_MZ[2] = SH_MAG[7] + SH_MAG[8] - 2.0f*magD*q2;
        SK_MZ[3] = 2.0f*q0*q1 - 2.0f*q2*q3;
        SK_MZ[4] = 2.0f*q0*q2 + 2.0f*q1*q3;

        if (update_all_states) {
            Kfusion(0) = SK_MZ[0]*(P(0,21) + P(0,0)*SH_MAG[1] - P(0,1)*SH_MAG[2] + P(0,3)*SH_MAG[0] + P(0,2)*SK_MZ[2] + P(0,18)*SK_MZ[1] + P(0,16)*SK_MZ[4] - P(0,17)*SK_MZ[3]);
            Kfusion(1) = SK_MZ[0]*(P(1,21) + P(1,0)*SH_MAG[1] - P(1,1)*SH_MAG[2] + P(1,3)*SH_MAG[0] + P(1,2)*SK_MZ[2] + P(1,18)*SK_MZ[1] + P(1,16)*SK_MZ[4] - P(1,17)*SK_MZ[3]);
            Kfusion(2) = SK_MZ[0]*(P(2,21) + P(2,0)*SH_MAG[1] - P(2,1)*SH_MAG[2] + P(2,3)*SH_MAG[0] + P(2,2)*SK_MZ[2] + P(2,18)*SK_MZ[1] + P(2,16)*SK_MZ[4] - P(2,17)*SK_MZ[3]);
            Kfusion(3) = SK_MZ[0]*(P(3,21) + P(3,0)*SH_MAG[1] - P(3,1)*SH_MAG[2] + P(3,3)*SH_MAG[0] + P(3,2)*SK_MZ[2] + P(3,18)*SK_MZ[1] + P(3,16)*SK_MZ[4] - P(3,17)*SK_MZ[3]);
            Kfusion(4) = SK_MZ[0]*(P(4,21) + P(4,0)*SH_MAG[1] - P(4,1)*SH_MAG[2] + P(4,3)*SH_MAG[0] + P(4,2)*SK_MZ[2] + P(4,18)*SK_MZ[1] + P(4,16)*SK_MZ[4] - P(4,17)*SK_MZ[3]);
            Kfusion(5) = SK_MZ[0]*(P(5,21) + P(5,0)*SH_MAG[1] - P(5,1)*SH_MAG[2] + P(5,3)*SH_MAG[0] + P(5,2)*SK_MZ[2] + P(5,18)*SK_MZ[1] + P(5,16)*SK_MZ[4] - P(5,17)*SK_MZ[3]);
            Kfusion(6) = SK_MZ[0]*(P(6,21) + P(6,0)*SH_MAG[1] - P(6,1)*SH_MAG[2] + P(6,3)*SH_MAG[0] + P(6,2)*SK_MZ[2] + P(6,18)*SK_MZ[1] + P(6,16)*SK_MZ[4] - P(6,17)*SK_MZ[3]);
            Kfusion(7) = SK_MZ[0]*(P(7,21) + P(7,0)*SH_MAG[1] - P(7,1)*SH_MAG[2] + P(7,3)*SH_MAG[0] + P(7,2)*SK_MZ[2] + P(7,18)*SK_MZ[1] + P(7,16)*SK_MZ[4] - P(7,17)*SK_MZ[3]);
            Kfusion(8) = SK_MZ[0]*(P(8,21) + P(8,0)*SH_MAG[1] - P(8,1)*SH_MAG[2] + P(8,3)*SH_MAG[0] + P(8,2)*SK_MZ[2] + P(8,18)*SK_MZ[1] + P(8,16)*SK_MZ[4] - P(8,17)*SK_MZ[3]);
            Kfusion(9) = SK_MZ[0]*(P(9,21) + P(9,0)*SH_MAG[1] - P(9,1)*SH_MAG[2] + P(9,3)*SH_MAG[0] + P(9,2)*SK_MZ[2] + P(9,18)*SK_MZ[1] + P(9,16)*SK_MZ[4] - P(9,17)*SK_MZ[3]);
            Kfusion(10) = SK_MZ[0]*(P(10,21) + P(10,0)*SH_MAG[1] - P(10,1)*SH_MAG[2] + P(10,3)*SH_MAG[0] + P(10,2)*SK_MZ[2] + P(10,18)*SK_MZ[1] + P(10,16)*SK_MZ[4] - P(10,17)*SK_MZ[3]);
            Kfusion(11) = SK_MZ[0]*(P(11,21) + P(11,0)*SH_MAG[1] - P(11,1)*SH_MAG[2] + P(11,3)*SH_MAG[0] + P(11,2)*SK_MZ[2] + P(11,18)*SK_MZ[1] + P(11,16)*SK_MZ[4] - P(11,17)*SK_MZ[3]);
            Kfusion(12) = SK_MZ[0]*(P(12,21) + P(12,0)*SH_MAG[1] - P(12,1)*SH_MAG[2] + P(12,3)*SH_MAG[0] + P(12,2)*SK_MZ[2] + P(12,18)*SK_MZ[1] + P(12,16)*SK_MZ[4] - P(12,17)*SK_MZ[3]);
            Kfusion(13) = SK_MZ[0]*(P(13,21) + P(13,0)*SH_MAG[1] - P(13,1)*SH_MAG[2] + P(13,3)*SH_MAG[0] + P(13,2)*SK_MZ[2] + P(13,18)*SK_MZ[1] + P(13,16)*SK_MZ[4] - P(13,17)*SK_MZ[3]);
            Kfusion(14) = SK_MZ[0]*(P(14,21) + P(14,0)*SH_MAG[1] - P(14,1)*SH_MAG[2] + P(14,3)*SH_MAG[0] + P(14,2)*SK_MZ[2] + P(14,18)*SK_MZ[1] + P(14,16)*SK_MZ[4] - P(14,17)*SK_MZ[3]);
            Kfusion(15) = SK_MZ[0]*(P(15,21) + P(15,0)*SH_MAG[1] - P(15,1)*SH_MAG[2] + P(15,3)*SH_MAG[0] + P(15,2)*SK_MZ[2] + P(15,18)*SK_MZ[1] + P(15,16)*SK_MZ[4] - P(15,17)*SK_MZ[3]);
            Kfusion(22) = SK_MZ[0]*(P(22,21) + P(22,0)*SH_MAG[1] - P(22,1)*SH_MAG[2] + P(22,3)*SH_MAG[0] + P(22,2)*SK_MZ[2] + P(22,18)*SK_MZ[1] + P(22,16)*SK_MZ[4] - P(22,17)*SK_MZ[3]);
            Kfusion(23) = SK_MZ[0]*(P(23,21) + P(23,0)*SH_MAG[1] - P(23,1)*SH_MAG[2] + P(23,3)*SH_MAG[0] + P(23,2)*SK_MZ[2] + P(23,18)*SK_MZ[1] + P(23,16)*SK_MZ[4] - P(23,17)*SK_MZ[3]);
        }

        Kfusion(16) = SK_MZ[0]*(P(16,21) + P(16,0)*SH_MAG[1] - P(16,1)*SH_MAG[2] + P(16,3)*SH_MAG[0] + P(16,2)*SK_MZ[2] + P(16,18)*SK_MZ[1] + P(16,16)*SK_MZ[4] - P(16,17)*SK_MZ[3]);
        Kfusion(17) = SK_MZ[0]*(P(17,21) + P(17,0)*SH_MAG[1] - P(17,1)*SH_MAG[2] + P(17,3)*SH_MAG[0] + P(17,2)*SK_MZ[2] + P(17,18)*SK_MZ[1] + P(17,16)*SK_MZ[4] - P(17,17)*SK_MZ[3]);
        Kfusion(18) = SK_MZ[0]*(P(18,21) + P(18,0)*SH_MAG[1] - P(18,1)*SH_MAG[2] + P(18,3)*SH_MAG[0] + P(18,2)*SK_MZ[2] + P(18,18)*SK_MZ[1] + P(18,16)*SK_MZ[4] - P(18,17)*SK_MZ[3]);
        Kfusion(19) = SK_MZ[0]*(P(19,21) + P(19,0)*SH_MAG[1] - P(19,1)*SH_MAG[2] + P(19,3)*SH_MAG[0] + P(19,2)*SK_MZ[2] + P(19,18)*SK_MZ[1] + P(19,16)*SK_MZ[4] - P(19,17)*SK_MZ[3]);
        Kfusion(20) = SK_MZ[0]*(P(20,21) + P(20,0)*SH_MAG[1] - P(20,1)*SH_MAG[2] + P(20,3)*SH_MAG[0] + P(20,2)*SK_MZ[2] + P(20,18)*SK_MZ[1] + P(20,16)*SK_MZ[4] - P(20,17)*SK_MZ[3]);
        Kfusion(21) = SK_MZ[0]*(P(21,21) + P(21,0)*SH_MAG[1] - P(21,1)*SH_MAG[2] + P(21,3)*SH_MAG[0] + P(21,2)*SK_MZ[2] + P(21,18)*SK_MZ[1] + P(21,16)*SK_MZ[4] - P(21,17)*SK_MZ[3]);

        // find largest observation variance difference as a fraction of the matlab value
        float max_diff_fraction = 0.0f;
        int max_row;
        float max_old, max_new;
        for (int row=0; row<24; row++) {
            float diff_fraction;
            if (H_MAG(row) != 0.0f) {
                diff_fraction = fabsf(Hfusion_sympy[row] - H_MAG(row)) / fabsf(H_MAG(row));
            } else if (Hfusion_sympy[row] != 0.0f) {
                diff_fraction = fabsf(Hfusion_sympy[row] - H_MAG(row)) / fabsf(Hfusion_sympy[row]);
            } else {
                diff_fraction = 0.0f;
            }
            if (diff_fraction > max_diff_fraction) {
                max_diff_fraction = diff_fraction;
                max_row = row;
                max_old = H_MAG(row);
                max_new = Hfusion_sympy[row];
            }
        }

        printf("Z axis Hfusion max_diff_fraction = %e , old = %e , new = %e , location index = %i\n",max_diff_fraction, max_old, max_new, max_row);

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
            if (diff_fraction > max_diff_fraction) {
                max_diff_fraction = diff_fraction;
                max_row = row;
                max_old = Kfusion(row);
                max_new = Kfusion_sympy(row);
            }
        }

        printf("Z axis Kfusion max_diff_fraction = %e , old = %e , new = %e , location index = %i\n",max_diff_fraction, max_old, max_new, max_row);

    }

    return 0;
}
