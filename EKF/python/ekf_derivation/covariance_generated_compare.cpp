#include <math.h>
#include <stdio.h>
#include <cstdlib>
#include "../../../../matrix/matrix/math.hpp"

typedef matrix::SquareMatrix<float, 24> SquareMatrix24f;

inline float sq(float in) {
    return in * in;
}

namespace ecl{
	inline float powf(float x, int exp)
	{
		float ret;
		if (exp > 0) {
			ret = x;
			for (int count = 1; count < exp; count++) {
				ret *= x;
			}
			return ret;
		} else if (exp < 0) {
			ret = x;
			for (int count = -1; count > exp; count--) {
				ret *= x;
			}
			return 1.0f / ret;
		} else {
			return 1.0f;
		}
	}
}

int main()
{
    // test integer power function
    const float x_in = 2.0f;
    printf("ecl::powf(0.5,n) = %e , %e , %e , %e , %e\n",(double)ecl::powf(x_in,-2),
                                    (double)ecl::powf(x_in,-1),
                                    (double)ecl::powf(x_in, 0),
                                    (double)ecl::powf(x_in, 1),
                                    (double)ecl::powf(x_in, 2));
    printf("     powf(0.5,n) = %e , %e , %e , %e , %e\n\n",(double)powf(x_in,-2),
                                    (double)powf(x_in,-1),
                                    (double)powf(x_in, 0),
                                    (double)powf(x_in, 1),
                                    (double)powf(x_in, 2));
    // create input data
    const float dt = 0.01f;

    for (uint8_t iteration=0; iteration<100; iteration++) {
        // quaternion states must be normalised
        float q0 = 2.0f * ((float)rand() - 0.5f);
        float q1 = 2.0f * ((float)rand() - 0.5f);
        float q2 = 2.0f * ((float)rand() - 0.5f);
        float q3 = 2.0f * ((float)rand() - 0.5f);
        const float length = sqrtf(sq(q0) + sq(q1) + sq(q2) + sq(q3));
        q0 /= length;
        q1 /= length;
        q2 /= length;
        q3 /= length;

        // up to 1 rad/sec of rate
        const float dax = 2.0f * dt * ((float)rand() - 0.5f);
        const float day = 2.0f * dt * ((float)rand() - 0.5f);
        const float daz = 2.0f * dt * ((float)rand() - 0.5f);

        // up to 2g of accel
        const float dvx = 2.0f * 20.0f * dt * ((float)rand() - 0.5f);
        const float dvy = 2.0f * 20.0f * dt * ((float)rand() - 0.5f);
        const float dvz = 2.0f * 20.0f * dt * ((float)rand() - 0.5f);

        // up to 0.1 rad/sec of gyro bias
        const float dax_b = 2.0f * 0.1f * dt * ((float)rand() - 0.5f);
        const float day_b = 2.0f * 0.1f * dt * ((float)rand() - 0.5f);
        const float daz_b = 2.0f * 0.1f * dt * ((float)rand() - 0.5f);

        // up to 0.5 m/s/s of accel bias
        const float dvx_b = 2.0f * 0.5f * dt * ((float)rand() - 0.5f);
        const float dvy_b = 2.0f * 0.5f * dt * ((float)rand() - 0.5f);
        const float dvz_b = 2.0f * 0.5f * dt * ((float)rand() - 0.5f);


        const float daxVar = sq(dt * 0.015f);
        const float dayVar = daxVar;
        const float dazVar = daxVar;

        const float dvxVar = sq(dt * 0.3f);
        const float dvyVar = dvxVar;
        const float dvzVar = dvxVar;

        // create a symmetrical positive dfinite matrix with off diagonals between -1 and 1 and diagonals between 0 and 1
        SquareMatrix24f P;
        for (int col=0; col<=23; col++) {
            for (int row=0; row<=col; row++) {
                if (row == col) {
                    P(row,col) = (float)rand();
                } else {
                    P(col,row) = P(row,col) = 0.01f * 2.0f * ((float)rand() - 0.5f);
                }
            }
        }

        // Equations for covariance matrix prediction, without process noise!
        const float PS0 = ecl::powf(q1, 2);
        const float PS1 = 0.25F*daxVar;
        const float PS2 = ecl::powf(q2, 2);
        const float PS3 = 0.25F*dayVar;
        const float PS4 = ecl::powf(q3, 2);
        const float PS5 = 0.25F*dazVar;
        const float PS6 = 0.5F*q1;
        const float PS7 = 0.5F*q2;
        const float PS8 = P(10,11)*PS7;
        const float PS9 = 0.5F*q3;
        const float PS10 = P(10,12)*PS9;
        const float PS11 = 0.5F*dax - 0.5F*dax_b;
        const float PS12 = 0.5F*day - 0.5F*day_b;
        const float PS13 = 0.5F*daz - 0.5F*daz_b;
        const float PS14 = P(0,10) - P(1,10)*PS11 + P(10,10)*PS6 - P(2,10)*PS12 - P(3,10)*PS13 + PS10 + PS8;
        const float PS15 = P(10,11)*PS6;
        const float PS16 = P(11,12)*PS9;
        const float PS17 = P(0,11) - P(1,11)*PS11 + P(11,11)*PS7 - P(2,11)*PS12 - P(3,11)*PS13 + PS15 + PS16;
        const float PS18 = P(10,12)*PS6;
        const float PS19 = P(11,12)*PS7;
        const float PS20 = P(0,12) - P(1,12)*PS11 + P(12,12)*PS9 - P(2,12)*PS12 - P(3,12)*PS13 + PS18 + PS19;
        const float PS21 = P(1,2)*PS12;
        const float PS22 = -P(1,3)*PS13;
        const float PS23 = P(0,1) - P(1,1)*PS11 + P(1,10)*PS6 + P(1,11)*PS7 + P(1,12)*PS9 - PS21 + PS22;
        const float PS24 = -P(1,2)*PS11;
        const float PS25 = P(2,3)*PS13;
        const float PS26 = P(0,2) + P(2,10)*PS6 + P(2,11)*PS7 + P(2,12)*PS9 - P(2,2)*PS12 + PS24 - PS25;
        const float PS27 = P(1,3)*PS11;
        const float PS28 = -P(2,3)*PS12;
        const float PS29 = P(0,3) + P(3,10)*PS6 + P(3,11)*PS7 + P(3,12)*PS9 - P(3,3)*PS13 - PS27 + PS28;
        const float PS30 = P(0,1)*PS11;
        const float PS31 = P(0,2)*PS12;
        const float PS32 = P(0,3)*PS13;
        const float PS33 = P(0,0) + P(0,10)*PS6 + P(0,11)*PS7 + P(0,12)*PS9 - PS30 - PS31 - PS32;
        const float PS34 = 0.5F*q0;
        const float PS35 = q2*q3;
        const float PS36 = q0*q1;
        const float PS37 = q1*q3;
        const float PS38 = q0*q2;
        const float PS39 = q1*q2;
        const float PS40 = q0*q3;
        const float PS41 = -PS2;
        const float PS42 = ecl::powf(q0, 2);
        const float PS43 = -PS4 + PS42;
        const float PS44 = PS0 + PS41 + PS43;
        const float PS45 = P(0,13) - P(1,13)*PS11 + P(10,13)*PS6 + P(11,13)*PS7 + P(12,13)*PS9 - P(2,13)*PS12 - P(3,13)*PS13;
        const float PS46 = PS37 + PS38;
        const float PS47 = P(0,15) - P(1,15)*PS11 + P(10,15)*PS6 + P(11,15)*PS7 + P(12,15)*PS9 - P(2,15)*PS12 - P(3,15)*PS13;
        const float PS48 = 2*PS47;
        const float PS49 = dvy - dvy_b;
        const float PS50 = dvx - dvx_b;
        const float PS51 = dvz - dvz_b;
        const float PS52 = PS49*q0 + PS50*q3 - PS51*q1;
        const float PS53 = 2*PS29;
        const float PS54 = -PS39 + PS40;
        const float PS55 = P(0,14) - P(1,14)*PS11 + P(10,14)*PS6 + P(11,14)*PS7 + P(12,14)*PS9 - P(2,14)*PS12 - P(3,14)*PS13;
        const float PS56 = 2*PS55;
        const float PS57 = -PS49*q3 + PS50*q0 + PS51*q2;
        const float PS58 = 2*PS33;
        const float PS59 = PS49*q1 - PS50*q2 + PS51*q0;
        const float PS60 = 2*PS59;
        const float PS61 = PS49*q2 + PS50*q1 + PS51*q3;
        const float PS62 = 2*PS61;
        const float PS63 = P(0,4) - P(1,4)*PS11 - P(2,4)*PS12 - P(3,4)*PS13 + P(4,10)*PS6 + P(4,11)*PS7 + P(4,12)*PS9;
        const float PS64 = -PS0;
        const float PS65 = PS2 + PS43 + PS64;
        const float PS66 = PS39 + PS40;
        const float PS67 = 2*PS45;
        const float PS68 = -PS35 + PS36;
        const float PS69 = P(0,5) - P(1,5)*PS11 - P(2,5)*PS12 - P(3,5)*PS13 + P(5,10)*PS6 + P(5,11)*PS7 + P(5,12)*PS9;
        const float PS70 = PS4 + PS41 + PS42 + PS64;
        const float PS71 = PS35 + PS36;
        const float PS72 = 2*PS57;
        const float PS73 = -PS37 + PS38;
        const float PS74 = 2*PS52;
        const float PS75 = P(0,6) - P(1,6)*PS11 - P(2,6)*PS12 - P(3,6)*PS13 + P(6,10)*PS6 + P(6,11)*PS7 + P(6,12)*PS9;
        const float PS76 = -P(10,11)*PS34;
        const float PS77 = P(0,11)*PS11 + P(1,11) + P(11,11)*PS9 + P(2,11)*PS13 - P(3,11)*PS12 - PS19 + PS76;
        const float PS78 = P(0,2)*PS13;
        const float PS79 = P(0,3)*PS12;
        const float PS80 = P(0,0)*PS11 + P(0,1) - P(0,10)*PS34 + P(0,11)*PS9 - P(0,12)*PS7 + PS78 - PS79;
        const float PS81 = P(0,2)*PS11;
        const float PS82 = P(1,2) - P(2,10)*PS34 + P(2,11)*PS9 - P(2,12)*PS7 + P(2,2)*PS13 + PS28 + PS81;
        const float PS83 = P(10,11)*PS9;
        const float PS84 = P(10,12)*PS7;
        const float PS85 = P(0,10)*PS11 + P(1,10) - P(10,10)*PS34 + P(2,10)*PS13 - P(3,10)*PS12 + PS83 - PS84;
        const float PS86 = -P(10,12)*PS34;
        const float PS87 = P(0,12)*PS11 + P(1,12) - P(12,12)*PS7 + P(2,12)*PS13 - P(3,12)*PS12 + PS16 + PS86;
        const float PS88 = P(0,3)*PS11;
        const float PS89 = P(1,3) - P(3,10)*PS34 + P(3,11)*PS9 - P(3,12)*PS7 - P(3,3)*PS12 + PS25 + PS88;
        const float PS90 = P(1,2)*PS13;
        const float PS91 = P(1,3)*PS12;
        const float PS92 = P(1,1) - P(1,10)*PS34 + P(1,11)*PS9 - P(1,12)*PS7 + PS30 + PS90 - PS91;
        const float PS93 = P(0,13)*PS11 + P(1,13) - P(10,13)*PS34 + P(11,13)*PS9 - P(12,13)*PS7 + P(2,13)*PS13 - P(3,13)*PS12;
        const float PS94 = P(0,15)*PS11 + P(1,15) - P(10,15)*PS34 + P(11,15)*PS9 - P(12,15)*PS7 + P(2,15)*PS13 - P(3,15)*PS12;
        const float PS95 = 2*PS94;
        const float PS96 = P(0,14)*PS11 + P(1,14) - P(10,14)*PS34 + P(11,14)*PS9 - P(12,14)*PS7 + P(2,14)*PS13 - P(3,14)*PS12;
        const float PS97 = 2*PS96;
        const float PS98 = P(0,4)*PS11 + P(1,4) + P(2,4)*PS13 - P(3,4)*PS12 - P(4,10)*PS34 + P(4,11)*PS9 - P(4,12)*PS7;
        const float PS99 = 2*PS93;
        const float PS100 = P(0,5)*PS11 + P(1,5) + P(2,5)*PS13 - P(3,5)*PS12 - P(5,10)*PS34 + P(5,11)*PS9 - P(5,12)*PS7;
        const float PS101 = P(0,6)*PS11 + P(1,6) + P(2,6)*PS13 - P(3,6)*PS12 - P(6,10)*PS34 + P(6,11)*PS9 - P(6,12)*PS7;
        const float PS102 = -P(11,12)*PS34;
        const float PS103 = P(0,12)*PS12 - P(1,12)*PS13 + P(12,12)*PS6 + P(2,12) + P(3,12)*PS11 - PS10 + PS102;
        const float PS104 = P(2,3) - P(3,10)*PS9 - P(3,11)*PS34 + P(3,12)*PS6 + P(3,3)*PS11 + PS22 + PS79;
        const float PS105 = P(0,1)*PS13;
        const float PS106 = P(0,0)*PS12 - P(0,10)*PS9 - P(0,11)*PS34 + P(0,12)*PS6 + P(0,2) - PS105 + PS88;
        const float PS107 = P(11,12)*PS6;
        const float PS108 = P(0,11)*PS12 - P(1,11)*PS13 - P(11,11)*PS34 + P(2,11) + P(3,11)*PS11 + PS107 - PS83;
        const float PS109 = P(0,10)*PS12 - P(1,10)*PS13 - P(10,10)*PS9 + P(2,10) + P(3,10)*PS11 + PS18 + PS76;
        const float PS110 = P(0,1)*PS12;
        const float PS111 = -P(1,1)*PS13 - P(1,10)*PS9 - P(1,11)*PS34 + P(1,12)*PS6 + P(1,2) + PS110 + PS27;
        const float PS112 = P(2,3)*PS11;
        const float PS113 = -P(2,10)*PS9 - P(2,11)*PS34 + P(2,12)*PS6 + P(2,2) + PS112 + PS31 - PS90;
        const float PS114 = P(0,13)*PS12 - P(1,13)*PS13 - P(10,13)*PS9 - P(11,13)*PS34 + P(12,13)*PS6 + P(2,13) + P(3,13)*PS11;
        const float PS115 = P(0,15)*PS12 - P(1,15)*PS13 - P(10,15)*PS9 - P(11,15)*PS34 + P(12,15)*PS6 + P(2,15) + P(3,15)*PS11;
        const float PS116 = 2*PS115;
        const float PS117 = P(0,14)*PS12 - P(1,14)*PS13 - P(10,14)*PS9 - P(11,14)*PS34 + P(12,14)*PS6 + P(2,14) + P(3,14)*PS11;
        const float PS118 = 2*PS117;
        const float PS119 = P(0,4)*PS12 - P(1,4)*PS13 + P(2,4) + P(3,4)*PS11 - P(4,10)*PS9 - P(4,11)*PS34 + P(4,12)*PS6;
        const float PS120 = 2*PS114;
        const float PS121 = P(0,5)*PS12 - P(1,5)*PS13 + P(2,5) + P(3,5)*PS11 - P(5,10)*PS9 - P(5,11)*PS34 + P(5,12)*PS6;
        const float PS122 = P(0,6)*PS12 - P(1,6)*PS13 + P(2,6) + P(3,6)*PS11 - P(6,10)*PS9 - P(6,11)*PS34 + P(6,12)*PS6;
        const float PS123 = P(0,10)*PS13 + P(1,10)*PS12 + P(10,10)*PS7 - P(2,10)*PS11 + P(3,10) - PS15 + PS86;
        const float PS124 = P(1,1)*PS12 + P(1,10)*PS7 - P(1,11)*PS6 - P(1,12)*PS34 + P(1,3) + PS105 + PS24;
        const float PS125 = P(0,0)*PS13 + P(0,10)*PS7 - P(0,11)*PS6 - P(0,12)*PS34 + P(0,3) + PS110 - PS81;
        const float PS126 = P(0,12)*PS13 + P(1,12)*PS12 - P(12,12)*PS34 - P(2,12)*PS11 + P(3,12) - PS107 + PS84;
        const float PS127 = P(0,11)*PS13 + P(1,11)*PS12 - P(11,11)*PS6 - P(2,11)*PS11 + P(3,11) + PS102 + PS8;
        const float PS128 = P(2,10)*PS7 - P(2,11)*PS6 - P(2,12)*PS34 - P(2,2)*PS11 + P(2,3) + PS21 + PS78;
        const float PS129 = P(3,10)*PS7 - P(3,11)*PS6 - P(3,12)*PS34 + P(3,3) - PS112 + PS32 + PS91;
        const float PS130 = P(0,13)*PS13 + P(1,13)*PS12 + P(10,13)*PS7 - P(11,13)*PS6 - P(12,13)*PS34 - P(2,13)*PS11 + P(3,13);
        const float PS131 = P(0,15)*PS13 + P(1,15)*PS12 + P(10,15)*PS7 - P(11,15)*PS6 - P(12,15)*PS34 - P(2,15)*PS11 + P(3,15);
        const float PS132 = 2*PS131;
        const float PS133 = P(0,14)*PS13 + P(1,14)*PS12 + P(10,14)*PS7 - P(11,14)*PS6 - P(12,14)*PS34 - P(2,14)*PS11 + P(3,14);
        const float PS134 = 2*PS133;
        const float PS135 = P(0,4)*PS13 + P(1,4)*PS12 - P(2,4)*PS11 + P(3,4) + P(4,10)*PS7 - P(4,11)*PS6 - P(4,12)*PS34;
        const float PS136 = 2*PS130;
        const float PS137 = P(0,5)*PS13 + P(1,5)*PS12 - P(2,5)*PS11 + P(3,5) + P(5,10)*PS7 - P(5,11)*PS6 - P(5,12)*PS34;
        const float PS138 = P(0,6)*PS13 + P(1,6)*PS12 - P(2,6)*PS11 + P(3,6) + P(6,10)*PS7 - P(6,11)*PS6 - P(6,12)*PS34;
        const float PS139 = 2*PS46;
        const float PS140 = 2*PS54;
        const float PS141 = P(0,13)*PS72 + P(1,13)*PS62 - P(13,13)*PS44 + P(13,14)*PS140 - P(13,15)*PS139 + P(2,13)*PS60 - P(3,13)*PS74 + P(4,13);
        const float PS142 = P(0,15)*PS72 + P(1,15)*PS62 - P(13,15)*PS44 + P(14,15)*PS140 - P(15,15)*PS139 + P(2,15)*PS60 - P(3,15)*PS74 + P(4,15);
        const float PS143 = P(1,3)*PS62;
        const float PS144 = P(0,3)*PS72;
        const float PS145 = P(2,3)*PS60 - P(3,13)*PS44 + P(3,14)*PS140 - P(3,15)*PS139 - P(3,3)*PS74 + P(3,4) + PS143 + PS144;
        const float PS146 = P(0,14)*PS72 + P(1,14)*PS62 - P(13,14)*PS44 + P(14,14)*PS140 - P(14,15)*PS139 + P(2,14)*PS60 - P(3,14)*PS74 + P(4,14);
        const float PS147 = P(0,2)*PS60;
        const float PS148 = P(0,3)*PS74;
        const float PS149 = P(0,0)*PS72 + P(0,1)*PS62 - P(0,13)*PS44 + P(0,14)*PS140 - P(0,15)*PS139 + P(0,4) + PS147 - PS148;
        const float PS150 = P(1,2)*PS62;
        const float PS151 = P(0,2)*PS72;
        const float PS152 = -P(2,13)*PS44 + P(2,14)*PS140 - P(2,15)*PS139 + P(2,2)*PS60 - P(2,3)*PS74 + P(2,4) + PS150 + PS151;
        const float PS153 = P(1,2)*PS60;
        const float PS154 = P(1,3)*PS74;
        const float PS155 = P(0,1)*PS72 + P(1,1)*PS62 - P(1,13)*PS44 + P(1,14)*PS140 - P(1,15)*PS139 + P(1,4) + PS153 - PS154;
        const float PS156 = 4*dvyVar;
        const float PS157 = 4*dvzVar;
        const float PS158 = P(0,4)*PS72 + P(1,4)*PS62 + P(2,4)*PS60 - P(3,4)*PS74 - P(4,13)*PS44 + P(4,14)*PS140 - P(4,15)*PS139 + P(4,4);
        const float PS159 = 2*PS141;
        const float PS160 = 2*PS68;
        const float PS161 = PS65*dvyVar;
        const float PS162 = 2*PS66;
        const float PS163 = PS44*dvxVar;
        const float PS164 = P(0,5)*PS72 + P(1,5)*PS62 + P(2,5)*PS60 - P(3,5)*PS74 + P(4,5) - P(5,13)*PS44 + P(5,14)*PS140 - P(5,15)*PS139;
        const float PS165 = 2*PS71;
        const float PS166 = 2*PS73;
        const float PS167 = PS70*dvzVar;
        const float PS168 = P(0,6)*PS72 + P(1,6)*PS62 + P(2,6)*PS60 - P(3,6)*PS74 + P(4,6) - P(6,13)*PS44 + P(6,14)*PS140 - P(6,15)*PS139;
        const float PS169 = P(0,14)*PS74 - P(1,14)*PS60 - P(13,14)*PS162 - P(14,14)*PS65 + P(14,15)*PS160 + P(2,14)*PS62 + P(3,14)*PS72 + P(5,14);
        const float PS170 = P(0,13)*PS74 - P(1,13)*PS60 - P(13,13)*PS162 - P(13,14)*PS65 + P(13,15)*PS160 + P(2,13)*PS62 + P(3,13)*PS72 + P(5,13);
        const float PS171 = P(0,1)*PS74;
        const float PS172 = -P(1,1)*PS60 - P(1,13)*PS162 - P(1,14)*PS65 + P(1,15)*PS160 + P(1,3)*PS72 + P(1,5) + PS150 + PS171;
        const float PS173 = P(0,15)*PS74 - P(1,15)*PS60 - P(13,15)*PS162 - P(14,15)*PS65 + P(15,15)*PS160 + P(2,15)*PS62 + P(3,15)*PS72 + P(5,15);
        const float PS174 = P(2,3)*PS62;
        const float PS175 = -P(1,3)*PS60 - P(3,13)*PS162 - P(3,14)*PS65 + P(3,15)*PS160 + P(3,3)*PS72 + P(3,5) + PS148 + PS174;
        const float PS176 = P(0,1)*PS60;
        const float PS177 = P(0,0)*PS74 - P(0,13)*PS162 - P(0,14)*PS65 + P(0,15)*PS160 + P(0,2)*PS62 + P(0,5) + PS144 - PS176;
        const float PS178 = P(2,3)*PS72;
        const float PS179 = P(0,2)*PS74 - P(2,13)*PS162 - P(2,14)*PS65 + P(2,15)*PS160 + P(2,2)*PS62 + P(2,5) - PS153 + PS178;
        const float PS180 = 4*dvxVar;
        const float PS181 = P(0,5)*PS74 - P(1,5)*PS60 + P(2,5)*PS62 + P(3,5)*PS72 - P(5,13)*PS162 - P(5,14)*PS65 + P(5,15)*PS160 + P(5,5);
        const float PS182 = P(0,6)*PS74 - P(1,6)*PS60 + P(2,6)*PS62 + P(3,6)*PS72 + P(5,6) - P(6,13)*PS162 - P(6,14)*PS65 + P(6,15)*PS160;
        const float PS183 = P(0,15)*PS60 + P(1,15)*PS74 + P(13,15)*PS166 - P(14,15)*PS165 - P(15,15)*PS70 - P(2,15)*PS72 + P(3,15)*PS62 + P(6,15);
        const float PS184 = P(0,14)*PS60 + P(1,14)*PS74 + P(13,14)*PS166 - P(14,14)*PS165 - P(14,15)*PS70 - P(2,14)*PS72 + P(3,14)*PS62 + P(6,14);
        const float PS185 = P(0,13)*PS60 + P(1,13)*PS74 + P(13,13)*PS166 - P(13,14)*PS165 - P(13,15)*PS70 - P(2,13)*PS72 + P(3,13)*PS62 + P(6,13);
        const float PS186 = P(0,6)*PS60 + P(1,6)*PS74 - P(2,6)*PS72 + P(3,6)*PS62 + P(6,13)*PS166 - P(6,14)*PS165 - P(6,15)*PS70 + P(6,6);

        SquareMatrix24f nextP;

        nextP(0,0) = PS0*PS1 - PS11*PS23 - PS12*PS26 - PS13*PS29 + PS14*PS6 + PS17*PS7 + PS2*PS3 + PS20*PS9 + PS33 + PS4*PS5;
        nextP(0,1) = -PS1*PS36 + PS11*PS33 - PS12*PS29 + PS13*PS26 - PS14*PS34 + PS17*PS9 - PS20*PS7 + PS23 + PS3*PS35 - PS35*PS5;
        nextP(1,1) = PS1*PS42 + PS11*PS80 - PS12*PS89 + PS13*PS82 + PS2*PS5 + PS3*PS4 - PS34*PS85 - PS7*PS87 + PS77*PS9 + PS92;
        nextP(0,2) = -PS1*PS37 + PS11*PS29 + PS12*PS33 - PS13*PS23 - PS14*PS9 - PS17*PS34 + PS20*PS6 + PS26 - PS3*PS38 + PS37*PS5;
        nextP(1,2) = PS1*PS40 + PS11*PS89 + PS12*PS80 - PS13*PS92 - PS3*PS40 - PS34*PS77 - PS39*PS5 + PS6*PS87 + PS82 - PS85*PS9;
        nextP(2,2) = PS0*PS5 + PS1*PS4 + PS103*PS6 + PS104*PS11 + PS106*PS12 - PS108*PS34 - PS109*PS9 - PS111*PS13 + PS113 + PS3*PS42;
        nextP(0,3) = PS1*PS39 - PS11*PS26 + PS12*PS23 + PS13*PS33 + PS14*PS7 - PS17*PS6 - PS20*PS34 + PS29 - PS3*PS39 - PS40*PS5;
        nextP(1,3) = -PS1*PS38 - PS11*PS82 + PS12*PS92 + PS13*PS80 - PS3*PS37 - PS34*PS87 + PS38*PS5 - PS6*PS77 + PS7*PS85 + PS89;
        nextP(2,3) = -PS1*PS35 - PS103*PS34 + PS104 + PS106*PS13 - PS108*PS6 + PS109*PS7 - PS11*PS113 + PS111*PS12 + PS3*PS36 - PS36*PS5;
        nextP(3,3) = PS0*PS3 + PS1*PS2 - PS11*PS128 + PS12*PS124 + PS123*PS7 + PS125*PS13 - PS126*PS34 - PS127*PS6 + PS129 + PS42*PS5;
        nextP(0,4) = PS23*PS62 + PS26*PS60 - PS44*PS45 - PS46*PS48 - PS52*PS53 + PS54*PS56 + PS57*PS58 + PS63;
        nextP(1,4) = -PS44*PS93 - PS46*PS95 + PS54*PS97 + PS60*PS82 + PS62*PS92 + PS72*PS80 - PS74*PS89 + PS98;
        nextP(2,4) = -PS104*PS74 + PS106*PS72 + PS111*PS62 + PS113*PS60 - PS114*PS44 - PS116*PS46 + PS118*PS54 + PS119;
        nextP(3,4) = PS124*PS62 + PS125*PS72 + PS128*PS60 - PS129*PS74 - PS130*PS44 - PS132*PS46 + PS134*PS54 + PS135;
        nextP(4,4) = -PS139*PS142 + PS140*PS146 - PS141*PS44 - PS145*PS74 + PS149*PS72 + PS152*PS60 + PS155*PS62 + PS156*ecl::powf(PS54, 2) + PS157*ecl::powf(PS46, 2) + PS158 + ecl::powf(PS44, 2)*dvxVar;
        nextP(0,5) = -PS23*PS60 + PS26*PS62 + PS48*PS68 + PS52*PS58 + PS53*PS57 - PS55*PS65 - PS66*PS67 + PS69;
        nextP(1,5) = PS100 - PS60*PS92 + PS62*PS82 - PS65*PS96 - PS66*PS99 + PS68*PS95 + PS72*PS89 + PS74*PS80;
        nextP(2,5) = PS104*PS72 + PS106*PS74 - PS111*PS60 + PS113*PS62 + PS116*PS68 - PS117*PS65 - PS120*PS66 + PS121;
        nextP(3,5) = -PS124*PS60 + PS125*PS74 + PS128*PS62 + PS129*PS72 + PS132*PS68 - PS133*PS65 - PS136*PS66 + PS137;
        nextP(4,5) = -PS140*PS161 + PS142*PS160 + PS145*PS72 - PS146*PS65 + PS149*PS74 + PS152*PS62 - PS155*PS60 - PS157*PS46*PS68 - PS159*PS66 + PS162*PS163 + PS164;
        nextP(5,5) = PS157*ecl::powf(PS68, 2) + PS160*PS173 - PS162*PS170 - PS169*PS65 - PS172*PS60 + PS175*PS72 + PS177*PS74 + PS179*PS62 + PS180*ecl::powf(PS66, 2) + PS181 + ecl::powf(PS65, 2)*dvyVar;
        nextP(0,6) = PS23*PS74 - PS26*PS72 - PS47*PS70 + PS53*PS61 - PS56*PS71 + PS58*PS59 + PS67*PS73 + PS75;
        nextP(1,6) = PS101 + PS60*PS80 + PS62*PS89 - PS70*PS94 - PS71*PS97 - PS72*PS82 + PS73*PS99 + PS74*PS92;
        nextP(2,6) = PS104*PS62 + PS106*PS60 + PS111*PS74 - PS113*PS72 - PS115*PS70 - PS118*PS71 + PS120*PS73 + PS122;
        nextP(3,6) = PS124*PS74 + PS125*PS60 - PS128*PS72 + PS129*PS62 - PS131*PS70 - PS134*PS71 + PS136*PS73 + PS138;
        nextP(4,6) = PS139*PS167 - PS142*PS70 + PS145*PS62 - PS146*PS165 + PS149*PS60 - PS152*PS72 + PS155*PS74 - PS156*PS54*PS71 + PS159*PS73 - PS163*PS166 + PS168;
        nextP(5,6) = -PS160*PS167 + PS161*PS165 - PS165*PS169 + PS166*PS170 + PS172*PS74 - PS173*PS70 + PS175*PS62 + PS177*PS60 - PS179*PS72 - PS180*PS66*PS73 + PS182;
        nextP(6,6) = PS156*ecl::powf(PS71, 2) - PS165*PS184 + PS166*PS185 + PS180*ecl::powf(PS73, 2) - PS183*PS70 + PS186 + PS60*(P(0,0)*PS60 + P(0,13)*PS166 - P(0,14)*PS165 - P(0,15)*PS70 + P(0,3)*PS62 + P(0,6) - PS151 + PS171) + PS62*(P(0,3)*PS60 + P(3,13)*PS166 - P(3,14)*PS165 - P(3,15)*PS70 + P(3,3)*PS62 + P(3,6) + PS154 - PS178) + ecl::powf(PS70, 2)*dvzVar - PS72*(P(1,2)*PS74 + P(2,13)*PS166 - P(2,14)*PS165 - P(2,15)*PS70 - P(2,2)*PS72 + P(2,6) + PS147 + PS174) + PS74*(P(1,1)*PS74 + P(1,13)*PS166 - P(1,14)*PS165 - P(1,15)*PS70 - P(1,2)*PS72 + P(1,6) + PS143 + PS176);
        nextP(0,7) = P(0,7) - P(1,7)*PS11 - P(2,7)*PS12 - P(3,7)*PS13 + P(7,10)*PS6 + P(7,11)*PS7 + P(7,12)*PS9 + PS63*dt;
        nextP(1,7) = P(0,7)*PS11 + P(1,7) + P(2,7)*PS13 - P(3,7)*PS12 - P(7,10)*PS34 + P(7,11)*PS9 - P(7,12)*PS7 + PS98*dt;
        nextP(2,7) = P(0,7)*PS12 - P(1,7)*PS13 + P(2,7) + P(3,7)*PS11 - P(7,10)*PS9 - P(7,11)*PS34 + P(7,12)*PS6 + PS119*dt;
        nextP(3,7) = P(0,7)*PS13 + P(1,7)*PS12 - P(2,7)*PS11 + P(3,7) + P(7,10)*PS7 - P(7,11)*PS6 - P(7,12)*PS34 + PS135*dt;
        nextP(4,7) = P(0,7)*PS72 + P(1,7)*PS62 + P(2,7)*PS60 - P(3,7)*PS74 + P(4,7) - P(7,13)*PS44 + P(7,14)*PS140 - P(7,15)*PS139 + PS158*dt;
        nextP(5,7) = P(0,7)*PS74 - P(1,7)*PS60 + P(2,7)*PS62 + P(3,7)*PS72 + P(5,7) - P(7,13)*PS162 - P(7,14)*PS65 + P(7,15)*PS160 + dt*(P(0,4)*PS74 - P(1,4)*PS60 + P(2,4)*PS62 + P(3,4)*PS72 - P(4,13)*PS162 - P(4,14)*PS65 + P(4,15)*PS160 + P(4,5));
        nextP(6,7) = P(0,7)*PS60 + P(1,7)*PS74 - P(2,7)*PS72 + P(3,7)*PS62 + P(6,7) + P(7,13)*PS166 - P(7,14)*PS165 - P(7,15)*PS70 + dt*(P(0,4)*PS60 + P(1,4)*PS74 - P(2,4)*PS72 + P(3,4)*PS62 + P(4,13)*PS166 - P(4,14)*PS165 - P(4,15)*PS70 + P(4,6));
        nextP(7,7) = P(4,7)*dt + P(7,7) + dt*(P(4,4)*dt + P(4,7));
        nextP(0,8) = P(0,8) - P(1,8)*PS11 - P(2,8)*PS12 - P(3,8)*PS13 + P(8,10)*PS6 + P(8,11)*PS7 + P(8,12)*PS9 + PS69*dt;
        nextP(1,8) = P(0,8)*PS11 + P(1,8) + P(2,8)*PS13 - P(3,8)*PS12 - P(8,10)*PS34 + P(8,11)*PS9 - P(8,12)*PS7 + PS100*dt;
        nextP(2,8) = P(0,8)*PS12 - P(1,8)*PS13 + P(2,8) + P(3,8)*PS11 - P(8,10)*PS9 - P(8,11)*PS34 + P(8,12)*PS6 + PS121*dt;
        nextP(3,8) = P(0,8)*PS13 + P(1,8)*PS12 - P(2,8)*PS11 + P(3,8) + P(8,10)*PS7 - P(8,11)*PS6 - P(8,12)*PS34 + PS137*dt;
        nextP(4,8) = P(0,8)*PS72 + P(1,8)*PS62 + P(2,8)*PS60 - P(3,8)*PS74 + P(4,8) - P(8,13)*PS44 + P(8,14)*PS140 - P(8,15)*PS139 + PS164*dt;
        nextP(5,8) = P(0,8)*PS74 - P(1,8)*PS60 + P(2,8)*PS62 + P(3,8)*PS72 + P(5,8) - P(8,13)*PS162 - P(8,14)*PS65 + P(8,15)*PS160 + PS181*dt;
        nextP(6,8) = P(0,8)*PS60 + P(1,8)*PS74 - P(2,8)*PS72 + P(3,8)*PS62 + P(6,8) + P(8,13)*PS166 - P(8,14)*PS165 - P(8,15)*PS70 + dt*(P(0,5)*PS60 + P(1,5)*PS74 - P(2,5)*PS72 + P(3,5)*PS62 + P(5,13)*PS166 - P(5,14)*PS165 - P(5,15)*PS70 + P(5,6));
        nextP(7,8) = P(4,8)*dt + P(7,8) + dt*(P(4,5)*dt + P(5,7));
        nextP(8,8) = P(5,8)*dt + P(8,8) + dt*(P(5,5)*dt + P(5,8));
        nextP(0,9) = P(0,9) - P(1,9)*PS11 - P(2,9)*PS12 - P(3,9)*PS13 + P(9,10)*PS6 + P(9,11)*PS7 + P(9,12)*PS9 + PS75*dt;
        nextP(1,9) = P(0,9)*PS11 + P(1,9) + P(2,9)*PS13 - P(3,9)*PS12 - P(9,10)*PS34 + P(9,11)*PS9 - P(9,12)*PS7 + PS101*dt;
        nextP(2,9) = P(0,9)*PS12 - P(1,9)*PS13 + P(2,9) + P(3,9)*PS11 - P(9,10)*PS9 - P(9,11)*PS34 + P(9,12)*PS6 + PS122*dt;
        nextP(3,9) = P(0,9)*PS13 + P(1,9)*PS12 - P(2,9)*PS11 + P(3,9) + P(9,10)*PS7 - P(9,11)*PS6 - P(9,12)*PS34 + PS138*dt;
        nextP(4,9) = P(0,9)*PS72 + P(1,9)*PS62 + P(2,9)*PS60 - P(3,9)*PS74 + P(4,9) - P(9,13)*PS44 + P(9,14)*PS140 - P(9,15)*PS139 + PS168*dt;
        nextP(5,9) = P(0,9)*PS74 - P(1,9)*PS60 + P(2,9)*PS62 + P(3,9)*PS72 + P(5,9) - P(9,13)*PS162 - P(9,14)*PS65 + P(9,15)*PS160 + PS182*dt;
        nextP(6,9) = P(0,9)*PS60 + P(1,9)*PS74 - P(2,9)*PS72 + P(3,9)*PS62 + P(6,9) + P(9,13)*PS166 - P(9,14)*PS165 - P(9,15)*PS70 + PS186*dt;
        nextP(7,9) = P(4,9)*dt + P(7,9) + dt*(P(4,6)*dt + P(6,7));
        nextP(8,9) = P(5,9)*dt + P(8,9) + dt*(P(5,6)*dt + P(6,8));
        nextP(9,9) = P(6,9)*dt + P(9,9) + dt*(P(6,6)*dt + P(6,9));
        nextP(0,10) = PS14;
        nextP(1,10) = PS85;
        nextP(2,10) = PS109;
        nextP(3,10) = PS123;
        nextP(4,10) = P(0,10)*PS72 + P(1,10)*PS62 - P(10,13)*PS44 + P(10,14)*PS140 - P(10,15)*PS139 + P(2,10)*PS60 - P(3,10)*PS74 + P(4,10);
        nextP(5,10) = P(0,10)*PS74 - P(1,10)*PS60 - P(10,13)*PS162 - P(10,14)*PS65 + P(10,15)*PS160 + P(2,10)*PS62 + P(3,10)*PS72 + P(5,10);
        nextP(6,10) = P(0,10)*PS60 + P(1,10)*PS74 + P(10,13)*PS166 - P(10,14)*PS165 - P(10,15)*PS70 - P(2,10)*PS72 + P(3,10)*PS62 + P(6,10);
        nextP(7,10) = P(4,10)*dt + P(7,10);
        nextP(8,10) = P(5,10)*dt + P(8,10);
        nextP(9,10) = P(6,10)*dt + P(9,10);
        nextP(10,10) = P(10,10);
        nextP(0,11) = PS17;
        nextP(1,11) = PS77;
        nextP(2,11) = PS108;
        nextP(3,11) = PS127;
        nextP(4,11) = P(0,11)*PS72 + P(1,11)*PS62 - P(11,13)*PS44 + P(11,14)*PS140 - P(11,15)*PS139 + P(2,11)*PS60 - P(3,11)*PS74 + P(4,11);
        nextP(5,11) = P(0,11)*PS74 - P(1,11)*PS60 - P(11,13)*PS162 - P(11,14)*PS65 + P(11,15)*PS160 + P(2,11)*PS62 + P(3,11)*PS72 + P(5,11);
        nextP(6,11) = P(0,11)*PS60 + P(1,11)*PS74 + P(11,13)*PS166 - P(11,14)*PS165 - P(11,15)*PS70 - P(2,11)*PS72 + P(3,11)*PS62 + P(6,11);
        nextP(7,11) = P(4,11)*dt + P(7,11);
        nextP(8,11) = P(5,11)*dt + P(8,11);
        nextP(9,11) = P(6,11)*dt + P(9,11);
        nextP(10,11) = P(10,11);
        nextP(11,11) = P(11,11);
        nextP(0,12) = PS20;
        nextP(1,12) = PS87;
        nextP(2,12) = PS103;
        nextP(3,12) = PS126;
        nextP(4,12) = P(0,12)*PS72 + P(1,12)*PS62 - P(12,13)*PS44 + P(12,14)*PS140 - P(12,15)*PS139 + P(2,12)*PS60 - P(3,12)*PS74 + P(4,12);
        nextP(5,12) = P(0,12)*PS74 - P(1,12)*PS60 - P(12,13)*PS162 - P(12,14)*PS65 + P(12,15)*PS160 + P(2,12)*PS62 + P(3,12)*PS72 + P(5,12);
        nextP(6,12) = P(0,12)*PS60 + P(1,12)*PS74 + P(12,13)*PS166 - P(12,14)*PS165 - P(12,15)*PS70 - P(2,12)*PS72 + P(3,12)*PS62 + P(6,12);
        nextP(7,12) = P(4,12)*dt + P(7,12);
        nextP(8,12) = P(5,12)*dt + P(8,12);
        nextP(9,12) = P(6,12)*dt + P(9,12);
        nextP(10,12) = P(10,12);
        nextP(11,12) = P(11,12);
        nextP(12,12) = P(12,12);
        nextP(0,13) = PS45;
        nextP(1,13) = PS93;
        nextP(2,13) = PS114;
        nextP(3,13) = PS130;
        nextP(4,13) = PS141;
        nextP(5,13) = PS170;
        nextP(6,13) = PS185;
        nextP(7,13) = P(4,13)*dt + P(7,13);
        nextP(8,13) = P(5,13)*dt + P(8,13);
        nextP(9,13) = P(6,13)*dt + P(9,13);
        nextP(10,13) = P(10,13);
        nextP(11,13) = P(11,13);
        nextP(12,13) = P(12,13);
        nextP(13,13) = P(13,13);
        nextP(0,14) = PS55;
        nextP(1,14) = PS96;
        nextP(2,14) = PS117;
        nextP(3,14) = PS133;
        nextP(4,14) = PS146;
        nextP(5,14) = PS169;
        nextP(6,14) = PS184;
        nextP(7,14) = P(4,14)*dt + P(7,14);
        nextP(8,14) = P(5,14)*dt + P(8,14);
        nextP(9,14) = P(6,14)*dt + P(9,14);
        nextP(10,14) = P(10,14);
        nextP(11,14) = P(11,14);
        nextP(12,14) = P(12,14);
        nextP(13,14) = P(13,14);
        nextP(14,14) = P(14,14);
        nextP(0,15) = PS47;
        nextP(1,15) = PS94;
        nextP(2,15) = PS115;
        nextP(3,15) = PS131;
        nextP(4,15) = PS142;
        nextP(5,15) = PS173;
        nextP(6,15) = PS183;
        nextP(7,15) = P(4,15)*dt + P(7,15);
        nextP(8,15) = P(5,15)*dt + P(8,15);
        nextP(9,15) = P(6,15)*dt + P(9,15);
        nextP(10,15) = P(10,15);
        nextP(11,15) = P(11,15);
        nextP(12,15) = P(12,15);
        nextP(13,15) = P(13,15);
        nextP(14,15) = P(14,15);
        nextP(15,15) = P(15,15);
        nextP(0,16) = P(0,16) - P(1,16)*PS11 + P(10,16)*PS6 + P(11,16)*PS7 + P(12,16)*PS9 - P(2,16)*PS12 - P(3,16)*PS13;
        nextP(1,16) = P(0,16)*PS11 + P(1,16) - P(10,16)*PS34 + P(11,16)*PS9 - P(12,16)*PS7 + P(2,16)*PS13 - P(3,16)*PS12;
        nextP(2,16) = P(0,16)*PS12 - P(1,16)*PS13 - P(10,16)*PS9 - P(11,16)*PS34 + P(12,16)*PS6 + P(2,16) + P(3,16)*PS11;
        nextP(3,16) = P(0,16)*PS13 + P(1,16)*PS12 + P(10,16)*PS7 - P(11,16)*PS6 - P(12,16)*PS34 - P(2,16)*PS11 + P(3,16);
        nextP(4,16) = P(0,16)*PS72 + P(1,16)*PS62 - P(13,16)*PS44 + P(14,16)*PS140 - P(15,16)*PS139 + P(2,16)*PS60 - P(3,16)*PS74 + P(4,16);
        nextP(5,16) = P(0,16)*PS74 - P(1,16)*PS60 - P(13,16)*PS162 - P(14,16)*PS65 + P(15,16)*PS160 + P(2,16)*PS62 + P(3,16)*PS72 + P(5,16);
        nextP(6,16) = P(0,16)*PS60 + P(1,16)*PS74 + P(13,16)*PS166 - P(14,16)*PS165 - P(15,16)*PS70 - P(2,16)*PS72 + P(3,16)*PS62 + P(6,16);
        nextP(7,16) = P(4,16)*dt + P(7,16);
        nextP(8,16) = P(5,16)*dt + P(8,16);
        nextP(9,16) = P(6,16)*dt + P(9,16);
        nextP(10,16) = P(10,16);
        nextP(11,16) = P(11,16);
        nextP(12,16) = P(12,16);
        nextP(13,16) = P(13,16);
        nextP(14,16) = P(14,16);
        nextP(15,16) = P(15,16);
        nextP(16,16) = P(16,16);
        nextP(0,17) = P(0,17) - P(1,17)*PS11 + P(10,17)*PS6 + P(11,17)*PS7 + P(12,17)*PS9 - P(2,17)*PS12 - P(3,17)*PS13;
        nextP(1,17) = P(0,17)*PS11 + P(1,17) - P(10,17)*PS34 + P(11,17)*PS9 - P(12,17)*PS7 + P(2,17)*PS13 - P(3,17)*PS12;
        nextP(2,17) = P(0,17)*PS12 - P(1,17)*PS13 - P(10,17)*PS9 - P(11,17)*PS34 + P(12,17)*PS6 + P(2,17) + P(3,17)*PS11;
        nextP(3,17) = P(0,17)*PS13 + P(1,17)*PS12 + P(10,17)*PS7 - P(11,17)*PS6 - P(12,17)*PS34 - P(2,17)*PS11 + P(3,17);
        nextP(4,17) = P(0,17)*PS72 + P(1,17)*PS62 - P(13,17)*PS44 + P(14,17)*PS140 - P(15,17)*PS139 + P(2,17)*PS60 - P(3,17)*PS74 + P(4,17);
        nextP(5,17) = P(0,17)*PS74 - P(1,17)*PS60 - P(13,17)*PS162 - P(14,17)*PS65 + P(15,17)*PS160 + P(2,17)*PS62 + P(3,17)*PS72 + P(5,17);
        nextP(6,17) = P(0,17)*PS60 + P(1,17)*PS74 + P(13,17)*PS166 - P(14,17)*PS165 - P(15,17)*PS70 - P(2,17)*PS72 + P(3,17)*PS62 + P(6,17);
        nextP(7,17) = P(4,17)*dt + P(7,17);
        nextP(8,17) = P(5,17)*dt + P(8,17);
        nextP(9,17) = P(6,17)*dt + P(9,17);
        nextP(10,17) = P(10,17);
        nextP(11,17) = P(11,17);
        nextP(12,17) = P(12,17);
        nextP(13,17) = P(13,17);
        nextP(14,17) = P(14,17);
        nextP(15,17) = P(15,17);
        nextP(16,17) = P(16,17);
        nextP(17,17) = P(17,17);
        nextP(0,18) = P(0,18) - P(1,18)*PS11 + P(10,18)*PS6 + P(11,18)*PS7 + P(12,18)*PS9 - P(2,18)*PS12 - P(3,18)*PS13;
        nextP(1,18) = P(0,18)*PS11 + P(1,18) - P(10,18)*PS34 + P(11,18)*PS9 - P(12,18)*PS7 + P(2,18)*PS13 - P(3,18)*PS12;
        nextP(2,18) = P(0,18)*PS12 - P(1,18)*PS13 - P(10,18)*PS9 - P(11,18)*PS34 + P(12,18)*PS6 + P(2,18) + P(3,18)*PS11;
        nextP(3,18) = P(0,18)*PS13 + P(1,18)*PS12 + P(10,18)*PS7 - P(11,18)*PS6 - P(12,18)*PS34 - P(2,18)*PS11 + P(3,18);
        nextP(4,18) = P(0,18)*PS72 + P(1,18)*PS62 - P(13,18)*PS44 + P(14,18)*PS140 - P(15,18)*PS139 + P(2,18)*PS60 - P(3,18)*PS74 + P(4,18);
        nextP(5,18) = P(0,18)*PS74 - P(1,18)*PS60 - P(13,18)*PS162 - P(14,18)*PS65 + P(15,18)*PS160 + P(2,18)*PS62 + P(3,18)*PS72 + P(5,18);
        nextP(6,18) = P(0,18)*PS60 + P(1,18)*PS74 + P(13,18)*PS166 - P(14,18)*PS165 - P(15,18)*PS70 - P(2,18)*PS72 + P(3,18)*PS62 + P(6,18);
        nextP(7,18) = P(4,18)*dt + P(7,18);
        nextP(8,18) = P(5,18)*dt + P(8,18);
        nextP(9,18) = P(6,18)*dt + P(9,18);
        nextP(10,18) = P(10,18);
        nextP(11,18) = P(11,18);
        nextP(12,18) = P(12,18);
        nextP(13,18) = P(13,18);
        nextP(14,18) = P(14,18);
        nextP(15,18) = P(15,18);
        nextP(16,18) = P(16,18);
        nextP(17,18) = P(17,18);
        nextP(18,18) = P(18,18);
        nextP(0,19) = P(0,19) - P(1,19)*PS11 + P(10,19)*PS6 + P(11,19)*PS7 + P(12,19)*PS9 - P(2,19)*PS12 - P(3,19)*PS13;
        nextP(1,19) = P(0,19)*PS11 + P(1,19) - P(10,19)*PS34 + P(11,19)*PS9 - P(12,19)*PS7 + P(2,19)*PS13 - P(3,19)*PS12;
        nextP(2,19) = P(0,19)*PS12 - P(1,19)*PS13 - P(10,19)*PS9 - P(11,19)*PS34 + P(12,19)*PS6 + P(2,19) + P(3,19)*PS11;
        nextP(3,19) = P(0,19)*PS13 + P(1,19)*PS12 + P(10,19)*PS7 - P(11,19)*PS6 - P(12,19)*PS34 - P(2,19)*PS11 + P(3,19);
        nextP(4,19) = P(0,19)*PS72 + P(1,19)*PS62 - P(13,19)*PS44 + P(14,19)*PS140 - P(15,19)*PS139 + P(2,19)*PS60 - P(3,19)*PS74 + P(4,19);
        nextP(5,19) = P(0,19)*PS74 - P(1,19)*PS60 - P(13,19)*PS162 - P(14,19)*PS65 + P(15,19)*PS160 + P(2,19)*PS62 + P(3,19)*PS72 + P(5,19);
        nextP(6,19) = P(0,19)*PS60 + P(1,19)*PS74 + P(13,19)*PS166 - P(14,19)*PS165 - P(15,19)*PS70 - P(2,19)*PS72 + P(3,19)*PS62 + P(6,19);
        nextP(7,19) = P(4,19)*dt + P(7,19);
        nextP(8,19) = P(5,19)*dt + P(8,19);
        nextP(9,19) = P(6,19)*dt + P(9,19);
        nextP(10,19) = P(10,19);
        nextP(11,19) = P(11,19);
        nextP(12,19) = P(12,19);
        nextP(13,19) = P(13,19);
        nextP(14,19) = P(14,19);
        nextP(15,19) = P(15,19);
        nextP(16,19) = P(16,19);
        nextP(17,19) = P(17,19);
        nextP(18,19) = P(18,19);
        nextP(19,19) = P(19,19);
        nextP(0,20) = P(0,20) - P(1,20)*PS11 + P(10,20)*PS6 + P(11,20)*PS7 + P(12,20)*PS9 - P(2,20)*PS12 - P(3,20)*PS13;
        nextP(1,20) = P(0,20)*PS11 + P(1,20) - P(10,20)*PS34 + P(11,20)*PS9 - P(12,20)*PS7 + P(2,20)*PS13 - P(3,20)*PS12;
        nextP(2,20) = P(0,20)*PS12 - P(1,20)*PS13 - P(10,20)*PS9 - P(11,20)*PS34 + P(12,20)*PS6 + P(2,20) + P(3,20)*PS11;
        nextP(3,20) = P(0,20)*PS13 + P(1,20)*PS12 + P(10,20)*PS7 - P(11,20)*PS6 - P(12,20)*PS34 - P(2,20)*PS11 + P(3,20);
        nextP(4,20) = P(0,20)*PS72 + P(1,20)*PS62 - P(13,20)*PS44 + P(14,20)*PS140 - P(15,20)*PS139 + P(2,20)*PS60 - P(3,20)*PS74 + P(4,20);
        nextP(5,20) = P(0,20)*PS74 - P(1,20)*PS60 - P(13,20)*PS162 - P(14,20)*PS65 + P(15,20)*PS160 + P(2,20)*PS62 + P(3,20)*PS72 + P(5,20);
        nextP(6,20) = P(0,20)*PS60 + P(1,20)*PS74 + P(13,20)*PS166 - P(14,20)*PS165 - P(15,20)*PS70 - P(2,20)*PS72 + P(3,20)*PS62 + P(6,20);
        nextP(7,20) = P(4,20)*dt + P(7,20);
        nextP(8,20) = P(5,20)*dt + P(8,20);
        nextP(9,20) = P(6,20)*dt + P(9,20);
        nextP(10,20) = P(10,20);
        nextP(11,20) = P(11,20);
        nextP(12,20) = P(12,20);
        nextP(13,20) = P(13,20);
        nextP(14,20) = P(14,20);
        nextP(15,20) = P(15,20);
        nextP(16,20) = P(16,20);
        nextP(17,20) = P(17,20);
        nextP(18,20) = P(18,20);
        nextP(19,20) = P(19,20);
        nextP(20,20) = P(20,20);
        nextP(0,21) = P(0,21) - P(1,21)*PS11 + P(10,21)*PS6 + P(11,21)*PS7 + P(12,21)*PS9 - P(2,21)*PS12 - P(3,21)*PS13;
        nextP(1,21) = P(0,21)*PS11 + P(1,21) - P(10,21)*PS34 + P(11,21)*PS9 - P(12,21)*PS7 + P(2,21)*PS13 - P(3,21)*PS12;
        nextP(2,21) = P(0,21)*PS12 - P(1,21)*PS13 - P(10,21)*PS9 - P(11,21)*PS34 + P(12,21)*PS6 + P(2,21) + P(3,21)*PS11;
        nextP(3,21) = P(0,21)*PS13 + P(1,21)*PS12 + P(10,21)*PS7 - P(11,21)*PS6 - P(12,21)*PS34 - P(2,21)*PS11 + P(3,21);
        nextP(4,21) = P(0,21)*PS72 + P(1,21)*PS62 - P(13,21)*PS44 + P(14,21)*PS140 - P(15,21)*PS139 + P(2,21)*PS60 - P(3,21)*PS74 + P(4,21);
        nextP(5,21) = P(0,21)*PS74 - P(1,21)*PS60 - P(13,21)*PS162 - P(14,21)*PS65 + P(15,21)*PS160 + P(2,21)*PS62 + P(3,21)*PS72 + P(5,21);
        nextP(6,21) = P(0,21)*PS60 + P(1,21)*PS74 + P(13,21)*PS166 - P(14,21)*PS165 - P(15,21)*PS70 - P(2,21)*PS72 + P(3,21)*PS62 + P(6,21);
        nextP(7,21) = P(4,21)*dt + P(7,21);
        nextP(8,21) = P(5,21)*dt + P(8,21);
        nextP(9,21) = P(6,21)*dt + P(9,21);
        nextP(10,21) = P(10,21);
        nextP(11,21) = P(11,21);
        nextP(12,21) = P(12,21);
        nextP(13,21) = P(13,21);
        nextP(14,21) = P(14,21);
        nextP(15,21) = P(15,21);
        nextP(16,21) = P(16,21);
        nextP(17,21) = P(17,21);
        nextP(18,21) = P(18,21);
        nextP(19,21) = P(19,21);
        nextP(20,21) = P(20,21);
        nextP(21,21) = P(21,21);
        nextP(0,22) = P(0,22) - P(1,22)*PS11 + P(10,22)*PS6 + P(11,22)*PS7 + P(12,22)*PS9 - P(2,22)*PS12 - P(3,22)*PS13;
        nextP(1,22) = P(0,22)*PS11 + P(1,22) - P(10,22)*PS34 + P(11,22)*PS9 - P(12,22)*PS7 + P(2,22)*PS13 - P(3,22)*PS12;
        nextP(2,22) = P(0,22)*PS12 - P(1,22)*PS13 - P(10,22)*PS9 - P(11,22)*PS34 + P(12,22)*PS6 + P(2,22) + P(3,22)*PS11;
        nextP(3,22) = P(0,22)*PS13 + P(1,22)*PS12 + P(10,22)*PS7 - P(11,22)*PS6 - P(12,22)*PS34 - P(2,22)*PS11 + P(3,22);
        nextP(4,22) = P(0,22)*PS72 + P(1,22)*PS62 - P(13,22)*PS44 + P(14,22)*PS140 - P(15,22)*PS139 + P(2,22)*PS60 - P(3,22)*PS74 + P(4,22);
        nextP(5,22) = P(0,22)*PS74 - P(1,22)*PS60 - P(13,22)*PS162 - P(14,22)*PS65 + P(15,22)*PS160 + P(2,22)*PS62 + P(3,22)*PS72 + P(5,22);
        nextP(6,22) = P(0,22)*PS60 + P(1,22)*PS74 + P(13,22)*PS166 - P(14,22)*PS165 - P(15,22)*PS70 - P(2,22)*PS72 + P(3,22)*PS62 + P(6,22);
        nextP(7,22) = P(4,22)*dt + P(7,22);
        nextP(8,22) = P(5,22)*dt + P(8,22);
        nextP(9,22) = P(6,22)*dt + P(9,22);
        nextP(10,22) = P(10,22);
        nextP(11,22) = P(11,22);
        nextP(12,22) = P(12,22);
        nextP(13,22) = P(13,22);
        nextP(14,22) = P(14,22);
        nextP(15,22) = P(15,22);
        nextP(16,22) = P(16,22);
        nextP(17,22) = P(17,22);
        nextP(18,22) = P(18,22);
        nextP(19,22) = P(19,22);
        nextP(20,22) = P(20,22);
        nextP(21,22) = P(21,22);
        nextP(22,22) = P(22,22);
        nextP(0,23) = P(0,23) - P(1,23)*PS11 + P(10,23)*PS6 + P(11,23)*PS7 + P(12,23)*PS9 - P(2,23)*PS12 - P(3,23)*PS13;
        nextP(1,23) = P(0,23)*PS11 + P(1,23) - P(10,23)*PS34 + P(11,23)*PS9 - P(12,23)*PS7 + P(2,23)*PS13 - P(3,23)*PS12;
        nextP(2,23) = P(0,23)*PS12 - P(1,23)*PS13 - P(10,23)*PS9 - P(11,23)*PS34 + P(12,23)*PS6 + P(2,23) + P(3,23)*PS11;
        nextP(3,23) = P(0,23)*PS13 + P(1,23)*PS12 + P(10,23)*PS7 - P(11,23)*PS6 - P(12,23)*PS34 - P(2,23)*PS11 + P(3,23);
        nextP(4,23) = P(0,23)*PS72 + P(1,23)*PS62 - P(13,23)*PS44 + P(14,23)*PS140 - P(15,23)*PS139 + P(2,23)*PS60 - P(3,23)*PS74 + P(4,23);
        nextP(5,23) = P(0,23)*PS74 - P(1,23)*PS60 - P(13,23)*PS162 - P(14,23)*PS65 + P(15,23)*PS160 + P(2,23)*PS62 + P(3,23)*PS72 + P(5,23);
        nextP(6,23) = P(0,23)*PS60 + P(1,23)*PS74 + P(13,23)*PS166 - P(14,23)*PS165 - P(15,23)*PS70 - P(2,23)*PS72 + P(3,23)*PS62 + P(6,23);
        nextP(7,23) = P(4,23)*dt + P(7,23);
        nextP(8,23) = P(5,23)*dt + P(8,23);
        nextP(9,23) = P(6,23)*dt + P(9,23);
        nextP(10,23) = P(10,23);
        nextP(11,23) = P(11,23);
        nextP(12,23) = P(12,23);
        nextP(13,23) = P(13,23);
        nextP(14,23) = P(14,23);
        nextP(15,23) = P(15,23);
        nextP(16,23) = P(16,23);
        nextP(17,23) = P(17,23);
        nextP(18,23) = P(18,23);
        nextP(19,23) = P(19,23);
        nextP(20,23) = P(20,23);
        nextP(21,23) = P(21,23);
        nextP(22,23) = P(22,23);
        nextP(23,23) = P(23,23);

        // save output and repeat calculation using legacy matlab generated code
        SquareMatrix24f nextP_sympy;
        for (int col=0; col<=23; col++) {
            for (int row=0; row<=col; row++) {
                nextP_sympy(row,col) = nextP(row,col);
            }
        }

        // intermediate calculations
        float SF[21];
        SF[0] = dvz - dvz_b;
        SF[1] = dvy - dvy_b;
        SF[2] = dvx - dvx_b;
        SF[3] = 2*q1*SF[2] + 2*q2*SF[1] + 2*q3*SF[0];
        SF[4] = 2*q0*SF[1] - 2*q1*SF[0] + 2*q3*SF[2];
        SF[5] = 2*q0*SF[2] + 2*q2*SF[0] - 2*q3*SF[1];
        SF[6] = day*0.5f - day_b*0.5f;
        SF[7] = daz*0.5f - daz_b*0.5f;
        SF[8] = dax*0.5f - dax_b*0.5f;
        SF[9] = dax_b*0.5f - dax*0.5f;
        SF[10] = daz_b*0.5f - daz*0.5f;
        SF[11] = day_b*0.5f - day*0.5f;
        SF[12] = 2*q1*SF[1];
        SF[13] = 2*q0*SF[0];
        SF[14] = q1*0.5f;
        SF[15] = q2*0.5f;
        SF[16] = q3*0.5f;
        SF[17] = sq(q3);
        SF[18] = sq(q2);
        SF[19] = sq(q1);
        SF[20] = sq(q0);

        float SG[8];
        SG[0] = q0*0.5f;
        SG[1] = sq(q3);
        SG[2] = sq(q2);
        SG[3] = sq(q1);
        SG[4] = sq(q0);
        SG[5] = 2*q2*q3;
        SG[6] = 2*q1*q3;
        SG[7] = 2*q1*q2;

        float SQ[11];
        SQ[0] = dvzVar*(SG[5] - 2*q0*q1)*(SG[1] - SG[2] - SG[3] + SG[4]) - dvyVar*(SG[5] + 2*q0*q1)*(SG[1] - SG[2] + SG[3] - SG[4]) + dvxVar*(SG[6] - 2*q0*q2)*(SG[7] + 2*q0*q3);
        SQ[1] = dvzVar*(SG[6] + 2*q0*q2)*(SG[1] - SG[2] - SG[3] + SG[4]) - dvxVar*(SG[6] - 2*q0*q2)*(SG[1] + SG[2] - SG[3] - SG[4]) + dvyVar*(SG[5] + 2*q0*q1)*(SG[7] - 2*q0*q3);
        SQ[2] = dvzVar*(SG[5] - 2*q0*q1)*(SG[6] + 2*q0*q2) - dvyVar*(SG[7] - 2*q0*q3)*(SG[1] - SG[2] + SG[3] - SG[4]) - dvxVar*(SG[7] + 2*q0*q3)*(SG[1] + SG[2] - SG[3] - SG[4]);
        SQ[3] = (dayVar*q1*SG[0])*0.5f - (dazVar*q1*SG[0])*0.5f - (daxVar*q2*q3)*0.25f;
        SQ[4] = (dazVar*q2*SG[0])*0.5f - (daxVar*q2*SG[0])*0.5f - (dayVar*q1*q3)*0.25f;
        SQ[5] = (daxVar*q3*SG[0])*0.5f - (dayVar*q3*SG[0])*0.5f - (dazVar*q1*q2)*0.25f;
        SQ[6] = (daxVar*q1*q2)*0.25f - (dazVar*q3*SG[0])*0.5f - (dayVar*q1*q2)*0.25f;
        SQ[7] = (dazVar*q1*q3)*0.25f - (daxVar*q1*q3)*0.25f - (dayVar*q2*SG[0])*0.5f;
        SQ[8] = (dayVar*q2*q3)*0.25f - (daxVar*q1*SG[0])*0.5f - (dazVar*q2*q3)*0.25f;
        SQ[9] = sq(SG[0]);
        SQ[10] = sq(q1);

        float SPP[11];
        SPP[0] = SF[12] + SF[13] - 2*q2*SF[2];
        SPP[1] = SF[17] - SF[18] - SF[19] + SF[20];
        SPP[2] = SF[17] - SF[18] + SF[19] - SF[20];
        SPP[3] = SF[17] + SF[18] - SF[19] - SF[20];
        SPP[4] = 2*q0*q2 - 2*q1*q3;
        SPP[5] = 2*q0*q1 - 2*q2*q3;
        SPP[6] = 2*q0*q3 - 2*q1*q2;
        SPP[7] = 2*q0*q1 + 2*q2*q3;
        SPP[8] = 2*q0*q3 + 2*q1*q2;
        SPP[9] = 2*q0*q2 + 2*q1*q3;
        SPP[10] = SF[16];

        // calculate variances and upper diagonal covariances for quaternion, velocity, position and gyro bias states
        nextP(0,0) = P(0,0) + P(1,0)*SF[9] + P(2,0)*SF[11] + P(3,0)*SF[10] + P(10,0)*SF[14] + P(11,0)*SF[15] + P(12,0)*SPP[10] + (daxVar*SQ[10])*0.25f + SF[9]*(P(0,1) + P(1,1)*SF[9] + P(2,1)*SF[11] + P(3,1)*SF[10] + P(10,1)*SF[14] + P(11,1)*SF[15] + P(12,1)*SPP[10]) + SF[11]*(P(0,2) + P(1,2)*SF[9] + P(2,2)*SF[11] + P(3,2)*SF[10] + P(10,2)*SF[14] + P(11,2)*SF[15] + P(12,2)*SPP[10]) + SF[10]*(P(0,3) + P(1,3)*SF[9] + P(2,3)*SF[11] + P(3,3)*SF[10] + P(10,3)*SF[14] + P(11,3)*SF[15] + P(12,3)*SPP[10]) + SF[14]*(P(0,10) + P(1,10)*SF[9] + P(2,10)*SF[11] + P(3,10)*SF[10] + P(10,10)*SF[14] + P(11,10)*SF[15] + P(12,10)*SPP[10]) + SF[15]*(P(0,11) + P(1,11)*SF[9] + P(2,11)*SF[11] + P(3,11)*SF[10] + P(10,11)*SF[14] + P(11,11)*SF[15] + P(12,11)*SPP[10]) + SPP[10]*(P(0,12) + P(1,12)*SF[9] + P(2,12)*SF[11] + P(3,12)*SF[10] + P(10,12)*SF[14] + P(11,12)*SF[15] + P(12,12)*SPP[10]) + (dayVar*sq(q2))*0.25f + (dazVar*sq(q3))*0.25f;
        nextP(0,1) = P(0,1) + SQ[8] + P(1,1)*SF[9] + P(2,1)*SF[11] + P(3,1)*SF[10] + P(10,1)*SF[14] + P(11,1)*SF[15] + P(12,1)*SPP[10] + SF[8]*(P(0,0) + P(1,0)*SF[9] + P(2,0)*SF[11] + P(3,0)*SF[10] + P(10,0)*SF[14] + P(11,0)*SF[15] + P(12,0)*SPP[10]) + SF[7]*(P(0,2) + P(1,2)*SF[9] + P(2,2)*SF[11] + P(3,2)*SF[10] + P(10,2)*SF[14] + P(11,2)*SF[15] + P(12,2)*SPP[10]) + SF[11]*(P(0,3) + P(1,3)*SF[9] + P(2,3)*SF[11] + P(3,3)*SF[10] + P(10,3)*SF[14] + P(11,3)*SF[15] + P(12,3)*SPP[10]) - SF[15]*(P(0,12) + P(1,12)*SF[9] + P(2,12)*SF[11] + P(3,12)*SF[10] + P(10,12)*SF[14] + P(11,12)*SF[15] + P(12,12)*SPP[10]) + SPP[10]*(P(0,11) + P(1,11)*SF[9] + P(2,11)*SF[11] + P(3,11)*SF[10] + P(10,11)*SF[14] + P(11,11)*SF[15] + P(12,11)*SPP[10]) - (q0*(P(0,10) + P(1,10)*SF[9] + P(2,10)*SF[11] + P(3,10)*SF[10] + P(10,10)*SF[14] + P(11,10)*SF[15] + P(12,10)*SPP[10]))*0.5f;
        nextP(1,1) = P(1,1) + P(0,1)*SF[8] + P(2,1)*SF[7] + P(3,1)*SF[11] - P(12,1)*SF[15] + P(11,1)*SPP[10] + daxVar*SQ[9] - (P(10,1)*q0)*0.5f + SF[8]*(P(1,0) + P(0,0)*SF[8] + P(2,0)*SF[7] + P(3,0)*SF[11] - P(12,0)*SF[15] + P(11,0)*SPP[10] - (P(10,0)*q0)*0.5f) + SF[7]*(P(1,2) + P(0,2)*SF[8] + P(2,2)*SF[7] + P(3,2)*SF[11] - P(12,2)*SF[15] + P(11,2)*SPP[10] - (P(10,2)*q0)*0.5f) + SF[11]*(P(1,3) + P(0,3)*SF[8] + P(2,3)*SF[7] + P(3,3)*SF[11] - P(12,3)*SF[15] + P(11,3)*SPP[10] - (P(10,3)*q0)*0.5f) - SF[15]*(P(1,12) + P(0,12)*SF[8] + P(2,12)*SF[7] + P(3,12)*SF[11] - P(12,12)*SF[15] + P(11,12)*SPP[10] - (P(10,12)*q0)*0.5f) + SPP[10]*(P(1,11) + P(0,11)*SF[8] + P(2,11)*SF[7] + P(3,11)*SF[11] - P(12,11)*SF[15] + P(11,11)*SPP[10] - (P(10,11)*q0)*0.5f) + (dayVar*sq(q3))*0.25f + (dazVar*sq(q2))*0.25f - (q0*(P(1,10) + P(0,10)*SF[8] + P(2,10)*SF[7] + P(3,10)*SF[11] - P(12,10)*SF[15] + P(11,10)*SPP[10] - (P(10,10)*q0)*0.5f))*0.5f;
        nextP(0,2) = P(0,2) + SQ[7] + P(1,2)*SF[9] + P(2,2)*SF[11] + P(3,2)*SF[10] + P(10,2)*SF[14] + P(11,2)*SF[15] + P(12,2)*SPP[10] + SF[6]*(P(0,0) + P(1,0)*SF[9] + P(2,0)*SF[11] + P(3,0)*SF[10] + P(10,0)*SF[14] + P(11,0)*SF[15] + P(12,0)*SPP[10]) + SF[10]*(P(0,1) + P(1,1)*SF[9] + P(2,1)*SF[11] + P(3,1)*SF[10] + P(10,1)*SF[14] + P(11,1)*SF[15] + P(12,1)*SPP[10]) + SF[8]*(P(0,3) + P(1,3)*SF[9] + P(2,3)*SF[11] + P(3,3)*SF[10] + P(10,3)*SF[14] + P(11,3)*SF[15] + P(12,3)*SPP[10]) + SF[14]*(P(0,12) + P(1,12)*SF[9] + P(2,12)*SF[11] + P(3,12)*SF[10] + P(10,12)*SF[14] + P(11,12)*SF[15] + P(12,12)*SPP[10]) - SPP[10]*(P(0,10) + P(1,10)*SF[9] + P(2,10)*SF[11] + P(3,10)*SF[10] + P(10,10)*SF[14] + P(11,10)*SF[15] + P(12,10)*SPP[10]) - (q0*(P(0,11) + P(1,11)*SF[9] + P(2,11)*SF[11] + P(3,11)*SF[10] + P(10,11)*SF[14] + P(11,11)*SF[15] + P(12,11)*SPP[10]))*0.5f;
        nextP(1,2) = P(1,2) + SQ[5] + P(0,2)*SF[8] + P(2,2)*SF[7] + P(3,2)*SF[11] - P(12,2)*SF[15] + P(11,2)*SPP[10] - (P(10,2)*q0)*0.5f + SF[6]*(P(1,0) + P(0,0)*SF[8] + P(2,0)*SF[7] + P(3,0)*SF[11] - P(12,0)*SF[15] + P(11,0)*SPP[10] - (P(10,0)*q0)*0.5f) + SF[10]*(P(1,1) + P(0,1)*SF[8] + P(2,1)*SF[7] + P(3,1)*SF[11] - P(12,1)*SF[15] + P(11,1)*SPP[10] - (P(10,1)*q0)*0.5f) + SF[8]*(P(1,3) + P(0,3)*SF[8] + P(2,3)*SF[7] + P(3,3)*SF[11] - P(12,3)*SF[15] + P(11,3)*SPP[10] - (P(10,3)*q0)*0.5f) + SF[14]*(P(1,12) + P(0,12)*SF[8] + P(2,12)*SF[7] + P(3,12)*SF[11] - P(12,12)*SF[15] + P(11,12)*SPP[10] - (P(10,12)*q0)*0.5f) - SPP[10]*(P(1,10) + P(0,10)*SF[8] + P(2,10)*SF[7] + P(3,10)*SF[11] - P(12,10)*SF[15] + P(11,10)*SPP[10] - (P(10,10)*q0)*0.5f) - (q0*(P(1,11) + P(0,11)*SF[8] + P(2,11)*SF[7] + P(3,11)*SF[11] - P(12,11)*SF[15] + P(11,11)*SPP[10] - (P(10,11)*q0)*0.5f))*0.5f;
        nextP(2,2) = P(2,2) + P(0,2)*SF[6] + P(1,2)*SF[10] + P(3,2)*SF[8] + P(12,2)*SF[14] - P(10,2)*SPP[10] + dayVar*SQ[9] + (dazVar*SQ[10])*0.25f - (P(11,2)*q0)*0.5f + SF[6]*(P(2,0) + P(0,0)*SF[6] + P(1,0)*SF[10] + P(3,0)*SF[8] + P(12,0)*SF[14] - P(10,0)*SPP[10] - (P(11,0)*q0)*0.5f) + SF[10]*(P(2,1) + P(0,1)*SF[6] + P(1,1)*SF[10] + P(3,1)*SF[8] + P(12,1)*SF[14] - P(10,1)*SPP[10] - (P(11,1)*q0)*0.5f) + SF[8]*(P(2,3) + P(0,3)*SF[6] + P(1,3)*SF[10] + P(3,3)*SF[8] + P(12,3)*SF[14] - P(10,3)*SPP[10] - (P(11,3)*q0)*0.5f) + SF[14]*(P(2,12) + P(0,12)*SF[6] + P(1,12)*SF[10] + P(3,12)*SF[8] + P(12,12)*SF[14] - P(10,12)*SPP[10] - (P(11,12)*q0)*0.5f) - SPP[10]*(P(2,10) + P(0,10)*SF[6] + P(1,10)*SF[10] + P(3,10)*SF[8] + P(12,10)*SF[14] - P(10,10)*SPP[10] - (P(11,10)*q0)*0.5f) + (daxVar*sq(q3))*0.25f - (q0*(P(2,11) + P(0,11)*SF[6] + P(1,11)*SF[10] + P(3,11)*SF[8] + P(12,11)*SF[14] - P(10,11)*SPP[10] - (P(11,11)*q0)*0.5f))*0.5f;
        nextP(0,3) = P(0,3) + SQ[6] + P(1,3)*SF[9] + P(2,3)*SF[11] + P(3,3)*SF[10] + P(10,3)*SF[14] + P(11,3)*SF[15] + P(12,3)*SPP[10] + SF[7]*(P(0,0) + P(1,0)*SF[9] + P(2,0)*SF[11] + P(3,0)*SF[10] + P(10,0)*SF[14] + P(11,0)*SF[15] + P(12,0)*SPP[10]) + SF[6]*(P(0,1) + P(1,1)*SF[9] + P(2,1)*SF[11] + P(3,1)*SF[10] + P(10,1)*SF[14] + P(11,1)*SF[15] + P(12,1)*SPP[10]) + SF[9]*(P(0,2) + P(1,2)*SF[9] + P(2,2)*SF[11] + P(3,2)*SF[10] + P(10,2)*SF[14] + P(11,2)*SF[15] + P(12,2)*SPP[10]) + SF[15]*(P(0,10) + P(1,10)*SF[9] + P(2,10)*SF[11] + P(3,10)*SF[10] + P(10,10)*SF[14] + P(11,10)*SF[15] + P(12,10)*SPP[10]) - SF[14]*(P(0,11) + P(1,11)*SF[9] + P(2,11)*SF[11] + P(3,11)*SF[10] + P(10,11)*SF[14] + P(11,11)*SF[15] + P(12,11)*SPP[10]) - (q0*(P(0,12) + P(1,12)*SF[9] + P(2,12)*SF[11] + P(3,12)*SF[10] + P(10,12)*SF[14] + P(11,12)*SF[15] + P(12,12)*SPP[10]))*0.5f;
        nextP(1,3) = P(1,3) + SQ[4] + P(0,3)*SF[8] + P(2,3)*SF[7] + P(3,3)*SF[11] - P(12,3)*SF[15] + P(11,3)*SPP[10] - (P(10,3)*q0)*0.5f + SF[7]*(P(1,0) + P(0,0)*SF[8] + P(2,0)*SF[7] + P(3,0)*SF[11] - P(12,0)*SF[15] + P(11,0)*SPP[10] - (P(10,0)*q0)*0.5f) + SF[6]*(P(1,1) + P(0,1)*SF[8] + P(2,1)*SF[7] + P(3,1)*SF[11] - P(12,1)*SF[15] + P(11,1)*SPP[10] - (P(10,1)*q0)*0.5f) + SF[9]*(P(1,2) + P(0,2)*SF[8] + P(2,2)*SF[7] + P(3,2)*SF[11] - P(12,2)*SF[15] + P(11,2)*SPP[10] - (P(10,2)*q0)*0.5f) + SF[15]*(P(1,10) + P(0,10)*SF[8] + P(2,10)*SF[7] + P(3,10)*SF[11] - P(12,10)*SF[15] + P(11,10)*SPP[10] - (P(10,10)*q0)*0.5f) - SF[14]*(P(1,11) + P(0,11)*SF[8] + P(2,11)*SF[7] + P(3,11)*SF[11] - P(12,11)*SF[15] + P(11,11)*SPP[10] - (P(10,11)*q0)*0.5f) - (q0*(P(1,12) + P(0,12)*SF[8] + P(2,12)*SF[7] + P(3,12)*SF[11] - P(12,12)*SF[15] + P(11,12)*SPP[10] - (P(10,12)*q0)*0.5f))*0.5f;
        nextP(2,3) = P(2,3) + SQ[3] + P(0,3)*SF[6] + P(1,3)*SF[10] + P(3,3)*SF[8] + P(12,3)*SF[14] - P(10,3)*SPP[10] - (P(11,3)*q0)*0.5f + SF[7]*(P(2,0) + P(0,0)*SF[6] + P(1,0)*SF[10] + P(3,0)*SF[8] + P(12,0)*SF[14] - P(10,0)*SPP[10] - (P(11,0)*q0)*0.5f) + SF[6]*(P(2,1) + P(0,1)*SF[6] + P(1,1)*SF[10] + P(3,1)*SF[8] + P(12,1)*SF[14] - P(10,1)*SPP[10] - (P(11,1)*q0)*0.5f) + SF[9]*(P(2,2) + P(0,2)*SF[6] + P(1,2)*SF[10] + P(3,2)*SF[8] + P(12,2)*SF[14] - P(10,2)*SPP[10] - (P(11,2)*q0)*0.5f) + SF[15]*(P(2,10) + P(0,10)*SF[6] + P(1,10)*SF[10] + P(3,10)*SF[8] + P(12,10)*SF[14] - P(10,10)*SPP[10] - (P(11,10)*q0)*0.5f) - SF[14]*(P(2,11) + P(0,11)*SF[6] + P(1,11)*SF[10] + P(3,11)*SF[8] + P(12,11)*SF[14] - P(10,11)*SPP[10] - (P(11,11)*q0)*0.5f) - (q0*(P(2,12) + P(0,12)*SF[6] + P(1,12)*SF[10] + P(3,12)*SF[8] + P(12,12)*SF[14] - P(10,12)*SPP[10] - (P(11,12)*q0)*0.5f))*0.5f;
        nextP(3,3) = P(3,3) + P(0,3)*SF[7] + P(1,3)*SF[6] + P(2,3)*SF[9] + P(10,3)*SF[15] - P(11,3)*SF[14] + (dayVar*SQ[10])*0.25f + dazVar*SQ[9] - (P(12,3)*q0)*0.5f + SF[7]*(P(3,0) + P(0,0)*SF[7] + P(1,0)*SF[6] + P(2,0)*SF[9] + P(10,0)*SF[15] - P(11,0)*SF[14] - (P(12,0)*q0)*0.5f) + SF[6]*(P(3,1) + P(0,1)*SF[7] + P(1,1)*SF[6] + P(2,1)*SF[9] + P(10,1)*SF[15] - P(11,1)*SF[14] - (P(12,1)*q0)*0.5f) + SF[9]*(P(3,2) + P(0,2)*SF[7] + P(1,2)*SF[6] + P(2,2)*SF[9] + P(10,2)*SF[15] - P(11,2)*SF[14] - (P(12,2)*q0)*0.5f) + SF[15]*(P(3,10) + P(0,10)*SF[7] + P(1,10)*SF[6] + P(2,10)*SF[9] + P(10,10)*SF[15] - P(11,10)*SF[14] - (P(12,10)*q0)*0.5f) - SF[14]*(P(3,11) + P(0,11)*SF[7] + P(1,11)*SF[6] + P(2,11)*SF[9] + P(10,11)*SF[15] - P(11,11)*SF[14] - (P(12,11)*q0)*0.5f) + (daxVar*sq(q2))*0.25f - (q0*(P(3,12) + P(0,12)*SF[7] + P(1,12)*SF[6] + P(2,12)*SF[9] + P(10,12)*SF[15] - P(11,12)*SF[14] - (P(12,12)*q0)*0.5f))*0.5f;
        nextP(0,4) = P(0,4) + P(1,4)*SF[9] + P(2,4)*SF[11] + P(3,4)*SF[10] + P(10,4)*SF[14] + P(11,4)*SF[15] + P(12,4)*SPP[10] + SF[5]*(P(0,0) + P(1,0)*SF[9] + P(2,0)*SF[11] + P(3,0)*SF[10] + P(10,0)*SF[14] + P(11,0)*SF[15] + P(12,0)*SPP[10]) + SF[3]*(P(0,1) + P(1,1)*SF[9] + P(2,1)*SF[11] + P(3,1)*SF[10] + P(10,1)*SF[14] + P(11,1)*SF[15] + P(12,1)*SPP[10]) - SF[4]*(P(0,3) + P(1,3)*SF[9] + P(2,3)*SF[11] + P(3,3)*SF[10] + P(10,3)*SF[14] + P(11,3)*SF[15] + P(12,3)*SPP[10]) + SPP[0]*(P(0,2) + P(1,2)*SF[9] + P(2,2)*SF[11] + P(3,2)*SF[10] + P(10,2)*SF[14] + P(11,2)*SF[15] + P(12,2)*SPP[10]) + SPP[3]*(P(0,13) + P(1,13)*SF[9] + P(2,13)*SF[11] + P(3,13)*SF[10] + P(10,13)*SF[14] + P(11,13)*SF[15] + P(12,13)*SPP[10]) + SPP[6]*(P(0,14) + P(1,14)*SF[9] + P(2,14)*SF[11] + P(3,14)*SF[10] + P(10,14)*SF[14] + P(11,14)*SF[15] + P(12,14)*SPP[10]) - SPP[9]*(P(0,15) + P(1,15)*SF[9] + P(2,15)*SF[11] + P(3,15)*SF[10] + P(10,15)*SF[14] + P(11,15)*SF[15] + P(12,15)*SPP[10]);
        nextP(1,4) = P(1,4) + P(0,4)*SF[8] + P(2,4)*SF[7] + P(3,4)*SF[11] - P(12,4)*SF[15] + P(11,4)*SPP[10] - (P(10,4)*q0)*0.5f + SF[5]*(P(1,0) + P(0,0)*SF[8] + P(2,0)*SF[7] + P(3,0)*SF[11] - P(12,0)*SF[15] + P(11,0)*SPP[10] - (P(10,0)*q0)*0.5f) + SF[3]*(P(1,1) + P(0,1)*SF[8] + P(2,1)*SF[7] + P(3,1)*SF[11] - P(12,1)*SF[15] + P(11,1)*SPP[10] - (P(10,1)*q0)*0.5f) - SF[4]*(P(1,3) + P(0,3)*SF[8] + P(2,3)*SF[7] + P(3,3)*SF[11] - P(12,3)*SF[15] + P(11,3)*SPP[10] - (P(10,3)*q0)*0.5f) + SPP[0]*(P(1,2) + P(0,2)*SF[8] + P(2,2)*SF[7] + P(3,2)*SF[11] - P(12,2)*SF[15] + P(11,2)*SPP[10] - (P(10,2)*q0)*0.5f) + SPP[3]*(P(1,13) + P(0,13)*SF[8] + P(2,13)*SF[7] + P(3,13)*SF[11] - P(12,13)*SF[15] + P(11,13)*SPP[10] - (P(10,13)*q0)*0.5f) + SPP[6]*(P(1,14) + P(0,14)*SF[8] + P(2,14)*SF[7] + P(3,14)*SF[11] - P(12,14)*SF[15] + P(11,14)*SPP[10] - (P(10,14)*q0)*0.5f) - SPP[9]*(P(1,15) + P(0,15)*SF[8] + P(2,15)*SF[7] + P(3,15)*SF[11] - P(12,15)*SF[15] + P(11,15)*SPP[10] - (P(10,15)*q0)*0.5f);
        nextP(2,4) = P(2,4) + P(0,4)*SF[6] + P(1,4)*SF[10] + P(3,4)*SF[8] + P(12,4)*SF[14] - P(10,4)*SPP[10] - (P(11,4)*q0)*0.5f + SF[5]*(P(2,0) + P(0,0)*SF[6] + P(1,0)*SF[10] + P(3,0)*SF[8] + P(12,0)*SF[14] - P(10,0)*SPP[10] - (P(11,0)*q0)*0.5f) + SF[3]*(P(2,1) + P(0,1)*SF[6] + P(1,1)*SF[10] + P(3,1)*SF[8] + P(12,1)*SF[14] - P(10,1)*SPP[10] - (P(11,1)*q0)*0.5f) - SF[4]*(P(2,3) + P(0,3)*SF[6] + P(1,3)*SF[10] + P(3,3)*SF[8] + P(12,3)*SF[14] - P(10,3)*SPP[10] - (P(11,3)*q0)*0.5f) + SPP[0]*(P(2,2) + P(0,2)*SF[6] + P(1,2)*SF[10] + P(3,2)*SF[8] + P(12,2)*SF[14] - P(10,2)*SPP[10] - (P(11,2)*q0)*0.5f) + SPP[3]*(P(2,13) + P(0,13)*SF[6] + P(1,13)*SF[10] + P(3,13)*SF[8] + P(12,13)*SF[14] - P(10,13)*SPP[10] - (P(11,13)*q0)*0.5f) + SPP[6]*(P(2,14) + P(0,14)*SF[6] + P(1,14)*SF[10] + P(3,14)*SF[8] + P(12,14)*SF[14] - P(10,14)*SPP[10] - (P(11,14)*q0)*0.5f) - SPP[9]*(P(2,15) + P(0,15)*SF[6] + P(1,15)*SF[10] + P(3,15)*SF[8] + P(12,15)*SF[14] - P(10,15)*SPP[10] - (P(11,15)*q0)*0.5f);
        nextP(3,4) = P(3,4) + P(0,4)*SF[7] + P(1,4)*SF[6] + P(2,4)*SF[9] + P(10,4)*SF[15] - P(11,4)*SF[14] - (P(12,4)*q0)*0.5f + SF[5]*(P(3,0) + P(0,0)*SF[7] + P(1,0)*SF[6] + P(2,0)*SF[9] + P(10,0)*SF[15] - P(11,0)*SF[14] - (P(12,0)*q0)*0.5f) + SF[3]*(P(3,1) + P(0,1)*SF[7] + P(1,1)*SF[6] + P(2,1)*SF[9] + P(10,1)*SF[15] - P(11,1)*SF[14] - (P(12,1)*q0)*0.5f) - SF[4]*(P(3,3) + P(0,3)*SF[7] + P(1,3)*SF[6] + P(2,3)*SF[9] + P(10,3)*SF[15] - P(11,3)*SF[14] - (P(12,3)*q0)*0.5f) + SPP[0]*(P(3,2) + P(0,2)*SF[7] + P(1,2)*SF[6] + P(2,2)*SF[9] + P(10,2)*SF[15] - P(11,2)*SF[14] - (P(12,2)*q0)*0.5f) + SPP[3]*(P(3,13) + P(0,13)*SF[7] + P(1,13)*SF[6] + P(2,13)*SF[9] + P(10,13)*SF[15] - P(11,13)*SF[14] - (P(12,13)*q0)*0.5f) + SPP[6]*(P(3,14) + P(0,14)*SF[7] + P(1,14)*SF[6] + P(2,14)*SF[9] + P(10,14)*SF[15] - P(11,14)*SF[14] - (P(12,14)*q0)*0.5f) - SPP[9]*(P(3,15) + P(0,15)*SF[7] + P(1,15)*SF[6] + P(2,15)*SF[9] + P(10,15)*SF[15] - P(11,15)*SF[14] - (P(12,15)*q0)*0.5f);
        nextP(4,4) = P(4,4) + P(0,4)*SF[5] + P(1,4)*SF[3] - P(3,4)*SF[4] + P(2,4)*SPP[0] + P(13,4)*SPP[3] + P(14,4)*SPP[6] - P(15,4)*SPP[9] + dvyVar*sq(SG[7] - 2*q0*q3) + dvzVar*sq(SG[6] + 2*q0*q2) + SF[5]*(P(4,0) + P(0,0)*SF[5] + P(1,0)*SF[3] - P(3,0)*SF[4] + P(2,0)*SPP[0] + P(13,0)*SPP[3] + P(14,0)*SPP[6] - P(15,0)*SPP[9]) + SF[3]*(P(4,1) + P(0,1)*SF[5] + P(1,1)*SF[3] - P(3,1)*SF[4] + P(2,1)*SPP[0] + P(13,1)*SPP[3] + P(14,1)*SPP[6] - P(15,1)*SPP[9]) - SF[4]*(P(4,3) + P(0,3)*SF[5] + P(1,3)*SF[3] - P(3,3)*SF[4] + P(2,3)*SPP[0] + P(13,3)*SPP[3] + P(14,3)*SPP[6] - P(15,3)*SPP[9]) + SPP[0]*(P(4,2) + P(0,2)*SF[5] + P(1,2)*SF[3] - P(3,2)*SF[4] + P(2,2)*SPP[0] + P(13,2)*SPP[3] + P(14,2)*SPP[6] - P(15,2)*SPP[9]) + SPP[3]*(P(4,13) + P(0,13)*SF[5] + P(1,13)*SF[3] - P(3,13)*SF[4] + P(2,13)*SPP[0] + P(13,13)*SPP[3] + P(14,13)*SPP[6] - P(15,13)*SPP[9]) + SPP[6]*(P(4,14) + P(0,14)*SF[5] + P(1,14)*SF[3] - P(3,14)*SF[4] + P(2,14)*SPP[0] + P(13,14)*SPP[3] + P(14,14)*SPP[6] - P(15,14)*SPP[9]) - SPP[9]*(P(4,15) + P(0,15)*SF[5] + P(1,15)*SF[3] - P(3,15)*SF[4] + P(2,15)*SPP[0] + P(13,15)*SPP[3] + P(14,15)*SPP[6] - P(15,15)*SPP[9]) + dvxVar*sq(SG[1] + SG[2] - SG[3] - SG[4]);
        nextP(0,5) = P(0,5) + P(1,5)*SF[9] + P(2,5)*SF[11] + P(3,5)*SF[10] + P(10,5)*SF[14] + P(11,5)*SF[15] + P(12,5)*SPP[10] + SF[4]*(P(0,0) + P(1,0)*SF[9] + P(2,0)*SF[11] + P(3,0)*SF[10] + P(10,0)*SF[14] + P(11,0)*SF[15] + P(12,0)*SPP[10]) + SF[3]*(P(0,2) + P(1,2)*SF[9] + P(2,2)*SF[11] + P(3,2)*SF[10] + P(10,2)*SF[14] + P(11,2)*SF[15] + P(12,2)*SPP[10]) + SF[5]*(P(0,3) + P(1,3)*SF[9] + P(2,3)*SF[11] + P(3,3)*SF[10] + P(10,3)*SF[14] + P(11,3)*SF[15] + P(12,3)*SPP[10]) - SPP[0]*(P(0,1) + P(1,1)*SF[9] + P(2,1)*SF[11] + P(3,1)*SF[10] + P(10,1)*SF[14] + P(11,1)*SF[15] + P(12,1)*SPP[10]) - SPP[8]*(P(0,13) + P(1,13)*SF[9] + P(2,13)*SF[11] + P(3,13)*SF[10] + P(10,13)*SF[14] + P(11,13)*SF[15] + P(12,13)*SPP[10]) + SPP[2]*(P(0,14) + P(1,14)*SF[9] + P(2,14)*SF[11] + P(3,14)*SF[10] + P(10,14)*SF[14] + P(11,14)*SF[15] + P(12,14)*SPP[10]) + SPP[5]*(P(0,15) + P(1,15)*SF[9] + P(2,15)*SF[11] + P(3,15)*SF[10] + P(10,15)*SF[14] + P(11,15)*SF[15] + P(12,15)*SPP[10]);
        nextP(1,5) = P(1,5) + P(0,5)*SF[8] + P(2,5)*SF[7] + P(3,5)*SF[11] - P(12,5)*SF[15] + P(11,5)*SPP[10] - (P(10,5)*q0)*0.5f + SF[4]*(P(1,0) + P(0,0)*SF[8] + P(2,0)*SF[7] + P(3,0)*SF[11] - P(12,0)*SF[15] + P(11,0)*SPP[10] - (P(10,0)*q0)*0.5f) + SF[3]*(P(1,2) + P(0,2)*SF[8] + P(2,2)*SF[7] + P(3,2)*SF[11] - P(12,2)*SF[15] + P(11,2)*SPP[10] - (P(10,2)*q0)*0.5f) + SF[5]*(P(1,3) + P(0,3)*SF[8] + P(2,3)*SF[7] + P(3,3)*SF[11] - P(12,3)*SF[15] + P(11,3)*SPP[10] - (P(10,3)*q0)*0.5f) - SPP[0]*(P(1,1) + P(0,1)*SF[8] + P(2,1)*SF[7] + P(3,1)*SF[11] - P(12,1)*SF[15] + P(11,1)*SPP[10] - (P(10,1)*q0)*0.5f) - SPP[8]*(P(1,13) + P(0,13)*SF[8] + P(2,13)*SF[7] + P(3,13)*SF[11] - P(12,13)*SF[15] + P(11,13)*SPP[10] - (P(10,13)*q0)*0.5f) + SPP[2]*(P(1,14) + P(0,14)*SF[8] + P(2,14)*SF[7] + P(3,14)*SF[11] - P(12,14)*SF[15] + P(11,14)*SPP[10] - (P(10,14)*q0)*0.5f) + SPP[5]*(P(1,15) + P(0,15)*SF[8] + P(2,15)*SF[7] + P(3,15)*SF[11] - P(12,15)*SF[15] + P(11,15)*SPP[10] - (P(10,15)*q0)*0.5f);
        nextP(2,5) = P(2,5) + P(0,5)*SF[6] + P(1,5)*SF[10] + P(3,5)*SF[8] + P(12,5)*SF[14] - P(10,5)*SPP[10] - (P(11,5)*q0)*0.5f + SF[4]*(P(2,0) + P(0,0)*SF[6] + P(1,0)*SF[10] + P(3,0)*SF[8] + P(12,0)*SF[14] - P(10,0)*SPP[10] - (P(11,0)*q0)*0.5f) + SF[3]*(P(2,2) + P(0,2)*SF[6] + P(1,2)*SF[10] + P(3,2)*SF[8] + P(12,2)*SF[14] - P(10,2)*SPP[10] - (P(11,2)*q0)*0.5f) + SF[5]*(P(2,3) + P(0,3)*SF[6] + P(1,3)*SF[10] + P(3,3)*SF[8] + P(12,3)*SF[14] - P(10,3)*SPP[10] - (P(11,3)*q0)*0.5f) - SPP[0]*(P(2,1) + P(0,1)*SF[6] + P(1,1)*SF[10] + P(3,1)*SF[8] + P(12,1)*SF[14] - P(10,1)*SPP[10] - (P(11,1)*q0)*0.5f) - SPP[8]*(P(2,13) + P(0,13)*SF[6] + P(1,13)*SF[10] + P(3,13)*SF[8] + P(12,13)*SF[14] - P(10,13)*SPP[10] - (P(11,13)*q0)*0.5f) + SPP[2]*(P(2,14) + P(0,14)*SF[6] + P(1,14)*SF[10] + P(3,14)*SF[8] + P(12,14)*SF[14] - P(10,14)*SPP[10] - (P(11,14)*q0)*0.5f) + SPP[5]*(P(2,15) + P(0,15)*SF[6] + P(1,15)*SF[10] + P(3,15)*SF[8] + P(12,15)*SF[14] - P(10,15)*SPP[10] - (P(11,15)*q0)*0.5f);
        nextP(3,5) = P(3,5) + P(0,5)*SF[7] + P(1,5)*SF[6] + P(2,5)*SF[9] + P(10,5)*SF[15] - P(11,5)*SF[14] - (P(12,5)*q0)*0.5f + SF[4]*(P(3,0) + P(0,0)*SF[7] + P(1,0)*SF[6] + P(2,0)*SF[9] + P(10,0)*SF[15] - P(11,0)*SF[14] - (P(12,0)*q0)*0.5f) + SF[3]*(P(3,2) + P(0,2)*SF[7] + P(1,2)*SF[6] + P(2,2)*SF[9] + P(10,2)*SF[15] - P(11,2)*SF[14] - (P(12,2)*q0)*0.5f) + SF[5]*(P(3,3) + P(0,3)*SF[7] + P(1,3)*SF[6] + P(2,3)*SF[9] + P(10,3)*SF[15] - P(11,3)*SF[14] - (P(12,3)*q0)*0.5f) - SPP[0]*(P(3,1) + P(0,1)*SF[7] + P(1,1)*SF[6] + P(2,1)*SF[9] + P(10,1)*SF[15] - P(11,1)*SF[14] - (P(12,1)*q0)*0.5f) - SPP[8]*(P(3,13) + P(0,13)*SF[7] + P(1,13)*SF[6] + P(2,13)*SF[9] + P(10,13)*SF[15] - P(11,13)*SF[14] - (P(12,13)*q0)*0.5f) + SPP[2]*(P(3,14) + P(0,14)*SF[7] + P(1,14)*SF[6] + P(2,14)*SF[9] + P(10,14)*SF[15] - P(11,14)*SF[14] - (P(12,14)*q0)*0.5f) + SPP[5]*(P(3,15) + P(0,15)*SF[7] + P(1,15)*SF[6] + P(2,15)*SF[9] + P(10,15)*SF[15] - P(11,15)*SF[14] - (P(12,15)*q0)*0.5f);
        nextP(4,5) = P(4,5) + SQ[2] + P(0,5)*SF[5] + P(1,5)*SF[3] - P(3,5)*SF[4] + P(2,5)*SPP[0] + P(13,5)*SPP[3] + P(14,5)*SPP[6] - P(15,5)*SPP[9] + SF[4]*(P(4,0) + P(0,0)*SF[5] + P(1,0)*SF[3] - P(3,0)*SF[4] + P(2,0)*SPP[0] + P(13,0)*SPP[3] + P(14,0)*SPP[6] - P(15,0)*SPP[9]) + SF[3]*(P(4,2) + P(0,2)*SF[5] + P(1,2)*SF[3] - P(3,2)*SF[4] + P(2,2)*SPP[0] + P(13,2)*SPP[3] + P(14,2)*SPP[6] - P(15,2)*SPP[9]) + SF[5]*(P(4,3) + P(0,3)*SF[5] + P(1,3)*SF[3] - P(3,3)*SF[4] + P(2,3)*SPP[0] + P(13,3)*SPP[3] + P(14,3)*SPP[6] - P(15,3)*SPP[9]) - SPP[0]*(P(4,1) + P(0,1)*SF[5] + P(1,1)*SF[3] - P(3,1)*SF[4] + P(2,1)*SPP[0] + P(13,1)*SPP[3] + P(14,1)*SPP[6] - P(15,1)*SPP[9]) - SPP[8]*(P(4,13) + P(0,13)*SF[5] + P(1,13)*SF[3] - P(3,13)*SF[4] + P(2,13)*SPP[0] + P(13,13)*SPP[3] + P(14,13)*SPP[6] - P(15,13)*SPP[9]) + SPP[2]*(P(4,14) + P(0,14)*SF[5] + P(1,14)*SF[3] - P(3,14)*SF[4] + P(2,14)*SPP[0] + P(13,14)*SPP[3] + P(14,14)*SPP[6] - P(15,14)*SPP[9]) + SPP[5]*(P(4,15) + P(0,15)*SF[5] + P(1,15)*SF[3] - P(3,15)*SF[4] + P(2,15)*SPP[0] + P(13,15)*SPP[3] + P(14,15)*SPP[6] - P(15,15)*SPP[9]);
        nextP(5,5) = P(5,5) + P(0,5)*SF[4] + P(2,5)*SF[3] + P(3,5)*SF[5] - P(1,5)*SPP[0] - P(13,5)*SPP[8] + P(14,5)*SPP[2] + P(15,5)*SPP[5] + dvxVar*sq(SG[7] + 2*q0*q3) + dvzVar*sq(SG[5] - 2*q0*q1) + SF[4]*(P(5,0) + P(0,0)*SF[4] + P(2,0)*SF[3] + P(3,0)*SF[5] - P(1,0)*SPP[0] - P(13,0)*SPP[8] + P(14,0)*SPP[2] + P(15,0)*SPP[5]) + SF[3]*(P(5,2) + P(0,2)*SF[4] + P(2,2)*SF[3] + P(3,2)*SF[5] - P(1,2)*SPP[0] - P(13,2)*SPP[8] + P(14,2)*SPP[2] + P(15,2)*SPP[5]) + SF[5]*(P(5,3) + P(0,3)*SF[4] + P(2,3)*SF[3] + P(3,3)*SF[5] - P(1,3)*SPP[0] - P(13,3)*SPP[8] + P(14,3)*SPP[2] + P(15,3)*SPP[5]) - SPP[0]*(P(5,1) + P(0,1)*SF[4] + P(2,1)*SF[3] + P(3,1)*SF[5] - P(1,1)*SPP[0] - P(13,1)*SPP[8] + P(14,1)*SPP[2] + P(15,1)*SPP[5]) - SPP[8]*(P(5,13) + P(0,13)*SF[4] + P(2,13)*SF[3] + P(3,13)*SF[5] - P(1,13)*SPP[0] - P(13,13)*SPP[8] + P(14,13)*SPP[2] + P(15,13)*SPP[5]) + SPP[2]*(P(5,14) + P(0,14)*SF[4] + P(2,14)*SF[3] + P(3,14)*SF[5] - P(1,14)*SPP[0] - P(13,14)*SPP[8] + P(14,14)*SPP[2] + P(15,14)*SPP[5]) + SPP[5]*(P(5,15) + P(0,15)*SF[4] + P(2,15)*SF[3] + P(3,15)*SF[5] - P(1,15)*SPP[0] - P(13,15)*SPP[8] + P(14,15)*SPP[2] + P(15,15)*SPP[5]) + dvyVar*sq(SG[1] - SG[2] + SG[3] - SG[4]);
        nextP(0,6) = P(0,6) + P(1,6)*SF[9] + P(2,6)*SF[11] + P(3,6)*SF[10] + P(10,6)*SF[14] + P(11,6)*SF[15] + P(12,6)*SPP[10] + SF[4]*(P(0,1) + P(1,1)*SF[9] + P(2,1)*SF[11] + P(3,1)*SF[10] + P(10,1)*SF[14] + P(11,1)*SF[15] + P(12,1)*SPP[10]) - SF[5]*(P(0,2) + P(1,2)*SF[9] + P(2,2)*SF[11] + P(3,2)*SF[10] + P(10,2)*SF[14] + P(11,2)*SF[15] + P(12,2)*SPP[10]) + SF[3]*(P(0,3) + P(1,3)*SF[9] + P(2,3)*SF[11] + P(3,3)*SF[10] + P(10,3)*SF[14] + P(11,3)*SF[15] + P(12,3)*SPP[10]) + SPP[0]*(P(0,0) + P(1,0)*SF[9] + P(2,0)*SF[11] + P(3,0)*SF[10] + P(10,0)*SF[14] + P(11,0)*SF[15] + P(12,0)*SPP[10]) + SPP[4]*(P(0,13) + P(1,13)*SF[9] + P(2,13)*SF[11] + P(3,13)*SF[10] + P(10,13)*SF[14] + P(11,13)*SF[15] + P(12,13)*SPP[10]) - SPP[7]*(P(0,14) + P(1,14)*SF[9] + P(2,14)*SF[11] + P(3,14)*SF[10] + P(10,14)*SF[14] + P(11,14)*SF[15] + P(12,14)*SPP[10]) - SPP[1]*(P(0,15) + P(1,15)*SF[9] + P(2,15)*SF[11] + P(3,15)*SF[10] + P(10,15)*SF[14] + P(11,15)*SF[15] + P(12,15)*SPP[10]);
        nextP(1,6) = P(1,6) + P(0,6)*SF[8] + P(2,6)*SF[7] + P(3,6)*SF[11] - P(12,6)*SF[15] + P(11,6)*SPP[10] - (P(10,6)*q0)*0.5f + SF[4]*(P(1,1) + P(0,1)*SF[8] + P(2,1)*SF[7] + P(3,1)*SF[11] - P(12,1)*SF[15] + P(11,1)*SPP[10] - (P(10,1)*q0)*0.5f) - SF[5]*(P(1,2) + P(0,2)*SF[8] + P(2,2)*SF[7] + P(3,2)*SF[11] - P(12,2)*SF[15] + P(11,2)*SPP[10] - (P(10,2)*q0)*0.5f) + SF[3]*(P(1,3) + P(0,3)*SF[8] + P(2,3)*SF[7] + P(3,3)*SF[11] - P(12,3)*SF[15] + P(11,3)*SPP[10] - (P(10,3)*q0)*0.5f) + SPP[0]*(P(1,0) + P(0,0)*SF[8] + P(2,0)*SF[7] + P(3,0)*SF[11] - P(12,0)*SF[15] + P(11,0)*SPP[10] - (P(10,0)*q0)*0.5f) + SPP[4]*(P(1,13) + P(0,13)*SF[8] + P(2,13)*SF[7] + P(3,13)*SF[11] - P(12,13)*SF[15] + P(11,13)*SPP[10] - (P(10,13)*q0)*0.5f) - SPP[7]*(P(1,14) + P(0,14)*SF[8] + P(2,14)*SF[7] + P(3,14)*SF[11] - P(12,14)*SF[15] + P(11,14)*SPP[10] - (P(10,14)*q0)*0.5f) - SPP[1]*(P(1,15) + P(0,15)*SF[8] + P(2,15)*SF[7] + P(3,15)*SF[11] - P(12,15)*SF[15] + P(11,15)*SPP[10] - (P(10,15)*q0)*0.5f);
        nextP(2,6) = P(2,6) + P(0,6)*SF[6] + P(1,6)*SF[10] + P(3,6)*SF[8] + P(12,6)*SF[14] - P(10,6)*SPP[10] - (P(11,6)*q0)*0.5f + SF[4]*(P(2,1) + P(0,1)*SF[6] + P(1,1)*SF[10] + P(3,1)*SF[8] + P(12,1)*SF[14] - P(10,1)*SPP[10] - (P(11,1)*q0)*0.5f) - SF[5]*(P(2,2) + P(0,2)*SF[6] + P(1,2)*SF[10] + P(3,2)*SF[8] + P(12,2)*SF[14] - P(10,2)*SPP[10] - (P(11,2)*q0)*0.5f) + SF[3]*(P(2,3) + P(0,3)*SF[6] + P(1,3)*SF[10] + P(3,3)*SF[8] + P(12,3)*SF[14] - P(10,3)*SPP[10] - (P(11,3)*q0)*0.5f) + SPP[0]*(P(2,0) + P(0,0)*SF[6] + P(1,0)*SF[10] + P(3,0)*SF[8] + P(12,0)*SF[14] - P(10,0)*SPP[10] - (P(11,0)*q0)*0.5f) + SPP[4]*(P(2,13) + P(0,13)*SF[6] + P(1,13)*SF[10] + P(3,13)*SF[8] + P(12,13)*SF[14] - P(10,13)*SPP[10] - (P(11,13)*q0)*0.5f) - SPP[7]*(P(2,14) + P(0,14)*SF[6] + P(1,14)*SF[10] + P(3,14)*SF[8] + P(12,14)*SF[14] - P(10,14)*SPP[10] - (P(11,14)*q0)*0.5f) - SPP[1]*(P(2,15) + P(0,15)*SF[6] + P(1,15)*SF[10] + P(3,15)*SF[8] + P(12,15)*SF[14] - P(10,15)*SPP[10] - (P(11,15)*q0)*0.5f);
        nextP(3,6) = P(3,6) + P(0,6)*SF[7] + P(1,6)*SF[6] + P(2,6)*SF[9] + P(10,6)*SF[15] - P(11,6)*SF[14] - (P(12,6)*q0)*0.5f + SF[4]*(P(3,1) + P(0,1)*SF[7] + P(1,1)*SF[6] + P(2,1)*SF[9] + P(10,1)*SF[15] - P(11,1)*SF[14] - (P(12,1)*q0)*0.5f) - SF[5]*(P(3,2) + P(0,2)*SF[7] + P(1,2)*SF[6] + P(2,2)*SF[9] + P(10,2)*SF[15] - P(11,2)*SF[14] - (P(12,2)*q0)*0.5f) + SF[3]*(P(3,3) + P(0,3)*SF[7] + P(1,3)*SF[6] + P(2,3)*SF[9] + P(10,3)*SF[15] - P(11,3)*SF[14] - (P(12,3)*q0)*0.5f) + SPP[0]*(P(3,0) + P(0,0)*SF[7] + P(1,0)*SF[6] + P(2,0)*SF[9] + P(10,0)*SF[15] - P(11,0)*SF[14] - (P(12,0)*q0)*0.5f) + SPP[4]*(P(3,13) + P(0,13)*SF[7] + P(1,13)*SF[6] + P(2,13)*SF[9] + P(10,13)*SF[15] - P(11,13)*SF[14] - (P(12,13)*q0)*0.5f) - SPP[7]*(P(3,14) + P(0,14)*SF[7] + P(1,14)*SF[6] + P(2,14)*SF[9] + P(10,14)*SF[15] - P(11,14)*SF[14] - (P(12,14)*q0)*0.5f) - SPP[1]*(P(3,15) + P(0,15)*SF[7] + P(1,15)*SF[6] + P(2,15)*SF[9] + P(10,15)*SF[15] - P(11,15)*SF[14] - (P(12,15)*q0)*0.5f);
        nextP(4,6) = P(4,6) + SQ[1] + P(0,6)*SF[5] + P(1,6)*SF[3] - P(3,6)*SF[4] + P(2,6)*SPP[0] + P(13,6)*SPP[3] + P(14,6)*SPP[6] - P(15,6)*SPP[9] + SF[4]*(P(4,1) + P(0,1)*SF[5] + P(1,1)*SF[3] - P(3,1)*SF[4] + P(2,1)*SPP[0] + P(13,1)*SPP[3] + P(14,1)*SPP[6] - P(15,1)*SPP[9]) - SF[5]*(P(4,2) + P(0,2)*SF[5] + P(1,2)*SF[3] - P(3,2)*SF[4] + P(2,2)*SPP[0] + P(13,2)*SPP[3] + P(14,2)*SPP[6] - P(15,2)*SPP[9]) + SF[3]*(P(4,3) + P(0,3)*SF[5] + P(1,3)*SF[3] - P(3,3)*SF[4] + P(2,3)*SPP[0] + P(13,3)*SPP[3] + P(14,3)*SPP[6] - P(15,3)*SPP[9]) + SPP[0]*(P(4,0) + P(0,0)*SF[5] + P(1,0)*SF[3] - P(3,0)*SF[4] + P(2,0)*SPP[0] + P(13,0)*SPP[3] + P(14,0)*SPP[6] - P(15,0)*SPP[9]) + SPP[4]*(P(4,13) + P(0,13)*SF[5] + P(1,13)*SF[3] - P(3,13)*SF[4] + P(2,13)*SPP[0] + P(13,13)*SPP[3] + P(14,13)*SPP[6] - P(15,13)*SPP[9]) - SPP[7]*(P(4,14) + P(0,14)*SF[5] + P(1,14)*SF[3] - P(3,14)*SF[4] + P(2,14)*SPP[0] + P(13,14)*SPP[3] + P(14,14)*SPP[6] - P(15,14)*SPP[9]) - SPP[1]*(P(4,15) + P(0,15)*SF[5] + P(1,15)*SF[3] - P(3,15)*SF[4] + P(2,15)*SPP[0] + P(13,15)*SPP[3] + P(14,15)*SPP[6] - P(15,15)*SPP[9]);
        nextP(5,6) = P(5,6) + SQ[0] + P(0,6)*SF[4] + P(2,6)*SF[3] + P(3,6)*SF[5] - P(1,6)*SPP[0] - P(13,6)*SPP[8] + P(14,6)*SPP[2] + P(15,6)*SPP[5] + SF[4]*(P(5,1) + P(0,1)*SF[4] + P(2,1)*SF[3] + P(3,1)*SF[5] - P(1,1)*SPP[0] - P(13,1)*SPP[8] + P(14,1)*SPP[2] + P(15,1)*SPP[5]) - SF[5]*(P(5,2) + P(0,2)*SF[4] + P(2,2)*SF[3] + P(3,2)*SF[5] - P(1,2)*SPP[0] - P(13,2)*SPP[8] + P(14,2)*SPP[2] + P(15,2)*SPP[5]) + SF[3]*(P(5,3) + P(0,3)*SF[4] + P(2,3)*SF[3] + P(3,3)*SF[5] - P(1,3)*SPP[0] - P(13,3)*SPP[8] + P(14,3)*SPP[2] + P(15,3)*SPP[5]) + SPP[0]*(P(5,0) + P(0,0)*SF[4] + P(2,0)*SF[3] + P(3,0)*SF[5] - P(1,0)*SPP[0] - P(13,0)*SPP[8] + P(14,0)*SPP[2] + P(15,0)*SPP[5]) + SPP[4]*(P(5,13) + P(0,13)*SF[4] + P(2,13)*SF[3] + P(3,13)*SF[5] - P(1,13)*SPP[0] - P(13,13)*SPP[8] + P(14,13)*SPP[2] + P(15,13)*SPP[5]) - SPP[7]*(P(5,14) + P(0,14)*SF[4] + P(2,14)*SF[3] + P(3,14)*SF[5] - P(1,14)*SPP[0] - P(13,14)*SPP[8] + P(14,14)*SPP[2] + P(15,14)*SPP[5]) - SPP[1]*(P(5,15) + P(0,15)*SF[4] + P(2,15)*SF[3] + P(3,15)*SF[5] - P(1,15)*SPP[0] - P(13,15)*SPP[8] + P(14,15)*SPP[2] + P(15,15)*SPP[5]);
        nextP(6,6) = P(6,6) + P(1,6)*SF[4] - P(2,6)*SF[5] + P(3,6)*SF[3] + P(0,6)*SPP[0] + P(13,6)*SPP[4] - P(14,6)*SPP[7] - P(15,6)*SPP[1] + dvxVar*sq(SG[6] - 2*q0*q2) + dvyVar*sq(SG[5] + 2*q0*q1) + SF[4]*(P(6,1) + P(1,1)*SF[4] - P(2,1)*SF[5] + P(3,1)*SF[3] + P(0,1)*SPP[0] + P(13,1)*SPP[4] - P(14,1)*SPP[7] - P(15,1)*SPP[1]) - SF[5]*(P(6,2) + P(1,2)*SF[4] - P(2,2)*SF[5] + P(3,2)*SF[3] + P(0,2)*SPP[0] + P(13,2)*SPP[4] - P(14,2)*SPP[7] - P(15,2)*SPP[1]) + SF[3]*(P(6,3) + P(1,3)*SF[4] - P(2,3)*SF[5] + P(3,3)*SF[3] + P(0,3)*SPP[0] + P(13,3)*SPP[4] - P(14,3)*SPP[7] - P(15,3)*SPP[1]) + SPP[0]*(P(6,0) + P(1,0)*SF[4] - P(2,0)*SF[5] + P(3,0)*SF[3] + P(0,0)*SPP[0] + P(13,0)*SPP[4] - P(14,0)*SPP[7] - P(15,0)*SPP[1]) + SPP[4]*(P(6,13) + P(1,13)*SF[4] - P(2,13)*SF[5] + P(3,13)*SF[3] + P(0,13)*SPP[0] + P(13,13)*SPP[4] - P(14,13)*SPP[7] - P(15,13)*SPP[1]) - SPP[7]*(P(6,14) + P(1,14)*SF[4] - P(2,14)*SF[5] + P(3,14)*SF[3] + P(0,14)*SPP[0] + P(13,14)*SPP[4] - P(14,14)*SPP[7] - P(15,14)*SPP[1]) - SPP[1]*(P(6,15) + P(1,15)*SF[4] - P(2,15)*SF[5] + P(3,15)*SF[3] + P(0,15)*SPP[0] + P(13,15)*SPP[4] - P(14,15)*SPP[7] - P(15,15)*SPP[1]) + dvzVar*sq(SG[1] - SG[2] - SG[3] + SG[4]);
        nextP(0,7) = P(0,7) + P(1,7)*SF[9] + P(2,7)*SF[11] + P(3,7)*SF[10] + P(10,7)*SF[14] + P(11,7)*SF[15] + P(12,7)*SPP[10] + dt*(P(0,4) + P(1,4)*SF[9] + P(2,4)*SF[11] + P(3,4)*SF[10] + P(10,4)*SF[14] + P(11,4)*SF[15] + P(12,4)*SPP[10]);
        nextP(1,7) = P(1,7) + P(0,7)*SF[8] + P(2,7)*SF[7] + P(3,7)*SF[11] - P(12,7)*SF[15] + P(11,7)*SPP[10] - (P(10,7)*q0)*0.5f + dt*(P(1,4) + P(0,4)*SF[8] + P(2,4)*SF[7] + P(3,4)*SF[11] - P(12,4)*SF[15] + P(11,4)*SPP[10] - (P(10,4)*q0)*0.5f);
        nextP(2,7) = P(2,7) + P(0,7)*SF[6] + P(1,7)*SF[10] + P(3,7)*SF[8] + P(12,7)*SF[14] - P(10,7)*SPP[10] - (P(11,7)*q0)*0.5f + dt*(P(2,4) + P(0,4)*SF[6] + P(1,4)*SF[10] + P(3,4)*SF[8] + P(12,4)*SF[14] - P(10,4)*SPP[10] - (P(11,4)*q0)*0.5f);
        nextP(3,7) = P(3,7) + P(0,7)*SF[7] + P(1,7)*SF[6] + P(2,7)*SF[9] + P(10,7)*SF[15] - P(11,7)*SF[14] - (P(12,7)*q0)*0.5f + dt*(P(3,4) + P(0,4)*SF[7] + P(1,4)*SF[6] + P(2,4)*SF[9] + P(10,4)*SF[15] - P(11,4)*SF[14] - (P(12,4)*q0)*0.5f);
        nextP(4,7) = P(4,7) + P(0,7)*SF[5] + P(1,7)*SF[3] - P(3,7)*SF[4] + P(2,7)*SPP[0] + P(13,7)*SPP[3] + P(14,7)*SPP[6] - P(15,7)*SPP[9] + dt*(P(4,4) + P(0,4)*SF[5] + P(1,4)*SF[3] - P(3,4)*SF[4] + P(2,4)*SPP[0] + P(13,4)*SPP[3] + P(14,4)*SPP[6] - P(15,4)*SPP[9]);
        nextP(5,7) = P(5,7) + P(0,7)*SF[4] + P(2,7)*SF[3] + P(3,7)*SF[5] - P(1,7)*SPP[0] - P(13,7)*SPP[8] + P(14,7)*SPP[2] + P(15,7)*SPP[5] + dt*(P(5,4) + P(0,4)*SF[4] + P(2,4)*SF[3] + P(3,4)*SF[5] - P(1,4)*SPP[0] - P(13,4)*SPP[8] + P(14,4)*SPP[2] + P(15,4)*SPP[5]);
        nextP(6,7) = P(6,7) + P(1,7)*SF[4] - P(2,7)*SF[5] + P(3,7)*SF[3] + P(0,7)*SPP[0] + P(13,7)*SPP[4] - P(14,7)*SPP[7] - P(15,7)*SPP[1] + dt*(P(6,4) + P(1,4)*SF[4] - P(2,4)*SF[5] + P(3,4)*SF[3] + P(0,4)*SPP[0] + P(13,4)*SPP[4] - P(14,4)*SPP[7] - P(15,4)*SPP[1]);
        nextP(7,7) = P(7,7) + P(4,7)*dt + dt*(P(7,4) + P(4,4)*dt);
        nextP(0,8) = P(0,8) + P(1,8)*SF[9] + P(2,8)*SF[11] + P(3,8)*SF[10] + P(10,8)*SF[14] + P(11,8)*SF[15] + P(12,8)*SPP[10] + dt*(P(0,5) + P(1,5)*SF[9] + P(2,5)*SF[11] + P(3,5)*SF[10] + P(10,5)*SF[14] + P(11,5)*SF[15] + P(12,5)*SPP[10]);
        nextP(1,8) = P(1,8) + P(0,8)*SF[8] + P(2,8)*SF[7] + P(3,8)*SF[11] - P(12,8)*SF[15] + P(11,8)*SPP[10] - (P(10,8)*q0)*0.5f + dt*(P(1,5) + P(0,5)*SF[8] + P(2,5)*SF[7] + P(3,5)*SF[11] - P(12,5)*SF[15] + P(11,5)*SPP[10] - (P(10,5)*q0)*0.5f);
        nextP(2,8) = P(2,8) + P(0,8)*SF[6] + P(1,8)*SF[10] + P(3,8)*SF[8] + P(12,8)*SF[14] - P(10,8)*SPP[10] - (P(11,8)*q0)*0.5f + dt*(P(2,5) + P(0,5)*SF[6] + P(1,5)*SF[10] + P(3,5)*SF[8] + P(12,5)*SF[14] - P(10,5)*SPP[10] - (P(11,5)*q0)*0.5f);
        nextP(3,8) = P(3,8) + P(0,8)*SF[7] + P(1,8)*SF[6] + P(2,8)*SF[9] + P(10,8)*SF[15] - P(11,8)*SF[14] - (P(12,8)*q0)*0.5f + dt*(P(3,5) + P(0,5)*SF[7] + P(1,5)*SF[6] + P(2,5)*SF[9] + P(10,5)*SF[15] - P(11,5)*SF[14] - (P(12,5)*q0)*0.5f);
        nextP(4,8) = P(4,8) + P(0,8)*SF[5] + P(1,8)*SF[3] - P(3,8)*SF[4] + P(2,8)*SPP[0] + P(13,8)*SPP[3] + P(14,8)*SPP[6] - P(15,8)*SPP[9] + dt*(P(4,5) + P(0,5)*SF[5] + P(1,5)*SF[3] - P(3,5)*SF[4] + P(2,5)*SPP[0] + P(13,5)*SPP[3] + P(14,5)*SPP[6] - P(15,5)*SPP[9]);
        nextP(5,8) = P(5,8) + P(0,8)*SF[4] + P(2,8)*SF[3] + P(3,8)*SF[5] - P(1,8)*SPP[0] - P(13,8)*SPP[8] + P(14,8)*SPP[2] + P(15,8)*SPP[5] + dt*(P(5,5) + P(0,5)*SF[4] + P(2,5)*SF[3] + P(3,5)*SF[5] - P(1,5)*SPP[0] - P(13,5)*SPP[8] + P(14,5)*SPP[2] + P(15,5)*SPP[5]);
        nextP(6,8) = P(6,8) + P(1,8)*SF[4] - P(2,8)*SF[5] + P(3,8)*SF[3] + P(0,8)*SPP[0] + P(13,8)*SPP[4] - P(14,8)*SPP[7] - P(15,8)*SPP[1] + dt*(P(6,5) + P(1,5)*SF[4] - P(2,5)*SF[5] + P(3,5)*SF[3] + P(0,5)*SPP[0] + P(13,5)*SPP[4] - P(14,5)*SPP[7] - P(15,5)*SPP[1]);
        nextP(7,8) = P(7,8) + P(4,8)*dt + dt*(P(7,5) + P(4,5)*dt);
        nextP(8,8) = P(8,8) + P(5,8)*dt + dt*(P(8,5) + P(5,5)*dt);
        nextP(0,9) = P(0,9) + P(1,9)*SF[9] + P(2,9)*SF[11] + P(3,9)*SF[10] + P(10,9)*SF[14] + P(11,9)*SF[15] + P(12,9)*SPP[10] + dt*(P(0,6) + P(1,6)*SF[9] + P(2,6)*SF[11] + P(3,6)*SF[10] + P(10,6)*SF[14] + P(11,6)*SF[15] + P(12,6)*SPP[10]);
        nextP(1,9) = P(1,9) + P(0,9)*SF[8] + P(2,9)*SF[7] + P(3,9)*SF[11] - P(12,9)*SF[15] + P(11,9)*SPP[10] - (P(10,9)*q0)*0.5f + dt*(P(1,6) + P(0,6)*SF[8] + P(2,6)*SF[7] + P(3,6)*SF[11] - P(12,6)*SF[15] + P(11,6)*SPP[10] - (P(10,6)*q0)*0.5f);
        nextP(2,9) = P(2,9) + P(0,9)*SF[6] + P(1,9)*SF[10] + P(3,9)*SF[8] + P(12,9)*SF[14] - P(10,9)*SPP[10] - (P(11,9)*q0)*0.5f + dt*(P(2,6) + P(0,6)*SF[6] + P(1,6)*SF[10] + P(3,6)*SF[8] + P(12,6)*SF[14] - P(10,6)*SPP[10] - (P(11,6)*q0)*0.5f);
        nextP(3,9) = P(3,9) + P(0,9)*SF[7] + P(1,9)*SF[6] + P(2,9)*SF[9] + P(10,9)*SF[15] - P(11,9)*SF[14] - (P(12,9)*q0)*0.5f + dt*(P(3,6) + P(0,6)*SF[7] + P(1,6)*SF[6] + P(2,6)*SF[9] + P(10,6)*SF[15] - P(11,6)*SF[14] - (P(12,6)*q0)*0.5f);
        nextP(4,9) = P(4,9) + P(0,9)*SF[5] + P(1,9)*SF[3] - P(3,9)*SF[4] + P(2,9)*SPP[0] + P(13,9)*SPP[3] + P(14,9)*SPP[6] - P(15,9)*SPP[9] + dt*(P(4,6) + P(0,6)*SF[5] + P(1,6)*SF[3] - P(3,6)*SF[4] + P(2,6)*SPP[0] + P(13,6)*SPP[3] + P(14,6)*SPP[6] - P(15,6)*SPP[9]);
        nextP(5,9) = P(5,9) + P(0,9)*SF[4] + P(2,9)*SF[3] + P(3,9)*SF[5] - P(1,9)*SPP[0] - P(13,9)*SPP[8] + P(14,9)*SPP[2] + P(15,9)*SPP[5] + dt*(P(5,6) + P(0,6)*SF[4] + P(2,6)*SF[3] + P(3,6)*SF[5] - P(1,6)*SPP[0] - P(13,6)*SPP[8] + P(14,6)*SPP[2] + P(15,6)*SPP[5]);
        nextP(6,9) = P(6,9) + P(1,9)*SF[4] - P(2,9)*SF[5] + P(3,9)*SF[3] + P(0,9)*SPP[0] + P(13,9)*SPP[4] - P(14,9)*SPP[7] - P(15,9)*SPP[1] + dt*(P(6,6) + P(1,6)*SF[4] - P(2,6)*SF[5] + P(3,6)*SF[3] + P(0,6)*SPP[0] + P(13,6)*SPP[4] - P(14,6)*SPP[7] - P(15,6)*SPP[1]);
        nextP(7,9) = P(7,9) + P(4,9)*dt + dt*(P(7,6) + P(4,6)*dt);
        nextP(8,9) = P(8,9) + P(5,9)*dt + dt*(P(8,6) + P(5,6)*dt);
        nextP(9,9) = P(9,9) + P(6,9)*dt + dt*(P(9,6) + P(6,6)*dt);
        nextP(0,10) = P(0,10) + P(1,10)*SF[9] + P(2,10)*SF[11] + P(3,10)*SF[10] + P(10,10)*SF[14] + P(11,10)*SF[15] + P(12,10)*SPP[10];
        nextP(1,10) = P(1,10) + P(0,10)*SF[8] + P(2,10)*SF[7] + P(3,10)*SF[11] - P(12,10)*SF[15] + P(11,10)*SPP[10] - (P(10,10)*q0)*0.5f;
        nextP(2,10) = P(2,10) + P(0,10)*SF[6] + P(1,10)*SF[10] + P(3,10)*SF[8] + P(12,10)*SF[14] - P(10,10)*SPP[10] - (P(11,10)*q0)*0.5f;
        nextP(3,10) = P(3,10) + P(0,10)*SF[7] + P(1,10)*SF[6] + P(2,10)*SF[9] + P(10,10)*SF[15] - P(11,10)*SF[14] - (P(12,10)*q0)*0.5f;
        nextP(4,10) = P(4,10) + P(0,10)*SF[5] + P(1,10)*SF[3] - P(3,10)*SF[4] + P(2,10)*SPP[0] + P(13,10)*SPP[3] + P(14,10)*SPP[6] - P(15,10)*SPP[9];
        nextP(5,10) = P(5,10) + P(0,10)*SF[4] + P(2,10)*SF[3] + P(3,10)*SF[5] - P(1,10)*SPP[0] - P(13,10)*SPP[8] + P(14,10)*SPP[2] + P(15,10)*SPP[5];
        nextP(6,10) = P(6,10) + P(1,10)*SF[4] - P(2,10)*SF[5] + P(3,10)*SF[3] + P(0,10)*SPP[0] + P(13,10)*SPP[4] - P(14,10)*SPP[7] - P(15,10)*SPP[1];
        nextP(7,10) = P(7,10) + P(4,10)*dt;
        nextP(8,10) = P(8,10) + P(5,10)*dt;
        nextP(9,10) = P(9,10) + P(6,10)*dt;
        nextP(10,10) = P(10,10);
        nextP(0,11) = P(0,11) + P(1,11)*SF[9] + P(2,11)*SF[11] + P(3,11)*SF[10] + P(10,11)*SF[14] + P(11,11)*SF[15] + P(12,11)*SPP[10];
        nextP(1,11) = P(1,11) + P(0,11)*SF[8] + P(2,11)*SF[7] + P(3,11)*SF[11] - P(12,11)*SF[15] + P(11,11)*SPP[10] - (P(10,11)*q0)*0.5f;
        nextP(2,11) = P(2,11) + P(0,11)*SF[6] + P(1,11)*SF[10] + P(3,11)*SF[8] + P(12,11)*SF[14] - P(10,11)*SPP[10] - (P(11,11)*q0)*0.5f;
        nextP(3,11) = P(3,11) + P(0,11)*SF[7] + P(1,11)*SF[6] + P(2,11)*SF[9] + P(10,11)*SF[15] - P(11,11)*SF[14] - (P(12,11)*q0)*0.5f;
        nextP(4,11) = P(4,11) + P(0,11)*SF[5] + P(1,11)*SF[3] - P(3,11)*SF[4] + P(2,11)*SPP[0] + P(13,11)*SPP[3] + P(14,11)*SPP[6] - P(15,11)*SPP[9];
        nextP(5,11) = P(5,11) + P(0,11)*SF[4] + P(2,11)*SF[3] + P(3,11)*SF[5] - P(1,11)*SPP[0] - P(13,11)*SPP[8] + P(14,11)*SPP[2] + P(15,11)*SPP[5];
        nextP(6,11) = P(6,11) + P(1,11)*SF[4] - P(2,11)*SF[5] + P(3,11)*SF[3] + P(0,11)*SPP[0] + P(13,11)*SPP[4] - P(14,11)*SPP[7] - P(15,11)*SPP[1];
        nextP(7,11) = P(7,11) + P(4,11)*dt;
        nextP(8,11) = P(8,11) + P(5,11)*dt;
        nextP(9,11) = P(9,11) + P(6,11)*dt;
        nextP(10,11) = P(10,11);
        nextP(11,11) = P(11,11);
        nextP(0,12) = P(0,12) + P(1,12)*SF[9] + P(2,12)*SF[11] + P(3,12)*SF[10] + P(10,12)*SF[14] + P(11,12)*SF[15] + P(12,12)*SPP[10];
        nextP(1,12) = P(1,12) + P(0,12)*SF[8] + P(2,12)*SF[7] + P(3,12)*SF[11] - P(12,12)*SF[15] + P(11,12)*SPP[10] - (P(10,12)*q0)*0.5f;
        nextP(2,12) = P(2,12) + P(0,12)*SF[6] + P(1,12)*SF[10] + P(3,12)*SF[8] + P(12,12)*SF[14] - P(10,12)*SPP[10] - (P(11,12)*q0)*0.5f;
        nextP(3,12) = P(3,12) + P(0,12)*SF[7] + P(1,12)*SF[6] + P(2,12)*SF[9] + P(10,12)*SF[15] - P(11,12)*SF[14] - (P(12,12)*q0)*0.5f;
        nextP(4,12) = P(4,12) + P(0,12)*SF[5] + P(1,12)*SF[3] - P(3,12)*SF[4] + P(2,12)*SPP[0] + P(13,12)*SPP[3] + P(14,12)*SPP[6] - P(15,12)*SPP[9];
        nextP(5,12) = P(5,12) + P(0,12)*SF[4] + P(2,12)*SF[3] + P(3,12)*SF[5] - P(1,12)*SPP[0] - P(13,12)*SPP[8] + P(14,12)*SPP[2] + P(15,12)*SPP[5];
        nextP(6,12) = P(6,12) + P(1,12)*SF[4] - P(2,12)*SF[5] + P(3,12)*SF[3] + P(0,12)*SPP[0] + P(13,12)*SPP[4] - P(14,12)*SPP[7] - P(15,12)*SPP[1];
        nextP(7,12) = P(7,12) + P(4,12)*dt;
        nextP(8,12) = P(8,12) + P(5,12)*dt;
        nextP(9,12) = P(9,12) + P(6,12)*dt;
        nextP(10,12) = P(10,12);
        nextP(11,12) = P(11,12);
        nextP(12,12) = P(12,12);

        for (unsigned i = 13; i <= 15; i++) {
            nextP(0,i) = P(0,i) + P(1,i)*SF[9] + P(2,i)*SF[11] + P(3,i)*SF[10] + P(10,i)*SF[14] + P(11,i)*SF[15] + P(12,i)*SPP[10];
            nextP(1,i) = P(1,i) + P(0,i)*SF[8] + P(2,i)*SF[7] + P(3,i)*SF[11] - P(12,i)*SF[15] + P(11,i)*SPP[10] - (P(10,i)*q0)*0.5f;
            nextP(2,i) = P(2,i) + P(0,i)*SF[6] + P(1,i)*SF[10] + P(3,i)*SF[8] + P(12,i)*SF[14] - P(10,i)*SPP[10] - (P(11,i)*q0)*0.5f;
            nextP(3,i) = P(3,i) + P(0,i)*SF[7] + P(1,i)*SF[6] + P(2,i)*SF[9] + P(10,i)*SF[15] - P(11,i)*SF[14] - (P(12,i)*q0)*0.5f;
            nextP(4,i) = P(4,i) + P(0,i)*SF[5] + P(1,i)*SF[3] - P(3,i)*SF[4] + P(2,i)*SPP[0] + P(13,i)*SPP[3] + P(14,i)*SPP[6] - P(15,i)*SPP[9];
            nextP(5,i) = P(5,i) + P(0,i)*SF[4] + P(2,i)*SF[3] + P(3,i)*SF[5] - P(1,i)*SPP[0] - P(13,i)*SPP[8] + P(14,i)*SPP[2] + P(15,i)*SPP[5];
            nextP(6,i) = P(6,i) + P(1,i)*SF[4] - P(2,i)*SF[5] + P(3,i)*SF[3] + P(0,i)*SPP[0] + P(13,i)*SPP[4] - P(14,i)*SPP[7] - P(15,i)*SPP[1];
            nextP(7,i) = P(7,i) + P(4,i)*dt;
            nextP(8,i) = P(8,i) + P(5,i)*dt;
            nextP(9,i) = P(9,i) + P(6,i)*dt;
            nextP(10,i) = P(10,i);
            nextP(11,i) = P(11,i);
            nextP(12,i) = P(12,i);
            nextP(13,i) = P(13,i);

            if (i > 13) {
                nextP(14,i) = P(14,i);
            }

            if (i > 14) {
                nextP(15,i) = P(15,i);
            }
        }

        nextP(0,16) = P(0,16) + P(1,16)*SF[9] + P(2,16)*SF[11] + P(3,16)*SF[10] + P(10,16)*SF[14] + P(11,16)*SF[15] + P(12,16)*SPP[10];
        nextP(1,16) = P(1,16) + P(0,16)*SF[8] + P(2,16)*SF[7] + P(3,16)*SF[11] - P(12,16)*SF[15] + P(11,16)*SPP[10] - (P(10,16)*q0)*0.5f;
        nextP(2,16) = P(2,16) + P(0,16)*SF[6] + P(1,16)*SF[10] + P(3,16)*SF[8] + P(12,16)*SF[14] - P(10,16)*SPP[10] - (P(11,16)*q0)*0.5f;
        nextP(3,16) = P(3,16) + P(0,16)*SF[7] + P(1,16)*SF[6] + P(2,16)*SF[9] + P(10,16)*SF[15] - P(11,16)*SF[14] - (P(12,16)*q0)*0.5f;
        nextP(4,16) = P(4,16) + P(0,16)*SF[5] + P(1,16)*SF[3] - P(3,16)*SF[4] + P(2,16)*SPP[0] + P(13,16)*SPP[3] + P(14,16)*SPP[6] - P(15,16)*SPP[9];
        nextP(5,16) = P(5,16) + P(0,16)*SF[4] + P(2,16)*SF[3] + P(3,16)*SF[5] - P(1,16)*SPP[0] - P(13,16)*SPP[8] + P(14,16)*SPP[2] + P(15,16)*SPP[5];
        nextP(6,16) = P(6,16) + P(1,16)*SF[4] - P(2,16)*SF[5] + P(3,16)*SF[3] + P(0,16)*SPP[0] + P(13,16)*SPP[4] - P(14,16)*SPP[7] - P(15,16)*SPP[1];
        nextP(7,16) = P(7,16) + P(4,16)*dt;
        nextP(8,16) = P(8,16) + P(5,16)*dt;
        nextP(9,16) = P(9,16) + P(6,16)*dt;
        nextP(10,16) = P(10,16);
        nextP(11,16) = P(11,16);
        nextP(12,16) = P(12,16);
        nextP(13,16) = P(13,16);
        nextP(14,16) = P(14,16);
        nextP(15,16) = P(15,16);
        nextP(16,16) = P(16,16);
        nextP(0,17) = P(0,17) + P(1,17)*SF[9] + P(2,17)*SF[11] + P(3,17)*SF[10] + P(10,17)*SF[14] + P(11,17)*SF[15] + P(12,17)*SPP[10];
        nextP(1,17) = P(1,17) + P(0,17)*SF[8] + P(2,17)*SF[7] + P(3,17)*SF[11] - P(12,17)*SF[15] + P(11,17)*SPP[10] - (P(10,17)*q0)*0.5f;
        nextP(2,17) = P(2,17) + P(0,17)*SF[6] + P(1,17)*SF[10] + P(3,17)*SF[8] + P(12,17)*SF[14] - P(10,17)*SPP[10] - (P(11,17)*q0)*0.5f;
        nextP(3,17) = P(3,17) + P(0,17)*SF[7] + P(1,17)*SF[6] + P(2,17)*SF[9] + P(10,17)*SF[15] - P(11,17)*SF[14] - (P(12,17)*q0)*0.5f;
        nextP(4,17) = P(4,17) + P(0,17)*SF[5] + P(1,17)*SF[3] - P(3,17)*SF[4] + P(2,17)*SPP[0] + P(13,17)*SPP[3] + P(14,17)*SPP[6] - P(15,17)*SPP[9];
        nextP(5,17) = P(5,17) + P(0,17)*SF[4] + P(2,17)*SF[3] + P(3,17)*SF[5] - P(1,17)*SPP[0] - P(13,17)*SPP[8] + P(14,17)*SPP[2] + P(15,17)*SPP[5];
        nextP(6,17) = P(6,17) + P(1,17)*SF[4] - P(2,17)*SF[5] + P(3,17)*SF[3] + P(0,17)*SPP[0] + P(13,17)*SPP[4] - P(14,17)*SPP[7] - P(15,17)*SPP[1];
        nextP(7,17) = P(7,17) + P(4,17)*dt;
        nextP(8,17) = P(8,17) + P(5,17)*dt;
        nextP(9,17) = P(9,17) + P(6,17)*dt;
        nextP(10,17) = P(10,17);
        nextP(11,17) = P(11,17);
        nextP(12,17) = P(12,17);
        nextP(13,17) = P(13,17);
        nextP(14,17) = P(14,17);
        nextP(15,17) = P(15,17);
        nextP(16,17) = P(16,17);
        nextP(17,17) = P(17,17);
        nextP(0,18) = P(0,18) + P(1,18)*SF[9] + P(2,18)*SF[11] + P(3,18)*SF[10] + P(10,18)*SF[14] + P(11,18)*SF[15] + P(12,18)*SPP[10];
        nextP(1,18) = P(1,18) + P(0,18)*SF[8] + P(2,18)*SF[7] + P(3,18)*SF[11] - P(12,18)*SF[15] + P(11,18)*SPP[10] - (P(10,18)*q0)*0.5f;
        nextP(2,18) = P(2,18) + P(0,18)*SF[6] + P(1,18)*SF[10] + P(3,18)*SF[8] + P(12,18)*SF[14] - P(10,18)*SPP[10] - (P(11,18)*q0)*0.5f;
        nextP(3,18) = P(3,18) + P(0,18)*SF[7] + P(1,18)*SF[6] + P(2,18)*SF[9] + P(10,18)*SF[15] - P(11,18)*SF[14] - (P(12,18)*q0)*0.5f;
        nextP(4,18) = P(4,18) + P(0,18)*SF[5] + P(1,18)*SF[3] - P(3,18)*SF[4] + P(2,18)*SPP[0] + P(13,18)*SPP[3] + P(14,18)*SPP[6] - P(15,18)*SPP[9];
        nextP(5,18) = P(5,18) + P(0,18)*SF[4] + P(2,18)*SF[3] + P(3,18)*SF[5] - P(1,18)*SPP[0] - P(13,18)*SPP[8] + P(14,18)*SPP[2] + P(15,18)*SPP[5];
        nextP(6,18) = P(6,18) + P(1,18)*SF[4] - P(2,18)*SF[5] + P(3,18)*SF[3] + P(0,18)*SPP[0] + P(13,18)*SPP[4] - P(14,18)*SPP[7] - P(15,18)*SPP[1];
        nextP(7,18) = P(7,18) + P(4,18)*dt;
        nextP(8,18) = P(8,18) + P(5,18)*dt;
        nextP(9,18) = P(9,18) + P(6,18)*dt;
        nextP(10,18) = P(10,18);
        nextP(11,18) = P(11,18);
        nextP(12,18) = P(12,18);
        nextP(13,18) = P(13,18);
        nextP(14,18) = P(14,18);
        nextP(15,18) = P(15,18);
        nextP(16,18) = P(16,18);
        nextP(17,18) = P(17,18);
        nextP(18,18) = P(18,18);
        nextP(0,19) = P(0,19) + P(1,19)*SF[9] + P(2,19)*SF[11] + P(3,19)*SF[10] + P(10,19)*SF[14] + P(11,19)*SF[15] + P(12,19)*SPP[10];
        nextP(1,19) = P(1,19) + P(0,19)*SF[8] + P(2,19)*SF[7] + P(3,19)*SF[11] - P(12,19)*SF[15] + P(11,19)*SPP[10] - (P(10,19)*q0)*0.5f;
        nextP(2,19) = P(2,19) + P(0,19)*SF[6] + P(1,19)*SF[10] + P(3,19)*SF[8] + P(12,19)*SF[14] - P(10,19)*SPP[10] - (P(11,19)*q0)*0.5f;
        nextP(3,19) = P(3,19) + P(0,19)*SF[7] + P(1,19)*SF[6] + P(2,19)*SF[9] + P(10,19)*SF[15] - P(11,19)*SF[14] - (P(12,19)*q0)*0.5f;
        nextP(4,19) = P(4,19) + P(0,19)*SF[5] + P(1,19)*SF[3] - P(3,19)*SF[4] + P(2,19)*SPP[0] + P(13,19)*SPP[3] + P(14,19)*SPP[6] - P(15,19)*SPP[9];
        nextP(5,19) = P(5,19) + P(0,19)*SF[4] + P(2,19)*SF[3] + P(3,19)*SF[5] - P(1,19)*SPP[0] - P(13,19)*SPP[8] + P(14,19)*SPP[2] + P(15,19)*SPP[5];
        nextP(6,19) = P(6,19) + P(1,19)*SF[4] - P(2,19)*SF[5] + P(3,19)*SF[3] + P(0,19)*SPP[0] + P(13,19)*SPP[4] - P(14,19)*SPP[7] - P(15,19)*SPP[1];
        nextP(7,19) = P(7,19) + P(4,19)*dt;
        nextP(8,19) = P(8,19) + P(5,19)*dt;
        nextP(9,19) = P(9,19) + P(6,19)*dt;
        nextP(10,19) = P(10,19);
        nextP(11,19) = P(11,19);
        nextP(12,19) = P(12,19);
        nextP(13,19) = P(13,19);
        nextP(14,19) = P(14,19);
        nextP(15,19) = P(15,19);
        nextP(16,19) = P(16,19);
        nextP(17,19) = P(17,19);
        nextP(18,19) = P(18,19);
        nextP(19,19) = P(19,19);
        nextP(0,20) = P(0,20) + P(1,20)*SF[9] + P(2,20)*SF[11] + P(3,20)*SF[10] + P(10,20)*SF[14] + P(11,20)*SF[15] + P(12,20)*SPP[10];
        nextP(1,20) = P(1,20) + P(0,20)*SF[8] + P(2,20)*SF[7] + P(3,20)*SF[11] - P(12,20)*SF[15] + P(11,20)*SPP[10] - (P(10,20)*q0)*0.5f;
        nextP(2,20) = P(2,20) + P(0,20)*SF[6] + P(1,20)*SF[10] + P(3,20)*SF[8] + P(12,20)*SF[14] - P(10,20)*SPP[10] - (P(11,20)*q0)*0.5f;
        nextP(3,20) = P(3,20) + P(0,20)*SF[7] + P(1,20)*SF[6] + P(2,20)*SF[9] + P(10,20)*SF[15] - P(11,20)*SF[14] - (P(12,20)*q0)*0.5f;
        nextP(4,20) = P(4,20) + P(0,20)*SF[5] + P(1,20)*SF[3] - P(3,20)*SF[4] + P(2,20)*SPP[0] + P(13,20)*SPP[3] + P(14,20)*SPP[6] - P(15,20)*SPP[9];
        nextP(5,20) = P(5,20) + P(0,20)*SF[4] + P(2,20)*SF[3] + P(3,20)*SF[5] - P(1,20)*SPP[0] - P(13,20)*SPP[8] + P(14,20)*SPP[2] + P(15,20)*SPP[5];
        nextP(6,20) = P(6,20) + P(1,20)*SF[4] - P(2,20)*SF[5] + P(3,20)*SF[3] + P(0,20)*SPP[0] + P(13,20)*SPP[4] - P(14,20)*SPP[7] - P(15,20)*SPP[1];
        nextP(7,20) = P(7,20) + P(4,20)*dt;
        nextP(8,20) = P(8,20) + P(5,20)*dt;
        nextP(9,20) = P(9,20) + P(6,20)*dt;
        nextP(10,20) = P(10,20);
        nextP(11,20) = P(11,20);
        nextP(12,20) = P(12,20);
        nextP(13,20) = P(13,20);
        nextP(14,20) = P(14,20);
        nextP(15,20) = P(15,20);
        nextP(16,20) = P(16,20);
        nextP(17,20) = P(17,20);
        nextP(18,20) = P(18,20);
        nextP(19,20) = P(19,20);
        nextP(20,20) = P(20,20);
        nextP(0,21) = P(0,21) + P(1,21)*SF[9] + P(2,21)*SF[11] + P(3,21)*SF[10] + P(10,21)*SF[14] + P(11,21)*SF[15] + P(12,21)*SPP[10];
        nextP(1,21) = P(1,21) + P(0,21)*SF[8] + P(2,21)*SF[7] + P(3,21)*SF[11] - P(12,21)*SF[15] + P(11,21)*SPP[10] - (P(10,21)*q0)*0.5f;
        nextP(2,21) = P(2,21) + P(0,21)*SF[6] + P(1,21)*SF[10] + P(3,21)*SF[8] + P(12,21)*SF[14] - P(10,21)*SPP[10] - (P(11,21)*q0)*0.5f;
        nextP(3,21) = P(3,21) + P(0,21)*SF[7] + P(1,21)*SF[6] + P(2,21)*SF[9] + P(10,21)*SF[15] - P(11,21)*SF[14] - (P(12,21)*q0)*0.5f;
        nextP(4,21) = P(4,21) + P(0,21)*SF[5] + P(1,21)*SF[3] - P(3,21)*SF[4] + P(2,21)*SPP[0] + P(13,21)*SPP[3] + P(14,21)*SPP[6] - P(15,21)*SPP[9];
        nextP(5,21) = P(5,21) + P(0,21)*SF[4] + P(2,21)*SF[3] + P(3,21)*SF[5] - P(1,21)*SPP[0] - P(13,21)*SPP[8] + P(14,21)*SPP[2] + P(15,21)*SPP[5];
        nextP(6,21) = P(6,21) + P(1,21)*SF[4] - P(2,21)*SF[5] + P(3,21)*SF[3] + P(0,21)*SPP[0] + P(13,21)*SPP[4] - P(14,21)*SPP[7] - P(15,21)*SPP[1];
        nextP(7,21) = P(7,21) + P(4,21)*dt;
        nextP(8,21) = P(8,21) + P(5,21)*dt;
        nextP(9,21) = P(9,21) + P(6,21)*dt;
        nextP(10,21) = P(10,21);
        nextP(11,21) = P(11,21);
        nextP(12,21) = P(12,21);
        nextP(13,21) = P(13,21);
        nextP(14,21) = P(14,21);
        nextP(15,21) = P(15,21);
        nextP(16,21) = P(16,21);
        nextP(17,21) = P(17,21);
        nextP(18,21) = P(18,21);
        nextP(19,21) = P(19,21);
        nextP(20,21) = P(20,21);
        nextP(21,21) = P(21,21);
        nextP(0,22) = P(0,22) + P(1,22)*SF[9] + P(2,22)*SF[11] + P(3,22)*SF[10] + P(10,22)*SF[14] + P(11,22)*SF[15] + P(12,22)*SPP[10];
        nextP(1,22) = P(1,22) + P(0,22)*SF[8] + P(2,22)*SF[7] + P(3,22)*SF[11] - P(12,22)*SF[15] + P(11,22)*SPP[10] - (P(10,22)*q0)*0.5f;
        nextP(2,22) = P(2,22) + P(0,22)*SF[6] + P(1,22)*SF[10] + P(3,22)*SF[8] + P(12,22)*SF[14] - P(10,22)*SPP[10] - (P(11,22)*q0)*0.5f;
        nextP(3,22) = P(3,22) + P(0,22)*SF[7] + P(1,22)*SF[6] + P(2,22)*SF[9] + P(10,22)*SF[15] - P(11,22)*SF[14] - (P(12,22)*q0)*0.5f;
        nextP(4,22) = P(4,22) + P(0,22)*SF[5] + P(1,22)*SF[3] - P(3,22)*SF[4] + P(2,22)*SPP[0] + P(13,22)*SPP[3] + P(14,22)*SPP[6] - P(15,22)*SPP[9];
        nextP(5,22) = P(5,22) + P(0,22)*SF[4] + P(2,22)*SF[3] + P(3,22)*SF[5] - P(1,22)*SPP[0] - P(13,22)*SPP[8] + P(14,22)*SPP[2] + P(15,22)*SPP[5];
        nextP(6,22) = P(6,22) + P(1,22)*SF[4] - P(2,22)*SF[5] + P(3,22)*SF[3] + P(0,22)*SPP[0] + P(13,22)*SPP[4] - P(14,22)*SPP[7] - P(15,22)*SPP[1];
        nextP(7,22) = P(7,22) + P(4,22)*dt;
        nextP(8,22) = P(8,22) + P(5,22)*dt;
        nextP(9,22) = P(9,22) + P(6,22)*dt;
        nextP(10,22) = P(10,22);
        nextP(11,22) = P(11,22);
        nextP(12,22) = P(12,22);
        nextP(13,22) = P(13,22);
        nextP(14,22) = P(14,22);
        nextP(15,22) = P(15,22);
        nextP(16,22) = P(16,22);
        nextP(17,22) = P(17,22);
        nextP(18,22) = P(18,22);
        nextP(19,22) = P(19,22);
        nextP(20,22) = P(20,22);
        nextP(21,22) = P(21,22);
        nextP(22,22) = P(22,22);
        nextP(0,23) = P(0,23) + P(1,23)*SF[9] + P(2,23)*SF[11] + P(3,23)*SF[10] + P(10,23)*SF[14] + P(11,23)*SF[15] + P(12,23)*SPP[10];
        nextP(1,23) = P(1,23) + P(0,23)*SF[8] + P(2,23)*SF[7] + P(3,23)*SF[11] - P(12,23)*SF[15] + P(11,23)*SPP[10] - (P(10,23)*q0)*0.5f;
        nextP(2,23) = P(2,23) + P(0,23)*SF[6] + P(1,23)*SF[10] + P(3,23)*SF[8] + P(12,23)*SF[14] - P(10,23)*SPP[10] - (P(11,23)*q0)*0.5f;
        nextP(3,23) = P(3,23) + P(0,23)*SF[7] + P(1,23)*SF[6] + P(2,23)*SF[9] + P(10,23)*SF[15] - P(11,23)*SF[14] - (P(12,23)*q0)*0.5f;
        nextP(4,23) = P(4,23) + P(0,23)*SF[5] + P(1,23)*SF[3] - P(3,23)*SF[4] + P(2,23)*SPP[0] + P(13,23)*SPP[3] + P(14,23)*SPP[6] - P(15,23)*SPP[9];
        nextP(5,23) = P(5,23) + P(0,23)*SF[4] + P(2,23)*SF[3] + P(3,23)*SF[5] - P(1,23)*SPP[0] - P(13,23)*SPP[8] + P(14,23)*SPP[2] + P(15,23)*SPP[5];
        nextP(6,23) = P(6,23) + P(1,23)*SF[4] - P(2,23)*SF[5] + P(3,23)*SF[3] + P(0,23)*SPP[0] + P(13,23)*SPP[4] - P(14,23)*SPP[7] - P(15,23)*SPP[1];
        nextP(7,23) = P(7,23) + P(4,23)*dt;
        nextP(8,23) = P(8,23) + P(5,23)*dt;
        nextP(9,23) = P(9,23) + P(6,23)*dt;
        nextP(10,23) = P(10,23);
        nextP(11,23) = P(11,23);
        nextP(12,23) = P(12,23);
        nextP(13,23) = P(13,23);
        nextP(14,23) = P(14,23);
        nextP(15,23) = P(15,23);
        nextP(16,23) = P(16,23);
        nextP(17,23) = P(17,23);
        nextP(18,23) = P(18,23);
        nextP(19,23) = P(19,23);
        nextP(20,23) = P(20,23);
        nextP(21,23) = P(21,23);
        nextP(22,23) = P(22,23);
        nextP(23,23) = P(23,23);

        // capture largest difference
        float max_diff_fraction = 0.0f;
        int max_row, max_col;
        float max_old, max_new;

        for (int col=0; col<=23; col++) {
            for (int row=0; row<=col; row++) {
                float diff_fraction = fabsf(nextP_sympy(row,col)-nextP(row,col)) / fabsf(nextP(row,col));
                if (diff_fraction > max_diff_fraction) {
                    max_diff_fraction = diff_fraction;
                    max_row = row;
                    max_col = col;
                    max_old = nextP(row,col);
                    max_new = nextP_sympy(row,col);
                }
            }
        }
        if (max_diff_fraction > 5E-5f) {
            printf("Fail: Covariance Prediction max diff fraction = %e , old = %e , new = %e , location index = %i,%i\n",max_diff_fraction, max_old, max_new, max_row, max_col);
        } else {
            printf("Pass: Covariance Prediction max diff fraction = %e , old = %e , new = %e , location index = %i,%i\n",max_diff_fraction, max_old, max_new, max_row, max_col);
        }
    }

    return 0;
}
