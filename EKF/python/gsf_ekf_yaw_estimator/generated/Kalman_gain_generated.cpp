// Equations for NE velocity Kalman gain
const float SK0 = powf(P01, 2);
const float SK1 = P11 + velObsVar;
const float SK2 = P00 + velObsVar;
const float SK3 = 1.0F/(SK0 - SK1*SK2);
const float SK4 = -P01*SK3*velObsVar;


K(0,0) = SK3*(-P00*SK1 + SK0);
K(1,0) = SK4;
K(2,0) = SK3*(P01*P12 - P02*SK1);
K(0,1) = SK4;
K(1,1) = SK3*(-P11*SK2 + SK0);
K(2,1) = SK3*(P01*P02 - P12*SK2);


