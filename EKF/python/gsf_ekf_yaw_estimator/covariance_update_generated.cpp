// Equations for covariance matrix update
const float SP0 = P11 + velObsVar;
const float SP1 = powf(P01, 2);
const float SP2 = -SP1;
const float SP3 = P00 + velObsVar;
const float SP4 = SP0*SP3;
const float SP5 = SP2 + SP4;
const float SP6 = 1.0F/SP5;
const float SP7 = SP6*velObsVar;
const float SP8 = (-P00*SP0 + SP1)/(SP1 - SP4);
const float SP9 = SP0*SP7 + SP8;
const float SP10 = SP1*velObsVar;
const float SP11 = SP10*SP6;
const float SP12 = SP11 + SP3*SP8;
const float SP13 = P11*SP3;
const float SP14 = SP1 - SP13;
const float SP15 = SP6*SP9;
const float SP16 = P01*P02 - P12*SP3;
const float SP17 = P01*SP16;
const float SP18 = P01*P12 - P02*SP0;
const float SP19 = powf(SP5, -2);
const float SP20 = SP19*(-SP0*SP14 + SP10);
const float SP21 = SP13 + SP2 + SP3*velObsVar;
const float SP22 = SP6*SP8;
const float SP23 = SP19*SP21;
const float SP24 = P01*SP18;
const float SP25 = SP19*(SP0*SP16 + SP24);
const float SP26 = P01*velObsVar;
const float SP27 = SP17 + SP18*SP3;
const float SP28 = SP19*SP27;


_ekf_gsf[model_index].P(0,0) = P00 - SP11*SP9 - SP12*SP8;
_ekf_gsf[model_index].P(0,1) = P01*(-SP12*SP7 + SP14*SP15 + 1);
_ekf_gsf[model_index].P(1,1) = P11 - SP10*SP23 + SP14*SP20;
_ekf_gsf[model_index].P(0,2) = P02 + SP12*SP18*SP6 + SP15*SP17;
_ekf_gsf[model_index].P(1,2) = P12 + SP16*SP20 + SP23*SP24;
_ekf_gsf[model_index].P(2,2) = P22 - SP16*SP25 - SP18*SP28;


