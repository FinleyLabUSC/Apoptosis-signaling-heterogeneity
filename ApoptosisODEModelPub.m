function dydt = odeModel(t,y,p)
% inputs: t - the timepoint
%         y - concentration of each species at timepoint t
%         p - parameters

% species
C8 = y(1);
C8a = y(2);
C3 = y(3);
C3a = y(4);
IAP = y(5);
C3a_IAP = y(6);
BAR = y(7);
C8a_BAR = y(8);
receptor = y(9);
complex = y(10);
ligand = y(11);

% parameters
k1f = p(1);
k1r = p(2);
k2f = p(3);
k2r = p(4);
k3f = p(5);
k3r = p(6);
k4f = p(7);
k4r = p(8);
k5f = p(9);
k5r = p(10);
k6f = p(11);
k6r = p(12);
k7f = p(13);
k7r = p(14);
k8f = p(15);
k8r = p(16);
k9f = p(17);
k9r = p(18);
k10f = p(19);
k10r = p(20);
k11f = p(21);
k11r = p(22);
k12f = p(23);
k12r = p(24);
k13f = p(25);
k13r = p(26);
k14f = p(27);
k14r = p(28);
ku = p(29);
kreceptor_production = p(30);
kreceptor_degredation = p(31);
kcomplex_internalization = p(32);


u = ku*complex;
%r_d = 0; %Placeholder for Parameter for degredation of receptor 
%r_p = 0; %Placeholder for Parameter for production of receptor

%Equation System
%v0 = r_p[receptor];
v1 = k1f*C8a*C3;
v2 = k2f*C3a*C8;
v3 = k3f*C3a*IAP-k3r*C3a_IAP;
v4 = k4f*C3a*IAP;
v5 = k5f*C8a;
v6 = k6f*C3a;
v7 = k7f*C3a_IAP;
v8 = k8f*IAP-k8r;
v9 = k9f*C8-k9r;
v10 = k10f*C3 - k10r;
v11 = k11f*C8a*BAR-k11r*C8a_BAR;
v12 = k12f*BAR - k12r;
v13 = k13f*C8a_BAR;
v14 = k14f*ligand*receptor - k14r*complex;
v15 = ku*complex*C8; %u*C8;
v16 = kreceptor_production - kreceptor_degredation*receptor;
v17 = kcomplex_internalization*complex;


% reactions
% receptorBinding = k1f*ligand*receptor - k1r*complex;

dydt(1,1)= -v2-v9-v15; %C8
dydt(2,1)= v2+v15-v5-v11; %C8a
dydt(3,1)= -v1-v10; %C3
dydt(4,1)= v1-v3-v6; %C3a
dydt(5,1)= -v3-v4-v8; %IAP
dydt(6,1)= v3-v7; %C3a_IAP
dydt(7,1)= -v11-v12; %BAR
dydt(8,1)= v11 - v13; %C8aBAR
dydt(9,1) = v16 - v14; % receptor
dydt(10,1) = v14 - v17; % complex
dydt(11,1) = -v14; % ligand
%dydt(12,1) = v0*v15 %production of 

end

