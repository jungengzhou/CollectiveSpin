clc;
clear;
close all;

N = 10;
sys = DickeTools(N);
theta_SCS = pi/2;
phi_SCS = 0;
SCS = sys.SCS(theta_SCS,phi_SCS);
[Q,h] = bloch(SCS);
%[Q,h] = bloch(SCS,"x");


chit = 0.1*pi;
U = sys.OAT(chit,"z");
SSS = U*SCS;
[Q1,h1] = bloch(SSS);
