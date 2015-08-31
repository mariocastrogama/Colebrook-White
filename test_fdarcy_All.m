% test friction factor by Colebrook-White
% Implicit Newton-Raphson vs explicit formulations
% solves one specific example
% 
% Mario Castro Gama
% MSc Hydroinformatics
% 2014.05.25
%
clc;
clear;
format long g


ks = 1.5e-6; % [m]   
D  = 0.1522;    % [m]
Q  = 42/1000;    % [m^3/s]
v  = 1.14e-6; % [m^2/s]
Re=numre(Q,D,v); % [adim]

% display the results for each formulation
disp('Results');
disp(['  f Newton-Raphson : ',num2str(fdarcynewton(ks,D,Q,v))]);
disp(['            f Avci : ',num2str(fdarcyavci(ks,D,Q,v))]);
disp(['            f Barr : ',num2str(fdarcybarr(ks,D,Q,v))]);
disp(['          f Brkic1 : ',num2str(fdarcybrkic1(ks,D,Q,v))]);
disp(['          f Brkic2 : ',num2str(fdarcybrkic2(ks,D,Q,v))]);
disp(['         f Buzelli : ',num2str(fdarcybuzzelli(ks,D,Q,v))]);
disp(['            f Chen : ',num2str(fdarcychen(ks,D,Q,v))]);
disp(['       f Churchill : ',num2str(fdarcychurchill(ks,D,Q,v))]);
disp(['        f Davidson : ',num2str(fdarcydavidson(ks,D,Q,v))]);
disp(['             f Eck : ',num2str(fdarcyeck(ks,D,Q,v))]);
disp(['            f Fang : ',num2str(fdarcyfang(ks,D,Q,v))]);
disp(['         f Goudar1 : ',num2str(fdarcygoudar1(ks,D,Q,v))]);
disp(['         f Goudar2 : ',num2str(fdarcygoudar2(ks,D,Q,v))]);
disp(['         f Haaland : ',num2str(fdarcyhaaland(ks,D,Q,v))]);
disp(['       f Manadilli : ',num2str(fdarcymanadilli(ks,D,Q,v))]);
disp(['           f Moody : ',num2str(fdarcymoody(ks,D,Q,v))]);
disp(['   f Papaevangelou : ',num2str(fdarcypapaevangelouetal(ks,D,Q,v))]);
disp(['             f Rao : ',num2str(fdarcyrao(ks,D,Q,v))]);
disp(['           f Romeo : ',num2str(fdarcyromeoetal(ks,D,Q,v))]);
disp(['           f Rough : ',num2str(fdarcyrough(ks,D))]);
disp(['           f Round : ',num2str(fdarcyround(ks,D,Q,v))]);
disp(['      f Serghides1 : ',num2str(fdarcyserghides1(ks,D,Q,v))]);
disp(['      f Serghides2 : ',num2str(fdarcyserghides2(ks,D,Q,v))]);
disp(['        f Slamashi : ',num2str(fdarcyslamashietal(ks,D,Q,v))]);
disp(['          f Sonnad : ',num2str(fdarcysonnad(ks,D,Q,v))]);
disp(['     f Swamee-Jain : ',num2str(fdarcyswameejain(ks,D,Q,v))]);
disp(['   f Swamee-Rathie : ',num2str(fdarcyswameerathie(ks,D,Q,v))]);
disp(['             f Vak : ',num2str(fdarcyvak(ks,D,Q,v))]);
disp(['            f Wood : ',num2str(fdarcywood(ks,D,Q,v))]);
disp(['        f Zigrang1 : ',num2str(fdarcyzigrang1(ks,D,Q,v))]);
disp(['        f Zigrang2 : ',num2str(fdarcyzigrang1(ks,D,Q,v))]);
