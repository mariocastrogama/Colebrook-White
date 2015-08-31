% test friction Colebrook-White
% several equations for distributed points in a meshgrid
%
% Shows the Absolute Error between the Implicit approximation (Newton-Raphson)
% and several equations gathered from internet.
% 
% The range of validity of the calculation is:
% ks/D = [1.0e-6,1.0e-01];
% Re   = [1.0e+04,1.0e+08];
% 
% Mario Castro Gama
% m.castrogama@unesc0-ihe.org
% MSc Hydroinformatics
% 2012.11.20
%
% Last Update:
% 2014.07.03
clc;
clear;
close all;
fclose all;
format long g;

% number of equations used for evaluation of friction factor
neqs = 31;

% Fixed values of Diameter [m] and Discharge [m^3/s]
D  = 1.0;    % [m]
Q  = 1.0;    % [m^3/s]
% Range of values of log10(ksi) is [1.0e-6, 1.0e-1] 
ksmin  = log10(1.0e-6);
ksmax  = log10(1.0e-1);
ksi    = ksmin:(ksmax-ksmin)/100:ksmax;
% Range of values of log10(vj)
c1     = 4/pi;
vmin   = log10(c1/0.4e4);
vmax   = log10(c1/1.0e8);
vi     = vmin:(vmax-vmin)/100:vmax;

% Empty matrices required to store values 
Re      = zeros([length(ksi),length(vi)]);
f      = zeros([length(ksi),length(vi),neqs]);
fobj = zeros([length(ksi),length(vi),neqs-1]);

% Main loop of calculation in every formula for friction of Darcy Weisbach
for ii=1:length(ksi);
  for jj=1:length(vi);
    ks = 10^ksi(ii); % [m]   
    v=10^vi(jj); % [m^2/s]
    Re(ii,jj)=numre(Q,D,v); % [adim]
    kk = 1;
    f(ii,jj, kk) = fdarcynewton(ks,D,Q,v);     kk=kk+1; if (ii+jj==2); list_names={'Newton-Raphson'}; end
    f(ii,jj, kk) = fdarcyavci(ks,D,Q,v);       kk=kk+1; if (ii+jj==2); list_names=cat(2,list_names,{'AVK'}); end
    f(ii,jj, kk) = fdarcybarr(ks,D,Q,v);       kk=kk+1; if (ii+jj==2); list_names=cat(2,list_names,{'BAR'}); end
    f(ii,jj, kk) = fdarcybrkic1(ks,D,Q,v);     kk=kk+1; if (ii+jj==2); list_names=cat(2,list_names,{'BRK1'}); end
    f(ii,jj, kk) = fdarcybrkic2(ks,D,Q,v);     kk=kk+1; if (ii+jj==2); list_names=cat(2,list_names,{'BRK2'}); end
    f(ii,jj, kk) = fdarcybuzzelli(ks,D,Q,v);   kk=kk+1; if (ii+jj==2); list_names=cat(2,list_names,{'BUZ'}); end
    f(ii,jj, kk) = fdarcychen(ks,D,Q,v);       kk=kk+1; if (ii+jj==2); list_names=cat(2,list_names,{'CHN'}); end
    f(ii,jj, kk) = fdarcychurchill(ks,D,Q,v);  kk=kk+1; if (ii+jj==2); list_names=cat(2,list_names,{'CHR'}); end
    f(ii,jj, kk) = fdarcydavidson(ks,D,Q,v);        kk=kk+1; if (ii+jj==2); list_names=cat(2,list_names,{'DAV'}); end
    f(ii,jj, kk) = fdarcyeck(ks,D,Q,v);        kk=kk+1; if (ii+jj==2); list_names=cat(2,list_names,{'ECK'}); end
    f(ii,jj, kk) = fdarcyfang(ks,D,Q,v);       kk=kk+1; if (ii+jj==2); list_names=cat(2,list_names,{'FNG'}); end
    f(ii,jj, kk) = fdarcygoudar1(ks,D,Q,v);    kk=kk+1; if (ii+jj==2); list_names=cat(2,list_names,{'GAS1'}); end
    f(ii,jj, kk) = fdarcygoudar2(ks,D,Q,v);    kk=kk+1; if (ii+jj==2); list_names=cat(2,list_names,{'GAS2'}); end
    f(ii,jj, kk) = fdarcyhaaland(ks,D,Q,v);    kk=kk+1; if (ii+jj==2); list_names=cat(2,list_names,{'HAL'}); end
    f(ii,jj, kk) = fdarcymanadilli(ks,D,Q,v);  kk=kk+1; if (ii+jj==2); list_names=cat(2,list_names,{'MAN'}); end
    f(ii,jj, kk) = fdarcymoody(ks,D,Q,v);      kk=kk+1; if (ii+jj==2); list_names=cat(2,list_names,{'MOO'}); end
    f(ii,jj, kk) = fdarcypapaevangelouetal(ks,D,Q,v);  kk=kk+1; if (ii+jj==2); list_names=cat(2,list_names,{'PET'}); end
    f(ii,jj, kk) = fdarcyrao(ks,D,Q,v);        kk=kk+1; if (ii+jj==2); list_names=cat(2,list_names,{'RAK'}); end
    f(ii,jj, kk) = fdarcyromeoetal(ks,D,Q,v);  kk=kk+1; if (ii+jj==2); list_names=cat(2,list_names,{'RRM'}); end
    f(ii,jj, kk) = fdarcyrough(ks,D);  kk=kk+1; if (ii+jj==2); list_names=cat(2,list_names,{'RGH'}); end
    f(ii,jj, kk) = fdarcyround(ks,D,Q,v);      kk=kk+1; if (ii+jj==2); list_names=cat(2,list_names,{'ROU'}); end
    f(ii,jj, kk) = fdarcyserghides1(ks,D,Q,v); kk=kk+1; if (ii+jj==2); list_names=cat(2,list_names,{'SER1'}); end
    f(ii,jj, kk) = fdarcyserghides2(ks,D,Q,v); kk=kk+1; if (ii+jj==2); list_names=cat(2,list_names,{'SER2'}); end
    f(ii,jj, kk) = fdarcyslamashietal(ks,D,Q,v);        kk=kk+1; if (ii+jj==2); list_names=cat(2,list_names,{'SKG'}); end
    f(ii,jj, kk) = fdarcysonnad(ks,D,Q,v);     kk=kk+1; if (ii+jj==2); list_names=cat(2,list_names,{'SAG'}); end
    f(ii,jj, kk) = fdarcyswameejain(ks,D,Q,v);      kk=kk+1; if (ii+jj==2); list_names=cat(2,list_names,{'SAJ'}); end
    f(ii,jj, kk) = fdarcyswameerathie(ks,D,Q,v);        kk=kk+1; if (ii+jj==2); list_names=cat(2,list_names,{'SAR'}); end
    f(ii,jj, kk) = fdarcyvak(ks,D,Q,v);        kk=kk+1; if (ii+jj==2); list_names=cat(2,list_names,{'VAK'}); end
    f(ii,jj, kk) = fdarcywood(ks,D,Q,v);       kk=kk+1; if (ii+jj==2); list_names=cat(2,list_names,{'WOO'}); end
    f(ii,jj, kk) = fdarcyzigrang1(ks,D,Q,v);   kk=kk+1; if (ii+jj==2); list_names=cat(2,list_names,{'ZIG1'}); end
    f(ii,jj, kk) = fdarcyzigrang2(ks,D,Q,v);   kk=kk+1; if (ii+jj==2); list_names=cat(2,list_names,{'ZIG2'}); end
  end
end


%% Figure of the errors / absolute errors

for ii=1:neqs-1;
  % Absolute error [%]
  fobj(:,:,ii)=100*abs(f(:,:,1)-f(:,:,ii+1))./f(:,:,1);
  % Error [%]
  %fobj(:,:,ii)=100*(f(:,:,1)-f(:,:,ii+1))./f(:,:,1);
  if ii<10;
    figure(1);
    set(gcf,'Position',[63 1 1538 833]);
  else
    if ii<19;
      figure(2);
      set(gcf,'Position',[63 1 1538 833]);
    else
      if ii<28;
        figure(3);
        set(gcf,'Position',[63 1 1538 833]);
      else
        figure(4);
        set(gcf,'Position',[63 1 1538 833]);
      end
    end
  end
  sp = mod(ii,9); if sp==0; sp=9; end;
  subplot(3,3,sp);
  %subplot(5,4,ii);
  surf(log10(4/pi)-vi,ksi,fobj(:,:,ii));
  %[c, h] = contour(log10(4/pi)-vi,ksi,fobj(:,:,ii)); clabel(c, h); %colorbar;
  t1=title(['[',num2str(ii),']   ',list_names{ii+1}]);
  set(t1,'FontName','Arial');
  set(t1,'FontSize',13);
  set(t1,'FontWeight','Bold');
  axis([3.602 8 -6 -1]); 
  set(gca,'FontName','Arial');
  set(gca,'FontSize',12);
  set(gca,'FontWeight','Bold');
  xlabel('Log_{10}( \it{Re} )');
  set(gca,'xtick',(4:1:8));
  set(gca,'xticklabel',(4:1:8));
  ylabel(['Log_{10}(\it{ k_{s} / D} )']);
  set(gca,'ytick',(-6:1:-1));
  set(gca,'yticklabel',(-6:1:-1));
  view(0,90); shading('flat'); 
  cm=colormap('jet'); 
  cm(1,:)=[0,0,0];
  cc=colorbar;
  set(cc,'FontName','Arial');
  set(cc,'FontSize',12);
  set(cc,'Fontweight','Bold');
  hold on;
end
hold off;
