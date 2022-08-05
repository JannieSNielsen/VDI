%calcULS - Edited version
% calculation of extreme load in DLC 1.1 using load extrapolation based on full distribution of wind speed
% and turbulence

clc
clear all
close all

Itype=2; %0: full turbulence distribution I,  2: characteristic value,
muV=10; %mean wind speed [m/s] class: 1:10, 2:8.5, 3:7.5
alphaV=2*muV/sqrt(pi);
npmax=0; %number of extremes to use for each windspeed - 0: all
%alphaV=11.28
%alphaV=9;
kV=2;%2.3;%2;
Iref=0.14; %turbulence intensity A: 0.16, B: 0.14, C: 0.12

%1.3761e+05 
%9.5221e+04 etm

%% Create file localPeaks_x using "extractPeaks.m"
% if Itype==0
    load id_I
    load localPeaks_I2
% elseif Itype==2
%     load id_Ic
%     load localPeaks_Ic2
% end

% c1=2;
% sigma_etm=c1*Iref*(0.072*(muV./c1+3).*(V./c1-4)+10);
% Iref_etm=sigma_etm./(0.75.*V+5.6)
% Iref=0.22;

%% Fit distribution function to obtain local distribution for peaks for each wind speed

plot_2=0;

[ALFA,BETA,GAMMA,NN,MAXX,max_ts]=fitLocal(localPeaks,V,Irefp,plot_2,npmax);

%i=5; j=5;
%[maxV,k_max]=max(squeeze(max_ts(i,j,:))) %find simulation with highest
%response

      %%      
figure
plot(V,MAXX)

figure
plot(Irefp,MAXX)

%% Create weights based on joint distribution of wind speed and turbulence

if Itype==0
    II=Irefp>0; %logical with ones
elseif Itype==2
   II=abs(Irefp-Iref)<0.00001; 
end

%Wind speed: Rayleight distribution
FV=1-exp(-(Vint./alphaV).^kV);
PV=diff(FV);

if Itype<=1
    %Turbulence: Weibull
    k=0.27*V+1.4;
    c=3.3; %[m/s]
    C=Iref.*(0.75.*V+c);
    Fsigma=1-exp(-(sigma_int./C).^k);
    Psigma=diff(Fsigma,1,2);
    weight=PV.*Psigma;
elseif Itype>1
    weight=PV;
end
TotPweight=sum(weight,'all');

%%
% CALCULATION OF PROBABILITY CORRESPONDING TO CHARACTERISTIC RESPONSE
Fchar = 1 - 10 / (60*24*365*50);
Prob=1-Fchar;
Lchar_hst = 1.6501e+05; %from HST program

%% Calculation of long term distribution
dL=100;
L=(0:dL:200000)';
Fshort=calcFshort(L,ALFA(:,II),BETA(:,II),GAMMA(:,II),NN(:,II),nsim); %short term distribution for all 10 min conditions

Fshorts=Fshort.*weight./TotPweight;
Flong=squeeze(sum(Fshorts,[1 2]));

Fshorts_log=log(Fshort).*weight;%./TotPweight;
Fshorts_log(isnan(Fshorts_log))=0;
Flong_log=exp(squeeze(sum(Fshorts_log,[1 2])));

Pshorts=diff(Fshorts,1,3);
Plong=diff(Flong);
Plong_log=diff(Flong_log);

Lmean=L(2:end)'*Plong;
Lstd=sqrt(Plong'*(L(2:end)-Lmean).^2);
Lcov=Lstd/Lmean

%Lchar=interp1(Flong(Flong>0 & Flong<1),L(Flong>0 & Flong<1),Fchar)
%Lchar=Lchar_hst;


% figure
% semilogy(L,1-Fshort)
% hold on
% semilogy([min(L) max(L)],Prob*[1 1],Lchar,Prob,'*k')
% grid on
% legend(num2str(windspeed'))
% xlabel('Maximum 10-min load (conditioned on wind speed)')
% ylabel('1-cdf')


figure(10)
semilogy(L,1-Flong)
hold on
semilogy(L,1-Flong_log,'--')
%semilogy([min(L) max(L)],Prob*[1 1],Lchar,Prob,'*k')
grid on
xlabel('Maximum 10-min load')
ylabel('1-cdf')

figure
if Itype==2
plot(L(2:end),squeeze(Pshorts)')
legend(num2str(V))
end
hold on

plot(L(2:end),Plong,'k','Displayname','Aggregate')
plot(L(2:end),Plong_log,'k--','Displayname','Log-based aggregate')
xlabel('Maximum 10-min load')
ylabel('pdf')


 figure
if Itype==2
semilogy(L(2:end),squeeze(Pshorts)')
legend(num2str(V))
end
hold on

semilogy(L(2:end),Plong,'k','Displayname','Aggregate')
semilogy(L(2:end),Plong_log,'k--','Displayname','Log-based aggregate')
xlabel('Maximum 10-min load')
ylabel('pdf')



if Itype==0
 figure
 for j=1:nI
 subplot(3,3,j)
 semilogy(L(2:end),squeeze(Pshorts(:,j,:))')
  hold on
 semilogy(L(2:end),Plong,'k')
  semilogy(L(2:end),Plong_log,'k--')
xlabel('Maximum 10-min load')
ylabel('pdf')
title([num2str(j) ': ' num2str(Irefp(j))])
set(gca,'Ylim',[1e-20 1])
 end
legend(num2str((1:nV)'./100+V))
end

%% Calculation of distribution for max annual load

Flong_yr=Flong.^(6*24*365);
Plong_yr=diff(Flong_yr);

Flong_log_yr=Flong_log.^(6*24*365);
Plong_log_yr=diff(Flong_log_yr);

figure(11)
plot(L,Flong_yr)
hold on
plot(L,Flong_log_yr,'--')
plot([min(L) max(L)],0.98*[1 1])%,Lchar,0.98,'*k')
grid on
xlabel('Maximum annual load')
ylabel('cdf')

figure(12)
hold on
plot(L(2:end),Plong_yr)
hold on
plot(L(2:end),Plong_log_yr,'--')
xlabel('Maximum annual load')
ylabel('pdf')

%% Approximate by Gumbel distribution

%use distributions found with log aggregation
Flong_yr=Flong_log_yr;%(isnan(Flong_log_yr)==0);
Plong_yr=Plong_log_yr;%(isnan(Plong_log_yr)==0);

Lmean_yr=L(2:end)'*Plong_yr
Lstd_yr=sqrt(Plong_yr'*(L(2:end)-Lmean_yr).^2);
Lcov_yr=Lstd_yr/Lmean_yr
%Lchar_yr=interp1(Flong_yr(Flong_yr>0),L(Flong_yr>0),0.98)

%L_sigma=Lmean_yr*L_cov;
L_a=pi/(Lstd_yr*sqrt(6));
L_b=Lmean_yr-0.5772/L_a;
Lc=L_b-1/L_a*log(-log(0.98))

Flong_yr_gumb=exp(-exp(-L_a*(L-L_b)));
Plong_yr_gumb=diff(Flong_yr_gumb);

figure(20)
plot(L,Flong_yr,'k')
hold on
%plot(L,Flong_log_yr,'b--')
plot(L,Flong_yr_gumb,'-r')
plot([min(L) max(L)],0.98*[1 1],Lc,0.98,'*r')%Lchar,1-0.98,'*k'
grid on
xlabel('Maximum annual load')
ylabel('cdf')

figure(21)
hold on
plot(L(2:end),Plong_yr,'k')
hold on
%plot(L(2:end),Plong_log_yr,'b--')
plot(L(2:end),Plong_yr_gumb,'-r')
xlabel('Maximum annual load')
ylabel('pdf')

figure(22)
semilogy(L,1-Flong_yr,'k')
hold on
%semilogy(L,1-Flong_log_yr,'b--')
semilogy(L,1-Flong_yr_gumb,'-r')
semilogy([min(L) max(L)],1-0.98*[1 1],Lc,1-0.98,'*r')%Lchar,1-0.98,'*k'
grid on
xlabel('Maximum annual load')
ylabel('1-cdf')

figure(100)
hold on
plot(L(2:end),Plong_yr_gumb,'-')
xlabel('Maximum annual load')
ylabel('pdf')

figure(101)
hold on
semilogy(L,1-Flong_yr_gumb)
grid on
xlabel('Maximum annual load')
ylabel('1-cdf')

%%
%figure(101)
%semilogy([min(L) max(L)],1-0.98*[1 1],'k')%,Lc,1-0.98,'*k')%Lchar,1-0.98,'*k'

% %% Aggregate first, when fit distribution for 10 min max
% 
% 
% 
% [nV,nI,nsim]=size(localPeaks);
% clear peaksL peaksV peaksI peaksW
% 
% np=1;
% for i=1:nV
%     for j=1:nI
%         for k=1:nsim
%             n=length(localPeaks{i,j,k});
%             peaksL(np:np+n-1)=localPeaks{i,j,k};
%             peaksV(np:np+n-1)=i;
%             peaksI(np:np+n-1)=j;
%             peaksW(np:np+n-1)=weight(i,j);
%             np=np+n;
%         end
%     end
% end
% 
% [peaksL,Is]=sort(peaksL);
% peaksV=peaksV(Is);
% peaksI=peaksI(Is);
% peaksW=peaksW(Is);
% % 
% %nsamp=10000;
% nsamp=length(peaksL);
% 
% peaksL=peaksL(end-nsamp+1:end)';
% peaksV=peaksV(end-nsamp+1:end)';
% peaksI=peaksI(end-nsamp+1:end)';
% peaksW=peaksW(end-nsamp+1:end)';
% peaksW=peaksW/max(peaksW);
% [peaksL peaksV peaksI];
% %peaksL=peaksL*1e-5;
% 
% %F(x)=exp(-exp(-a(x-b)))
% %x=1:0.001:2;
% x=0:1000:3e5;
% %a=5e4; b=5e4;
% %f=a.*exp(-((a.*(x-b))+exp(-(a.*(x-b)))))
% %f=gumbelpdf(x,mu,beta)
% 
% parm=fmincon(@(par) -sum(peaksW.*log(gumbelpdf(peaksL,par(1),par(2)))),[10000 20000],[],[],[],[],[],[])
% %parm=fmincon(@(par) -sum(peaksW.*log(gumbelpdf(peaksL,par(1),par(2)))),[0.1 0.2],[],[],[],[],[],[])
% 
% %parm=fmincon(@(par) sum(log(gumbelpdf(peaksL,par(1),par(2)))),[mean(peaksL) mean(peaksL)],[],[],[],[],[],[])
% 
% Fshort2=gumbelcdf(x,parm(1),parm(2)).^(length(peaksW)/(nsim*nV*nI));
% Pshort2=diff(Fshort2);
% 
% figure(100)
% plot(x,gumbelpdf(x,parm(1),parm(2))/sum(gumbelpdf(x,parm(1),parm(2))))
% hold on
% plot(peaksL,zeros(size(peaksL)),'*')
% plot(x(2:end),Pshort2)
% 
% figure(101)
% plot(peaksL,peaksW,'*')
% 
% figure(102)
% semilogy(x,1-gumbelcdf(x,parm(1),parm(2)))
% hold on
% grid on
% 
% figure(103)
% plot(x,Fshort2)
% hold on
% grid on
% 
% figure(104)
% plot(x(2:end),Pshort2)
% hold on
% grid on

%%
% 
% Lcov_yr =
% 
%     0.0626
% 
% Lchar_yr =
% 
%    1.6131e+05
% 
% 
% 
% Lmean_yr =
% 
%    1.3036e+05
% 
% 
% Lcov_yr =
% 
%     0.0648
% 
% 
% Lc =
% 
%    1.5227e+05