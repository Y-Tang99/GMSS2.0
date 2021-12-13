%=========================================================================%
% Ground Motion Simulation System (GMSS) Version 2.0 -- Stochastic Finite
% Fault Source Simulation
%=========================================================================%

% THIS IS FOR SINGLE-SCENARIO SIMULATION

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       ^^ INPUT PARAMETERS ^^                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-------------------------------------------------------------------------%
% (1) General Input Parameters                                            %
%-------------------------------------------------------------------------%

Nhyp=2;           % number of hypocenters
NS=3;             % number of simulations

M=6.1;             % Moment magnitude
stress=28;         % Stress drop, bars
Vs30=0.76;         % time-averaged shear wave velocity over top 30 m, km/s
kappa=0.06;        % high-frequency decay slope, s

beta0=3.5;         % Source shear wave velocity (SWV), km/s
stress_ref=70;     % reference stress drop, used for finding FL & FW
roll0=2.8;         % Source density, g/cm^3 
y=0.8;             % ratio of rupture velocity and Source SWV
vrup=y*beta0;      % rupture velocity
z=1.68;            % x=0.5
             
dt=0.005;          % Time step, no larger than 0.02 s

%%%---------------------------------------------------------------------%%%
% Note: To be more strict, the time pad should be computed using the following 
% function: tpad=1.5*n*f0
% where n is the order of Butterworth filter, usually using 4 or 8; and f0 
% is the corner-frequency (Boore, BSSA, 2005).
%%%---------------------------------------------------------------------%%%

tpadl=20;          % Time pad before simulated series
tpadt=20;          % Time pad after simulated series
npadl=tpadl/dt;
npadt=tpadt/dt;

%-------------------------------------------------------------------------%
% (2) Finite-fault Input Parameters                                       %
%-------------------------------------------------------------------------%

FaultLat=27.11;    % latitude of upper edge of fault
FaultLon=103.35;   % longitude of upper edge of fault
Fstrike=162;       % fault strike,degree (°)
Fdip=86;           % fault dip, degree (°)
Rake=45;           % rake angle, degree (°)

% Fault Dimensions 

%FDfactor=1;       % use customdefined fault dimensions
% FL0=55;           % fault length
% FW0=40;           % fault width

% FDfactor=2;      % use relation dveloped by Wells & Coppersmith (1994)
% FDfactor=3;        % use relation dveloped by Leonard (2010)
FDfactor=4;         % use relation dveloped by Kumar et al. (2017)
% FDfactor=5;      % use relation dveloped by Cheng et al. (2019)

% The following parameters are needed to locate the origin point

h_ref=10;          % fault depth to reference point
h_min=3.0;         % Campbell depth to seismogenic region, usually set as 3.0

% Subfault Dimension
dl=2;              % subfault length, no less than 1.5 km
dw=2;              % subfault width, no less than 1.5 km

pp=0.5;            % pulsing percentage

%-------------------------------------------------------------------------%
% (3) Site Input Parameters                                               %
%-------------------------------------------------------------------------%
% Site location

% Two options for determining site location: lattitude & longitude (LatLon),
% and distance & azimuth (DistAz). Users need to choose one for their purposes.

 SLfactor=1;    % Latitude and longitude
% SLfactor=2;   % Distance and Azimuth

SL1=27.299;           % site/station latitude
SL2=103.699;          % site/station longitude

%    SL1=5;           % Input values to get site location
%    SL2=0;

%-------------------------------------------------------------------------%
% (4) Input Parameters for Fourier domain and Response Spectra Output     %
%-------------------------------------------------------------------------%

fomin=0.1;        % lower output frequency
fomax=100;        % upper output frequency

Ntn=100;          % number of natural periods
psi=0.05;         % Viscous damping value

%-------------------------------------------------------------------------%
% (5) Input Parameters for selecting the simulation methods               %
%-------------------------------------------------------------------------%

f0factor=1; % use the rupture velocity dependent corner-frequency and the 
            % corresponding rise time (Tang, BSSA, 2021)
% f0factor=2; % Original EXSIM

% Suggest to use 1 or 2 for BLfactor
BLfactor=1; % use the baseline correction method proposed by Boore (BSSA, 2005)
% BLfactor=2; % use the baseline correction method proposed by Chiu (BSSA, 1997)
% BLfactor=3; % use the baseline correction method proposed by Graizer,V.M.(1979)
% BLfactor=4; % use the baseline correction method proposed by Jiang (IEM, 2010)
% BLfactor=5; % use the baseline correction method proposed by Papazafeiropoulos & Plevris (OpenSeismoMatlab,2018)

%--------------------------------------------------------------------------
%-------This is the end of INPUT PARAMETERS section-----------------------
%--------------------------------------------------------------------------

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       ^^ BASIC CALCULATIONS ^^
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if FDfactor==1       % usercustom fault dimension   
    FL=FL0;           % fault length
    FW=FW0;           % fault width
elseif FDfactor==2   % Wells & Coppersmith (1994)
    if Rake == 0 || Rake == 180                                  %% Strike Slip
        FL=10^(-2.57+0.62*M)*(stress_ref/stress)^(1/3);
        FW=10^(-0.76+0.27*M)*(stress_ref/stress)^(1/3);
    elseif (Rake > 0) && (Rake < 180) && Fdip ~=0 && Fdip ~= 90  %% Reverse
        FL=10^(-2.42+0.58*M)*(stress_ref/stress)^(1/3);
        FW=10^(-1.61+0.41*M)*(stress_ref/stress)^(1/3);
    elseif (Rake > -180) && (Rake < 0) && Fdip ~=0 && Fdip ~= 90 %% Normal
        FL=10^(-1.88+0.50*M)*(stress_ref/stress)^(1/3);
        FW=10^(-1.14+0.35*M)*(stress_ref/stress)^(1/3);
    else                                                         %% Undefined
        FL=10^(-2.44+0.59*M)*(stress_ref/stress)^(1/3);
        FW=10^(-1.01+0.32*M)*(stress_ref/stress)^(1/3);
    end
elseif FDfactor==3  % Leonard's correlation (2010)
    if Rake == 0 || Rake == 180                                  %% Strike Slip
        FL=10^(-2.5+0.6*M)*(stress_ref/stress)^(1/3);
        FW=10^(-1.49+0.4*M)*(stress_ref/stress)^(1/3);
    elseif (Rake > 0) && (Rake < 180) && Fdip ~=0 && Fdip ~= 90  %% Reverse
        FL=10^(-2.54+0.6*M)*(stress_ref/stress)^(1/3);
        FW=10^(-1.46+0.4*M)*(stress_ref/stress)^(1/3);
    elseif (Rake > -180) && (Rake < 0) && Fdip ~=0 && Fdip ~= 90 %% Normal
        FL=10^(-2.54+0.60*M)*(stress_ref/stress)^(1/3);
        FW=10^(-1.46+0.4*M)*(stress_ref/stress)^(1/3);
    else                                                         %% SCR
        FL=10^(-2.59+0.6*M)*(stress_ref/stress)^(1/3);
        FW=10^(-1.6+0.4*M)*(stress_ref/stress)^(1/3);
    end
elseif FDfactor==4  % Kumar et al.'s correlation (2017)
    if Rake == 0 || Rake == 180                                  %% Strike Slip
        FL=10^(-2.943+0.681*M)*(stress_ref/stress)^(1/3);
        FW=10^(-0.543+0.261*M)*(stress_ref/stress)^(1/3);
    elseif (Rake > 0) && (Rake < 180) && Fdip ~=0 && Fdip ~= 90  %% Reverse
        FL=10^(-2.693+0.614*M)*(stress_ref/stress)^(1/3);
        FW=10^(-1.669+0.435*M)*(stress_ref/stress)^(1/3);
    elseif (Rake > -180) && (Rake < 0) && Fdip ~=0 && Fdip ~= 90 %% Normal
        FL=10^(-1.722+0.485*M)*(stress_ref/stress)^(1/3);
        FW=10^(-0.829+0.323*M)*(stress_ref/stress)^(1/3);
    else                                                         %% Subduction interface
        FL=10^(-2.412+0.583*M)*(stress_ref/stress)^(1/3);
        FW=10^(-0.88+0.366*M)*(stress_ref/stress)^(1/3);
    end 
else                % Cheng et al.'s correlation (2019) (for mainland China)
    if Rake == 0 || Rake == 180                                  %% Strike Slip
        FL=10^(-2.45+0.61*M)*(stress_ref/stress)^(1/3);
        FW=10^(-1.38+0.41*M)*(stress_ref/stress)^(1/3);
    elseif (Rake > 0) && (Rake < 180) && Fdip ~=0 && Fdip ~= 90  %% Reverse
        FL=10^(-3.27+0.72*M)*(stress_ref/stress)^(1/3);
        FW=10^(-1.67+0.44*M)*(stress_ref/stress)^(1/3);
    elseif (Rake > -180) && (Rake < 0) && Fdip ~=0 && Fdip ~= 90 %% Normal
        FL=10^(-4.02+0.83*M)*(stress_ref/stress)^(1/3);
        FW=10^(-2.13+0.51*M)*(stress_ref/stress)^(1/3);
    else                                                         %% Undefined
        FL=10^(-2.67+0.63*M)*(stress_ref/stress)^(1/3);
        FW=10^(-1.38+0.4*M)*(stress_ref/stress)^(1/3);
    end
end

if Fstrike==0
    s1f=-FL/2;     % Along strike near edge
    s2f=FL/2;      % Along strike far edge 
    w1f=-FW/2;     % Down dip near edge
    w2f=FW/2;      % Down dip far edge
else
    s1f=0;         % Along strike near edge
    s2f=FL;        % Along strike far edge 
    w1f=0;         % Down dip near edge
    w2f=FW;        % Down dip far edge
end

nl=round(FL/dl);   % number of subfaults along strike
nw=round(FW/dw);   % number of subfaults along dip
if nl<=1
    nl=1;
    dl=FL;
end
if nw<=1
    nw=1;
    dw=FW;
end

NF=nl*nw;           % total number of subfaults

fo=logspace(log10(fomin),log10(fomax),Ntn);   % Output frequency
Tn=1./fo;
Tn=Tn(end:-1:1);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       ^^ MAIN PROGRAM ^^
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-------------------------------------------------------------------------%
% (1) Dtermination of Various Distances                                   %
%-------------------------------------------------------------------------%

% R: Distance with respect to the origin point
% Az: azimuth
% Rrup: Closest distance to fault rupture
% Rjb: Closet distance to the surface projection of the fault rupture
% Rseis: Closest distance to rupture portion in seismogenic crust
% Rx: Horizontal distance perpendicular to the rupture strike (site
% coordinate, positive on hanging wall side)
% Ry0: A distance along strike, used by Abrahamson et al. (2014)to taper
% their hanging wall effect.  It can only be greater than or equal to 0.0. 

[SiteLat,SiteLon,R,Az]=FUNSL(SLfactor,SL1,SL2,FaultLat,FaultLon);

[Rrup,Rjb,Rseis,Rx,Ry0,Azrup,Azjb,Azseis]=FUNDist(SiteLat,SiteLon,FaultLat,FaultLon,Fstrike,Fdip,h_ref,s1f,s2f,w1f,w2f,h_min);

%-------------------------------------------------------------------------%
% (2) Dtermination of Number of Time Steps                                %
%-------------------------------------------------------------------------%

C=10^(-20)*(0.55*2.0*0.707)/(4*pi*roll0*beta0^3);  % be careful of the unit of C
M0=10^(1.5*M+16.05);                               % total seismic moment
aveM0=M0/NF;

if f0factor==1
    f0=4.2*(10^6)*y*z*beta0*((stress/M0)^(1/3));         % total corner frequency
    firstf0=4.2*(10^6)*y*z*beta0*((stress/aveM0)^(1/3));
else
    f0=4.906*(10^6)*beta0*((stress/M0)^(1/3));           % total corner frequency
    firstf0=4.906*(10^6)*beta0*((stress/aveM0)^(1/3));
end


% (2.1) Locate hypocenter randomly

i0(Nhyp)=zeros();
j0(Nhyp)=zeros();
Rhyp(Nhyp)=zeros();

for nn=1:1:Nhyp
%     s=rng;
    randomN=rand(1);
    i0(nn)=round(randomN*nl);
    if i0(nn)<1
        i0(nn)=1;
    elseif i0(nn)>nl
        i0(nn)=nl;
    end
    j0(nn)=round(randomN*nw);
    if j0(nn)<1
        j0(nn)=1;
    elseif j0(nn)>nw
        j0(nn)=nw;
    end
    % hypocentral distance
    Rhyp(nn)=FUNsubR(R,h_ref,Fdip,Fstrike,Az,dl,dw,i0(nn),j0(nn)); 
end

% (2.2) Find the subfault moment

slipweight(Nhyp,nl,nw)=zeros();
totalweight(Nhyp)=zeros();
subM0(Nhyp,nl,nw)=zeros();
for nn=1:1:Nhyp
    for i=1:1:nl
        for j=1:1:nw
        slipweight(nn,i,j)=rand(1);
        end
    end
     totalweight(nn)=sum(slipweight(nn,:));
end

for nn=1:1:Nhyp
    for i=1:1:nl
        for j=1:1:nw
            subM0(nn,i,j)=M0*slipweight(nn,i,j)/totalweight(nn);
        end
    end
end


% (2.3) Find the number of active subfaults

NoEsubfault=round(nl*pp/2);  % number of effective subfaults
if NoEsubfault<=1
    NoEsubfault=1;
end

% number of active subfaults, a function of t
NR(Nhyp,nl,nw)=zeros();   
for nn=1:1:Nhyp
    for i=1:1:nl
        for j=1:1:nw
            NR(nn,i,j)=FUNNP(i,j,i0(nn),j0(nn),nl,nw,NoEsubfault);
            if NR(nn,i,j)==0
                NR(nn,i,j)=1;
            end
        end
    end
end

% (2.4) The following is for calculaing: subfault corner frequency, subfault 
% rise time subfault actual slip, subfault distance, subfault duration, 
% subfault delay time, arrive time, and end time.


subf0(Nhyp,nl,nw)=zeros();               % subfault corner frequency
subTrise(Nhyp,nl,nw)=zeros();            % subfault rise time
subTpath(Nhyp,nl,nw)=zeros();            % subfault path duration
subRadius=sqrt((dl*dw)/pi);              % subfault radius
subSlip(Nhyp,nl,nw)=zeros();             % actual subfault slip
subR(Nhyp,nl,nw)=zeros();                % subfault distance
subdur(Nhyp,nl,nw)=zeros();              % subfault duration
subdelay(Nhyp,nl,nw)=zeros();            % subfault delay time
subTarrive(Nhyp,nl,nw)=zeros();          % subfault arrive time
subTend(Nhyp,nl,nw)=zeros();             % subfault end time


for nn=1:1:Nhyp
    for i=1:1:nl
        for j=1:1:nw
            subf0(nn,i,j)=firstf0*NR(nn,i,j)^(-1/3); 
            subSlip(nn,i,j)=10^(-22)*subM0(nn,i,j)/(roll0*(beta0^2)*dl*dw);
            subR(nn,i,j)=FUNsubR(R,h_ref,Fdip,Fstrike,Az,dl,dw,i,j);
            subTpath(nn,i,j)=FUNTpath(subR(nn,i,j));
            if f0factor==1
                subTrise(nn,i,j)=subRadius/(2*vrup);    % rise time
            else
                subTrise(nn,i,j)=subRadius/(vrup);     % EXSIM original rise time
            end
            subdur(nn,i,j)=subTrise(nn,i,j)+subTpath(nn,i,j);  
            if i-i0(nn)~=0 || j-j0(nn)~=0
                subdelay(nn,i,j)=sqrt((dl*(i-i0(nn)))^2+(dw*(j-j0(nn)))^2)/vrup;
            else
                subdelay(nn,i,j)=0;
            end
            subTarrive(nn,i,j)=subdelay(nn,i,j)+subR(nn,i,j)/beta0;
            subTend(nn,i,j)=subTarrive(nn,i,j)+subdur(nn,i,j);
        end
    end
end

Trise_max=max(max(max(subTrise)));
Tarrive_min=min(min(min(subTarrive)));
Tend_max=max(max(max(subTend)));

NtotalWave=(Tend_max-Tarrive_min+tpadl+tpadt+Trise_max)/dt; % total number of dt

nseed=round(log2(NtotalWave));
if 2^nseed >= NtotalWave
    NT=2^nseed;
else
    NT=2^(nseed+1);
end

%-------------------------------------------------------------------------%
% (3) Main Procedures                                                     %
%-------------------------------------------------------------------------%

T=NT*dt;
t=dt:dt:T;
df=1/T;
f=(1:1:NT/2)*df;

At(Nhyp,NS,NT)=zeros();     % Time series of acceleration,       cm/s^2
Vt(Nhyp,NS,NT)=zeros();     % Time series of velocity,           cm/s
Dt(Nhyp,NS,NT)=zeros();     % Time series of displacement,       cm
AI0(Nhyp,NS,NT)=zeros();
AI(Nhyp,NS)=zeros();        % Arias Intensity,                   cm/s
Husid(Nhyp,NS,NT)=zeros();
PGA(Nhyp,NS)=zeros();       % Peak ground acceleration,          g     
PGV(Nhyp,NS)=zeros();       % Peak ground velocity,              cm/s
PGD(Nhyp,NS)=zeros();       % Peak ground displacement,          cm
FAS0(Nhyp,NS,NT)=zeros();
FAS(Nhyp,NS,NT/2)=zeros();  % Fourier amplitude spectrum,cm/s
Sa(Nhyp,NS,Ntn)=zeros();    % Response spectrum of acceleration, g
Sv(Nhyp,NS,Ntn)=zeros();    % Response spectrum of velocity,     cm/s
Sd(Nhyp,NS,Ntn)=zeros();    % Response spectrum of displacement, cm

for nn=1:1:Nhyp
    for m=1:1:NS
        % Time Series
        At(nn,m,:) = FUNMainPro(i0(nn),j0(nn),M0,nl,nw,pp,firstf0,vrup,roll0,beta0,C,dl,dw,R,h_ref,Fdip,Fstrike,Az,npadl,npadt,dt,f0,NF,Vs30,kappa,Tarrive_min,NT,f0factor,BLfactor);
        Vt(nn,m,:) = cumsum(At(nn,m,:))*dt; 
        Dt(nn,m,:) = cumsum(Vt(nn,m,:))*dt; 
        % Peak Ground Motions
        PGA(nn,m) = max(abs(At(nn,m,:)))/980; % unit "g"
        PGV(nn,m) = max(abs(Vt(nn,m,:)));     % unit "cm/s"
        PGD(nn,m) = max(abs(Dt(nn,m,:)));     % unit "cm"
        % Arias Intensity
        AI0(nn,m,:) = cumsum(At(nn,m,:).^2)*pi*dt/(2*980);  % unit "cm/s"
        AI(nn,m) = AI0(nn,m,end);
        Husid(nn,m,:) = AI0(nn,m,:)/AI(nn,m); 
        % Response Spectra
        [Sd(nn,m,:),Sv(nn,m,:),Sa(nn,m,:)]=CDM(At(nn,m,:),psi,dt,Tn);
        % Fourier Amplitude Spectra
        FAS0(nn,m,:) = abs(fft(At(nn,m,:),NT))*dt; % Fourier amplitudes 
        FAS(nn,m,:) = FAS0(nn,m,1:NT/2);                 % single sided FAS
    end
end

%-------------------------------------------------------------------------%
% (4) Calculation of Averages                                             %
%-------------------------------------------------------------------------%

% Arithmetic average over number of simulations (NS), 
% and geometric average over number of hypocenters (Nhyp)

if Nhyp==1 && NS==1
    meanPGA=PGA;
    meanPGV=PGV;
    meanAI0=AI0;
    meanAI=AI;
    meanHusid=Husid;
    meanFAS=FAS;
    meanSa=Sa/980;
    meanSv=Sv;
    meanSd=Sd;
else
  if Nhyp==1 && NS~=1
      meanPGA=mean(PGA);
      meanPGV=mean(PGV);
      meanAI0=mean(AI0);
      meanAI=mean(AI);
      meanHusid=rms(Husid);
      meanFAS=rms(FAS);
      meanSa=mean(Sa)/980;
      meanSv=mean(Sv);
      meanSd=mean(Sd);
  else
      if Nhyp~=1 && NS==1
          meanPGA=geomean(PGA);
          meanPGV=geomean(PGV);
          meanAI0=geomean(AI0);
          meanAI=geomean(AI);
          meanHusid=geomean(Husid);
          meanFAS=geomean(FAS);
          meanSa=geomean(Sa)/980;
          meanSv=geomean(Sv);
          meanSd=geomean(Sd);
      else
          meanPGA=geomean(mean(PGA));
          meanPGV=geomean(mean(PGV));
          meanAI0=geomean(mean(AI0));
          meanAI=geomean(mean(AI));
          meanHusid=geomean(rms(Husid));
          meanFAS=geomean(rms(FAS));
          meanSa=geomean(mean(Sa))/980;
          meanSv=geomean(mean(Sv));
          meanSd=geomean(mean(Sd));
      end
  end
end


%-------------------------------------------------------------------------%
% (5) Outputs                                                             %
%-------------------------------------------------------------------------%

% Duration
timed1 = t(meanAI0>=0.05*meanAI & meanAI0<=0.75*meanAI);
Output.t5_75 = [timed1(1),timed1(end)];
Output.D5_75 = timed1(end)-timed1(1);
timed2 = t(meanAI0>=0.05*meanAI & meanAI0<=0.95*meanAI);
Output.t5_95 = [timed2(1),timed2(end)];
Output.D5_95 = timed2(end)-timed2(1);

% Mean Period
fi = f(f>0.25 & f<20);
Ci = meanFAS(f>0.25 & f<20);
Tm = ((Ci(:)'.^2)*(1./fi(:)))/(Ci(:)'*Ci(:));
Output.Tm = Tm;

% Decrease the dimension of final ground motions
meanFAS=reshape(meanFAS,1,NT/2);
meanHusid=reshape(meanHusid,1,NT);
meanSa=reshape(meanSa,1,Ntn);
meanSv=reshape(meanSv,1,Ntn);
meanSd=reshape(meanSd,1,Ntn);
% Time series of first hypocenter and first simulation
At11=reshape(At(1,1,:),1,NT);
Vt11=reshape(Vt(1,1,:),1,NT);
Dt11=reshape(Dt(1,1,:),1,NT);

% Output FAS according to the output frequencies
FASo=interp1(f,meanFAS,fo,'spline');

% Output parameters

Output.Dist_Az.Repi=R;
Output.Dist_Az.Az=Az;
Output.Dist_Az.Rhyp=Rhyp;
Output.Dist_Az.Rrup=Rrup;
Output.Dist_Az.Rjb=Rjb;
Output.Dist_Az.Rseis=Rseis;
Output.Dist_Az.Rx=Rx;
Output.Dist_Az.Ry0=Ry0;

Output.T=t';
Output.F=fo';
Output.Tn=Tn';
Output.At=At11';
Output.Vt=Vt11';
Output.Dt=Dt11';
Output.AI=meanAI';
Output.Husid=meanHusid';
Output.PGA=meanPGA;
Output.PGV=meanPGV; 
Output.FAS=FASo';
Output.Sa=meanSa';
Output.Sv=meanSv';
Output.Sd=meanSd';
Output.Sa11=reshape(Sa(1,1,:)/980,Ntn,1);

%-------------------------------------------------------------------------%
% (6) Outputs Display                                                     %
%-------------------------------------------------------------------------%

figure (1)
subplot(311)
plot(t,At11)
subplot(312)
plot(t,Vt11)
subplot(313)
plot(t,Dt11)

figure (2)
plot(t,meanHusid)

figure(3)
loglog(fo,FASo)

figure(4)
subplot(311)
loglog(Tn,meanSa)
subplot(312)
loglog(Tn,meanSv)
subplot(313)
loglog(Tn,meanSd)

%-------------------------------------------------------------------------%
%----------------This is the end of MAIN PROGRAM section------------------%
%-------------------------------------------------------------------------%

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       ^^ SUBROUTINE FUNCTIONS ^^
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-------------------------------------------------------------------------%
%(1) Functions CANNOT be changed                                          %
%-------------------------------------------------------------------------%

function [SiteLat,SiteLon,R,Az] = FUNSL(SLfactor,SL1,SL2,FaultLat,FaultLon)
% This function is used to determine the site locations

if SLfactor == 1
    SiteLat=SL1;
    SiteLon=SL2;
    % calculate the distance and azimuth of the site with rspect to the origin
    [ArcLen,~]=distance(SiteLat,SiteLon,FaultLat,FaultLon);
    R=ArcLen*6371*pi/180;    % epicentral distance
    R1=6371*(SiteLat-FaultLat)*pi/180;
    if R1>=R
        Az=180;
    else
        Az=acos(R1/R)*180/pi;
    end
    if SiteLon-FaultLon<=0
        Az=360-Az;
    end
else
    if SLfactor == 2
        R=SL1;
        Az=SL2;
        ArcLen=(R/6371)*180/pi;
        [SiteLat,SiteLon]=reckon(FaultLat,FaultLon,ArcLen,Az);
    end
end

end

function [Rrup,Rjb,Rseis,Rx,Ry0,Azrup,Azjb,Azseis] = FUNDist(SiteLat,SiteLon,RefLat,RefLon,Fstrike,Fdip,h_Ref,s1,s2,w1,w2,h_min)
% this function is used for computing various distance measures
% Rrup: Closest distance to fault rupture
% Rjb: Closet distance to the surface projection of the fault rupture
% Rseis: Closest distance to rupture portion in seismogenic crust
% Rx: Horizontal distance perpendicular to the rupture strike (site
% coordinate, positive on hanging wall side)
% Ry0: A distance along strike, used by Abrahamson et al. (2014)to taper
% their hanging wall effect.  It can only be greater than or equal to 0.0. 
% See p. 1040 of Abrahamson et al. (2014), EqSpectra 30, 1025-1055.
% Az: azimuths


d2r=pi()/180;    % convert degree into radians

fstrike=Fstrike*d2r;
fdip=Fdip*d2r;


% compute unit vectors

ix_n=cos(fstrike);
ix_e=sin(fstrike);
ix_d=0;

iy_n=-sin(fstrike)*sin(fdip);
iy_e=cos(fstrike)*sin(fdip);
iy_d=-cos(fdip);

iz_n=-sin(fstrike)*cos(fdip);
iz_e=cos(fstrike)*cos(fdip);
iz_d=sin(fdip);

% convert site lat\lon into northern and estern distances
% simply assume earth is a sphere approximately

dist_n=(SiteLat-RefLat)*d2r*6371;
dist_e=(SiteLon-RefLon)*d2r*cos(0.5*d2r*(SiteLat+RefLat))*6371;
dist_d=-h_Ref;

% Convert coordinates of reference-to-site vector from n,e,d coordinates
% into fault coordinates

rx=dist_n*ix_n+dist_e*ix_e+dist_d*ix_d;
ry=dist_n*iy_n+dist_e*iy_e+dist_d*iy_d;
rz=dist_n*iz_n+dist_e*iz_e+dist_d*iz_d;


% find Rrup and Azrup

[hrup,icaserup]=FUNh(rx,rz,w1,w2,s1,s2);

Rrup=sqrt(hrup^2+ry^2);

% azimuth for Rrup

if icaserup ==1 || icaserup == 2 || icaserup == 3 
    Azrup = atan((rz - w1)/(rx - s2))/d2r;
elseif icaserup == 4
        Azrup = -90.0;
elseif icaserup == 5 || icaserup == 6
        Azrup = 90.0;
elseif icaserup == 7 || icaserup == 8 || icaserup == 9
        Azrup = atan((rz - w1)/(rx - s1))/d2r;
        if Azrup<0
            Azrup=180+Azrup;
        end
end

% find Rseis and Azseis

d2top_c=h_min;
d2top=h_Ref+w1*iz_d;
if d2top<d2top_c && iz_d~=0
    w1_c=(d2top_c - h_Ref)/iz_d;
else
    w1_c=w1;
end
[hseis,icaseseis]=FUNh(rx,rz,w1_c,w2,s1,s2);
Rseis=sqrt(hseis^2+ry^2);

% azimuth for Rseis

if icaseseis ==1 || icaseseis == 2 || icaseseis == 3 
    Azseis = atan((rz - w1_c)/(rx - s2))/d2r;
elseif icaseseis == 4
        Azseis = -90.0;
elseif icaseseis == 5 || icaseseis == 6
        Azseis = 90.0;
elseif icaseseis == 7 || icaseseis == 8 || icaseseis == 9
        Azseis = atan((rz - w1_c)/(rx - s1))/d2r;
        if Azseis<0
            Azseis=180+Azseis;
        end
end

% fine Rjb and Azjb

% find rx, ry, rz, w1, w2, s1, s2 in terms of coordinates of the fault 
% plane projected onto the surface

s1jb=s1;
s2jb=s2;
w1jb=w1*cos(fdip);
w2jb=w2*cos(fdip);
rxjb=rx;
rzjb=-sin(fstrike)*dist_n+cos(fstrike)*dist_e;

[hjb,icasejb]=FUNh(rxjb,rzjb,w1jb,w2jb,s1jb,s2jb);
Rjb=hjb;

% azimuth for Rjb

if icasejb ==1 || icasejb == 2 || icasejb == 3 
    Azjb = atan((rzjb - w1jb)/(rxjb - s2jb))/d2r;
elseif icasejb == 4
        Azjb = -90.0;
elseif icasejb == 5 || icasejb == 6
        Azjb = +90.0;
elseif icasejb == 7 || icasejb == 8 || icasejb == 9
        Azjb = atan((rzjb - w1jb)/(rxjb - s1jb))/d2r;
        if Azjb<0
            Azjb=180+Azjb;
        end
end

% find Rx
     
Rx=rzjb-w1jb;

% find Ry0

if (icasejb == 2 || icasejb == 3 )
        Ry0 = rxjb - s2jb;
elseif (icasejb == 8 || icasejb == 9 ) 
        Ry0 = abs(rxjb - s1jb);
else
        Ry0 = 0;
end

end

function [funh,icase] = FUNh(rx,rz,w1,w2,s1,s2)

 if rx<=s1 && rz<=w1
     icase=7;
     funh = sqrt((s1-rx)^2+(w1-rz)^2);
 elseif rz<=w1 && rx>=s1 && rx<=s2
     icase=4;
     funh = w1-rz;
     elseif rx>=s2 && rz<=w1
         icase=1;
         funh = sqrt((rx-s2)^2+(w1-rz)^2);
         elseif rx>=s2 && rz>=w1 && rz<=w2
             icase=2;
             funh = rx-s2;
             elseif rx>=s2 && rz>=w2
                 icase=3;
                 funh = sqrt((rx-s2)^2+(rz-w2)^2);
                 elseif rz>=w2 && rx>=s1 && rx<=s2
                     icase=6;
                     funh = rz-w2;
                     elseif rz>=w2 && rx<=s1
                         icase=9;
                         funh=sqrt((s1-rx)^2+(rz-w2)^2);
                         elseif rx<=s1 && rz>=w1 && rz<=w2
                             icase=8;
                             funh = s1-rx;
                             elseif rx>=s1 && rx<=s2 && rz>=w1 && rz<=w2
                                 icase=5;
                                 funh = 0.0;
 end
                       
end
                             
function [At] = FUNMainPro(i0,j0,M0,nl,nw,pp,firstf0,vrup,roll0,beta0,C,dl,dw,R,h_Ref,Fdip,Fstrike,Az,npadl,npadt,dt,f0,NF,Vs30,kappa,Tarrive_min,NT,f0factor,BLfactor)

% Find the subfault moment

slipweight(nl,nw)=zeros();
for i=1:1:nl
    for j=1:1:nw
        slipweight(i,j)=rand();
    end
end
totalweight=sum(slipweight(:));
subM0=M0*slipweight./totalweight;

% find the number of active subfaults

NoEsubfault=round(nl*pp/2);  % number of effective subfaults
if NoEsubfault<=1
    NoEsubfault=1;
end

NR(nl,nw)=zeros();   % number of active subfaults, a function of t.
for i=1:1:nl
    for j=1:1:nw
        NR(i,j)=FUNNP(i,j,i0,j0,nl,nw,NoEsubfault);
        if NR(i,j)==0
            NR(i,j)=1;
        end
    end
end


% The following is for calculaing: subfault corner frequency, subfault rise time
% subfault actual slip, subfault distance, subfault duration, subfault
% delay time, arrive time, and end time.

subf0(nl,nw)=zeros();               % subfault corner frequency
subTrise(nl,nw)=zeros();            % subfault rise time
subTpath(nl,nw)=zeros();            % subfault path duration
subRadius=sqrt((dl*dw)/pi);         % subfault radius
subSlip(nl,nw)=zeros();             % actual subfault slip
subR(nl,nw)=zeros();                % subfault distance
subdur(nl,nw)=zeros();              % subfault duration
subdelay(nl,nw)=zeros();            % subfault delay time
subTarrive(nl,nw)=zeros();          % subfault arrive time
subTend(nl,nw)=zeros();             % subfault end time

for i=1:1:nl
    for j=1:1:nw
        subf0(i,j)=firstf0*NR(i,j)^(-1/3); 
        subSlip(i,j)=10^(-22)*subM0(i,j)/(roll0*(beta0^2)*dl*dw);
        subR(i,j)=FUNsubR(R,h_Ref,Fdip,Fstrike,Az,dl,dw,i,j);
        subTpath(i,j)=FUNTpath(subR(i,j));
        if f0factor==1
             subTrise(i,j)=subRadius/(2*vrup);     % rise time
        else
            subTrise(i,j)=subRadius/vrup;     % EXSIM original rise time
        end
        subdur(i,j)=subTrise(i,j)+subTpath(i,j);  
        if i-i0~=0||j-j0~=0
            subdelay(i,j)=sqrt((dl*(i-i0))^2+(dw*(j-j0))^2)/vrup;
        else
            subdelay(i,j)=0;
        end
        subTarrive(i,j)=subdelay(i,j)+subR(i,j)/beta0;
        subTend(i,j)=subTarrive(i,j)+subdur(i,j);
    end
end

% The following is to generate and sum the subfault time series

subNWave(nl,nw)=zeros();
subnseed(nl,nw)=zeros();
subN(nl,nw)=zeros();
stutter(nl,nw)=zeros();
nshift(nl,nw)=zeros();
subAt(nl,nw,NT)=zeros();

for i=1:1:nl
    for j=1:1:nw
        if nl==1 && nw==1
            subN(i,j)=NT;
        else
            subNWave(i,j)=npadl+subdur(i,j)/dt+npadt;
            subnseed(i,j)=round(log2(subNWave(i,j)));
            if 2^subnseed(i,j)>=subNWave(i,j)
                subN(i,j)=2^subnseed(i,j);
            else
                subN(i,j)=2^(subnseed(i,j)+1);
            end
        end
    end
end

for i=1:1:nl
    for j=1:1:nw
        % generation of time series of single subfault, most important!!!
        subAt(i,j,1:subN(i,j))=FUNsubAt(beta0,roll0,C,subN(i,j),dt,subdur(i,j),npadl,subR(i,j),subM0(i,j),f0,subf0(i,j),NF,Vs30,kappa);
        stutter(i,j)=rand(1)*subTrise(i,j);     % random delay of complex slip process
        nshift(i,j)=round((subTarrive(i,j)-Tarrive_min+stutter(i,j))/dt);
        subAt(i,j,:)=circshift(subAt(i,j,:),[0,0,nshift(i,j)]);
    end
end

if nl == 1 && nw == 1
    At00=reshape(subAt,1,subN);
else
    At00=sum(sum(subAt));
end

At00=(At00(:))';

%=========================================================================%
% % Baseline correction
%=========================================================================%
t=(1:length(At00))*dt;

%-------------------------------------------------------------------------%
if BLfactor==1       % Baseline correction (Boore,BSSA 2005)
    df=1/(length(At00)*dt);
    f=(1:length(At00))*df;
    At00_f=fft(At00);
    % Apply a low-cut filter (high-pass)
    fcut=0.05;
    norder=8;
    blcf=1./(1+(fcut./f).^(2*norder));
    At00_f=blcf.*At00_f;
    At00=real(ifft(At00_f));

    [p,~,mu] = polyfit(t,At00,1);   % 1 level polynomial fitting
    baseline = polyval(p,t,[],mu);   
    At=At00-baseline;
%-------------------------------------------------------------------------%    
elseif BLfactor==2   % Baseline correction (Chiu,BSSA 1997)
    fcut=0.05;
    % STEP 1 : least-square fitting in acceleration
    [p,~,mu] = polyfit(t,At00,2);   % 2 level polynomial fitting
    bl0 = polyval(p,t,[],mu);       % baseline, equals to p(1)*t.^2 + p(2)*t+p(3)
    At0=At00-bl0;
    % STEP 2 : high-pass filtering in acceleration
    df=1/(dt*length(At0)); 
    f=(1:length(At0))*df;
    At0_f=fft(At0);
    [b,a] = butter(4,fcut/(f(end)/2),'high'); % 4th order butterworth filter
    At_f = filter(b,a,At0_f);
    At = real(ifft(At_f));

    % % STEP 3 : subtract the initial value in the velocity
    % %          subtract a linear trend from the displacement
    % Dt0 = cumtrapz(cumtrapz(At))*t(2)^2; % disp from acc
    % p_Dt = polyfit(t,Dt0,1);
    % bl_Dt = p_Dt(1)*t + p_Dt(2);
    % Dt = Dt0 - bl_Dt;
%-------------------------------------------------------------------------%    
elseif BLfactor==3   % Baseline correction (Graizer,1979)
    nfit=3;    % nfit level polynomial fitting, nfit is delf-determined
    Vt00=cumsum(At00)*dt;
    fitVt00 = polyfit(t,Vt00,nfit); 
    baselineVt00 = polyval(fitVt00,t);   
    bl=diff(baselineVt00,1);
    baseline=[0,bl];
    At=At00-baseline;
%-------------------------------------------------------------------------%    
elseif BLfactor==4  % Baseline correction used by Institute of Engineering Mechanics,
                    % China Earthquake Administration (Jiang, 2010)
    al=length(t);
    if al<10000
        al=al/2;
    else
        al=al/3;
    end
    Vt00=cumsum(At00)*dt;
    [p,~,mu]=polyfit(t(end-al:end),Vt00(end-al:end),1);
    sr=p(1);  % slope rate
    spt=NT-mean(Vt00(end-400:end))/sr; % the data point begin to tilt
    spt=round(spt);
    [C,I]=max(abs(At00));
    if spt<=I
        spt=I;   % tilt at the peak value
    end
    if spt>NT
        spt=NT;
    end
    if spt>0 && spt<NT
        At=At00-[zeros();sr*ones(NT-spt,1)];
    else
        At=At00;
    end
    Vt=cumsum(At)*dt;
    Dt=cumsum(Vt)*dt;
    [C1,I1]=max(abs(Vt)); % tilt at te beginning
    [p1,~,mu1]=polyfit(t(I1:end),Dt(I1:end),1);
    At(I:I+200)=At(I:I+200)-p1(1)/20;

%-------------------------------------------------------------------------%
else   % Baseline correction provided by OpenSeismoMatlab(George Papazafeiropoulos 2018)
    fitAt00 = polyfit(t,At00,1);   % 1 level polynomial fitting
    baseline0 = polyval(fitAt00,t);   
    At0=At00-baseline0;
    Vt0=cumsum(At0)*dt;
    [p,~,mu] = polyfit(t,Vt0,1);   % 1 level polynomial fitting
    baseline=p(1);
    At=At0-baseline;
end

end

function [NP] = FUNNP(i,j,i0,j0,nl,nw,NoEsubfault)  
% This function is used for determining the number of pulsing subfaults
Rmax=max(abs(i-i0)+1,abs(j-j0)+1);
Rmin=Rmax-NoEsubfault;
if Rmin<0
    Rmin=0;
end
n=0;
for ii=1:1:nl
    for jj=1:1:nw
        r=max(abs(ii-i0)+1,abs(jj-j0)+1);
        if r>Rmin && r<Rmax
            n=n+1;
        end
    end
end
NP=n;
end

function [subR] = FUNsubR(R,h_Ref,Fdip,Fstrike,Az,dl,dw,i,j)
 % This function is used for finding the subfault distance
Fstrike_radians=(Az-Fstrike)*pi/180;
Fdip_radians=(90-Fdip)*pi/180;

t1=R*cos(Fstrike_radians)-(2*i-1)*dl/2;
t2=R*sin(Fstrike_radians)-((2*j-1)*dw/2)*sin(Fdip_radians);
t3=-h_Ref-((2*j-1)*dw/2)*cos(Fdip_radians);
subR=sqrt(t1^2+t2^2+t3^2);
end

function [subAt] = FUNsubAt(beta0,roll0,C,subN,dt,subdur,npadl,subR,subM0,f0,subf0,NF,Vs30,kappa)
% This function is the main parogram for obtaining the time series of
% each subfault

taper=0.05;        
tmax=subN*dt;
subdf=1/tmax;
subNf=subN/2;
subf=(1:1:subNf)*subdf;
ndur=round(subdur/dt);
ntaper=round(taper*ndur);
nstop=ndur+2*ntaper;
%--------------------------------------------------------------------------
% Scaling factor at high frequencies
ScH=sqrt(NF)*(f0/subf0).^2;       % (Boore, 2009)

% Scaling factor at low frequencies
Csc=sqrt(NF)./ScH;
f0eff=subf0./sqrt(Csc);
L1=1+(subf./subf0).^2;
L2=1+(subf./f0eff).^2;
ScL=Csc.*L1./L2;

% total scaling fator 
Sc=ScL.*ScH;   
%--------------------------------------------------------------------------
% Gaussian white noise

nt0=wgn(1,nstop,1);            

% Window function (Sargoni & Hart, 1974)

eps=0.2;eta=0.2;  
nstart=1;

b=-eps*log10(eta)/(1+eps*(log10(eps)-1));
c=b/(eps*subdur);
a=(2.7182818/(eps*subdur))^b;

win(nstop)=zeros();
twin(nstop)=zeros();
twin1(nstop)=zeros();
wf(nstop)=zeros();
st(subN)=zeros();

for k=1:1:nstop
    if k<nstart || k>nstop
        win(k)=1;
    else
        if k>(nstart+ntaper) && k<(nstop-ntaper)
            twin(k)=(subdur/(nstop-nstart-2*ntaper))*(k-nstart-ntaper+1);
            win(k)=a*(twin(k)^b)*exp(-c*twin(k)); 
        else
            if k<=(nstart+ntaper)
                twin1(k)=subdur/(nstop-nstart-2*ntaper);
                wf(k)=a*(twin1(k))^b*exp(-c*twin1(k));
                win(k)=abs(sin((k-nstart)/(ntaper)*pi/2))*wf(k);
            else
                if k>=(nstop-ntaper)
                    wf(k)=a*(subdur)^b*exp(-c*subdur);
                    win(k)=abs(sin((nstop-k)/(ntaper)*pi/2))*wf(k);
                end
            end
        end
    end
    st(k+npadl)=win(k)*nt0(k);
end
% Fourier transform

As=fft(st);
Angle=angle(As(1:subNf));
Asf=abs(As(1:subNf));
Asf=Asf*dt;

% normalize to unit amplitude

Asfsum=sum(Asf.^2);
AveAsf=sqrt(Asfsum/(subNf));
Adj=Asf./AveAsf;                

% Fourier amplitude spectrum

subFAS=InputFAS(beta0,roll0,subf,subR,C,subM0,subf0,Vs30,kappa);

%--------------------------------------------------------------------------
% Apply a low-cut filter (high-pass)
fcut=0.05;
norder=8;
blcf=1./(1+(fcut./subf).^(2*norder));
subFAS=blcf.*subFAS;
% Apply the sclaing factor at high frequencies and the taper facotr at low
% frequencies
subFAS=subFAS.*Sc;
%--------------------------------------------------------------------------
 
Aaf=Adj.*subFAS;

% Inverse Fourier transform

subAt0=ifft(Aaf.*exp(1i*(Angle)),subN,'nonsymmetric')*subN;
subAt=real(subAt0)*2*subdf;     % Be careful of the factor "2"

%==========================================================================
% The Fourier transfrom adopted here is not suggested by D.Boore, which 
% requires to preserve the symmetric component of the Fourier spectrum.
% But I usually remove the symmetric component of the Fourier signal, and
% this does not change the final amplitudes. The only difference is
% that we need to add a factor "2" during the inverse Fourier transform
% procedure.
%--------------------------------------------------------------------------
% The following is an alternative way of Fourier transform (by Gideon)
%--------------------------------------------------------------------------
% % Fourier transform
% 
% As=fft(st);
% Angle=angle(As);
% Asf=dt*abs(As);
% 
% % normalize to unit amplitude
% 
% Asfsum=sum(Asf.^2);
% AveAsf=sqrt(Asfsum/(subN));
% Adj=Asf./AveAsf;      
% 
% % Fourier amplitude spectrum
% 
% subFAS=COMPARE(subf,subR,C,subM0,subf0,Vs30,kappa);
% 
% % Apply the scaling factor
% subFAS=subFAS.*Sc;
% 
% 
% % Make complex
%  Adjcomp=complex(Adj,0);
%  subFAScomp=complex(subFAS,0);
% 
%  Aaf(subN)=zeros();
%  for i=1:1:subN
%     if i<=subNf
%         Aaf(i)=Adjcomp(i)*subFAScomp(i);
%     else
%         Aaf(i)=conj(Aaf(subN+1-i));
%     end
% end

% % Inverse Fourier transform

% subAt0=ifft(Aaf.*exp(1i*(Angle)),subN,'nonsymmetric')*subN;
% subAt=real(subAt0)*subdf;    


% %--------------------------------------------------------------------------
% % Baseline correction (Boore,BSSA 2005)
% %--------------------------------------------------------------------------
% t=(1:length(subAt))*dt;
% [p,~,mu] = polyfit(t,subAt,1);   % 1 level polynomial fitting
% baseline = polyval(p,t,[],mu);   
% subAt=subAt-baseline;

%==========================================================================
end

function [Sd,Sv,Sa] = CDM( at,psi,dt,Tn)
% The following is the implementation of time-step integration method
% Central Difference Method
fr=1./Tn;
wn=fr.*(2*pi);
A=-wn.^2+2/(dt^2);
B=psi*wn./dt-1/(dt^2);
C=1/(dt^2)+psi*wn./dt;
Nr=length(Tn);
NT=length(at);
u(Nr,NT)=zeros();
u(Nr,1)=-at(2)*dt^2/2;
for i=1:1:Nr
    for j=2:1:NT
        u(i,j+1)=(-(at(j))+A(i)*u(i,j)+B(i)*u(i,j-1))/C(i);
    end
end
Sd=max(abs(u'));
Sv=wn.*Sd;
Sa=wn.*Sv;
end

%-------------------------------------------------------------------------%
%(2) Functions NEED to be changed                                         %
%-------------------------------------------------------------------------%

function [subFAS] = InputFAS(beta0,roll0,subf,subR,C,subM0,subf0,~,kappa)    
% This function gives the Fourier amplitude spectrum

% source spectrum
w=(2*subf.*pi).^2;
E1=C*w.*subM0;    
E2=1+(subf./subf0).^2;
E=E1./E2;
 
% geometric spreading

if subR<=50
    G=subR^(-1.0);
elseif subR<=90
    G=((50^0.3)/(50^1.0))*((subR)^(-0.3));
elseif subR<=120
    G=((50^0.3)/(50^1.0))*((90^1.1)/(90^0.3))*(subR^(-1.1));
else
    G=((50^0.3)/(50^1.0))*((90^1.1)/(90^0.3))*((120^0.5)/(120^1.1))*(subR^(-0.5));
end

   
% anelastic attenuation 
% Q0=180;
% nq=0.45;
% Qmin=60;

cq=3.5;        % note this is different from beta0!!
Q0=180;
nq=0.5;
%Qmin=60;
%Q=max(Qmin,Q0*(subf.^(nq)));
Q=Q0*(subf.^nq);

Ae1=-pi*subf.*subR;
Ae2=Q.*cq;
Ae=exp(Ae1./Ae2);

% crustal amplification

% Am=FUNAmfBJ(subf,beta0, roll0, Vs30);

Am=FUNAmfInput(subf,beta0, roll0);

% Note: Users can add any amplification factors, including site
% amplification, empirical amplification, etc.

% High-frequency attenuation

% note kappa=kappa0+dk*(M-Mkappa), in this prpgram dk is set to be 0, thus
% kappa=kappa0. Users can change it if needed. 

An=exp(-pi*subf.*kappa);

% Butterworth low-cut (high-pass) filter, which allows the signals
% at high frequencies and remove the signals at low frequencies

% final subFAS

subFAS=E.*Ae.*Am.*An.*G;     % unit cm/s;

end
 
function [Am] = FUNAmfBJ(f,beta0, roll0, Vs30)

H=[0,0.001,0.002,0.003,0.004,0.005,0.006,0.007,0.008,0.009,0.01,0.012,...
0.014,0.016,0.018,0.02,0.022,0.024,0.026,0.028,0.03,0.032,0.034,0.036,...
0.038,0.04,0.042,0.044,0.046,0.048,0.05,0.052,0.054,0.056,0.058,0.06,...
0.062,0.064,0.066,0.068,0.07,0.072,0.074,0.076,0.078,0.08,0.082,0.104,...
0.126,0.147,0.175,0.202,0.23,0.257,0.289,0.321,0.353,0.385,0.42,0.455,...
0.49,0.525,0.562,0.599,0.637,0.674,0.712,0.751,0.789,0.827,0.866,0.906,...
0.945,0.984,1.02,1.06,1.1,1.14,1.2,1.4,1.6,1.8,2,2.2,2.4,2.6,2.8,3,3.5,...
4,4.5,5,5.5,6,6.5,7,7.5,8,10,15,20,30,50];

% generic very hard rock
V1=[2.768,2.7688,2.7696,2.7704,2.7712,2.772,2.7728,2.7736,2.7744,2.7752...
2.776,2.7776,2.7792,2.7808,2.7824,2.784,2.7856,2.7872,2.7888,2.7904,...
2.792,2.7936,2.7952,2.7968,2.7984,2.8,2.8016,2.8032,2.8048,2.8064,2.808,...
2.80956,2.81112,2.81268,2.81424,2.8158,2.81736,2.81892,2.82048,2.82204,...
2.8236,2.82516,2.82672,2.82828,2.82984,2.8314,2.83296,2.85004,2.86676,...
2.88272,2.9035,2.92344,2.9436,2.9629,2.9853,3.00686,3.06098,3.0821,3.0718,...
3.0941,3.1158,3.1365,3.15796,3.17942,3.20072,3.22048,3.24024,3.260835504,...
3.271637476,3.281964539,3.29211285,3.302087707,3.311425176,3.320409798,...
3.32841313,3.337002316,3.345294257,3.353309507,3.364853485,3.399786096,...
3.430339103,3.457516593,3.482010086,3.50431659,3.510651329,3.516529187,...
3.521980007,3.527062193,3.538443827,3.54833273,3.557078298,3.564919739,...
3.572028075,3.578529853,3.584521359,3.590077571,3.595258021,3.600110775,...
3.616939825,3.647720809,3.669718999,3.700949146,3.740673099];


% generic rock
V2=[0.245,0.245,0.406898105,0.454341586,0.491321564,0.522065927,0.548608647,...
0.572100299,0.593261253,0.612575285,0.630384474,0.662434293,0.690800008,...
0.716351448,0.739672767,0.761177017,0.781168059,0.799876552,0.817482113,...
0.834127602,0.84992869,0.872659052,0.894459078,0.915511227,0.935880677,...
0.955623835,0.974789894,0.993422052,1.011558478,1.02923309,1.046476183,...
1.063314944,1.079773879,1.095875162,1.111638937,1.127083562,1.142225823,...
1.157081117,1.171663602,1.185986333,1.200061379,1.21389992,1.22751234,...
1.240908299,1.254096801,1.26708626,1.279884545,1.409876676,1.524401504,...
1.623105361,1.742468929,1.822101865,1.869784562,1.91154454,1.956709713,...
1.998031044,2.036174042,2.071640608,2.10782397,2.141667263,2.173485515,...
2.203532354,2.23359926,2.262120148,2.289978882,2.315853324,2.341268645,...
2.366247007,2.389604645,2.412077883,2.434298272,2.456270778,2.4769581,...
2.496972474,2.514890987,2.534215818,2.55296508,2.571175941,2.597555276,...
2.678472608,2.750601065,2.815833418,2.875495547,2.930554778,2.981739972,...
3.029614889,3.074625172,3.117129605,3.214232374,3.300788279,3.33118689,...
3.361507951,3.389174372,3.414630616,3.438216903,3.460199618,3.48079135,...
3.500164545,3.567982556,3.694592725,3.787139494,3.921526466,4.097643229];

S1=1./V1;
S2=1./V2;
beta1=(1/Vs30-1/0.618)/(1/2.780-1/0.618);
Nh=length(H);
S(Nh)=zeros();
for i=1:1:Nh
    S(i)=(beta1)*S1(i)+(1-beta1)*S2(i);
end
V=1./S;

thick(Nh-1)=zeros();   % thickness of the layer
wtt(Nh-1)=zeros();     % wave travelling time
acct(Nh-1)=zeros();    % accumulated time
period(Nh-1)=zeros();  % period
avev(Nh-1)=zeros();    % average velocity
fn(Nh-1)=zeros();      % accumulated time

for m=2:1:Nh
    thick(m)=H(m)-H(m-1);
    wtt(m)=thick(m)/((V(m)+V(m-1))/2);
    acct(m)=sum(wtt(2:m));
    period(m)=4*acct(m);
    avev(m)=(H(m)-H(1))/acct(m);
    fn(m)=1/period(m);
end

Density(Nh)=zeros();
Vp(Nh)=zeros();
for i=1:1:Nh
    if V(i)<0.3
        Density(i)=1+(1.53*V(i)^0.85)/(0.35+1.889*V(i)^1.7);
    else
        if V(i)<3.55
            Vp(i)=0.9409+V(i)*2.0947-0.8206*V(i)^2+0.2683*V(i)^3-0.0251*V(i)^4;
            Density(i)=1.74*Vp(i)^0.25;
        else
            Vp(i)=0.9409+V(i)*2.0947-0.8206*V(i)^2+0.2683*V(i)^3-0.0251*V(i)^4;
            Density(i)=1.6612*Vp(i)-0.4721*Vp(i)^2+0.0671*Vp(i)^3-0.0043*Vp(i)^4+0.000106*Vp(i)^5;
        end
    end
end

fn=fn(2:Nh);
avev=avev(2:Nh);   
Density=Density(2:Nh);
            
amp=sqrt((beta0*roll0)./(avev.*Density.*1));

fn=fliplr(fn);
amp=fliplr(amp);

vn=log(amp);   %using the log amplification-linear frequency interpolation 
vv=interp1(fn,vn,f,'linear','extrap');
Am=exp(vv);
end

function [Am] = FUNAmfInput(f,beta0, roll0)

fn=[0.05
0.053990051
0.058298513
0.062950795
0.067974333
0.073398755
0.079256051
0.085580765
0.092410198
0.099784627
0.107747544
0.116345908
0.125630432
0.13565587
0.146481348
0.15817071
0.170792896
0.184422345
0.199139438
0.21503097
0.232190663
0.250719717
0.270727408
0.292331734
0.315660108
0.340850109
0.368050299
0.397421092
0.429135704
0.463381175
0.50035947
0.540288671
0.583404264
0.629960525
0.680232024
0.73451524
0.793130312
0.856422928
0.924766359
0.998563667
1.078250076
1.164295543
1.257207526
1.357533981
1.465866591
1.582844255
1.709156856
1.845549334
1.992826071
2.151855644
2.32357594
2.508999693
2.709220452
2.925419034
3.158870486
3.410951604
3.683149055
3.977068142
4.294442276
4.637143192
5.007191994
5.406771071
5.838236971
6.304134294
6.807210702
7.350433127
7.93700526
8.570386453
9.254312118
9.992815756
10.79025274
11.65132602
12.58111384
13.58509968
14.66920463
15.83982226
17.10385639
18.46876174
19.94258795
21.53402701
23.25246454
25.10803516
27.11168222
29.27522238
31.61141527
34.13403877
36.85797021
39.79927419
42.97529726
46.40477024
50.10791869
54.1065822
58.42434319
63.08666594
68.12104685
73.55717654
79.42711498
85.76548055
92.60965422
100
];

amp=[0.899626613
0.986182296
1.072737979
1.159293662
1.245849345
1.332405028
1.418960711
1.506727319
1.59092177
1.667360037
1.732206785
1.786999892
1.834491807
1.891614397
1.970394605
2.070359402
2.168144816
2.255759594
2.337522731
2.417521455
2.493021774
2.560664197
2.603359368
2.60805883
2.573748035
2.516011608
2.439677917
2.347329362
2.244305143
2.16245732
2.075275315
1.980700642
1.88340598
1.786607021
1.704272072
1.648820619
1.613797682
1.590538323
1.547651381
1.551104214
1.567480606
1.580528226
1.614417765
1.65288853
1.689034543
1.730297724
1.814994141
1.907910622
1.975675339
2.03435957
2.095554793
2.110043805
2.117043669
2.123643655
2.117790136
2.087002409
2.052156891
2.001607586
1.925669035
1.847148346
1.790262898
1.737407844
1.700875991
1.673659646
1.629151615
1.576543925
1.521090673
1.487583699
1.471437404
1.539927451
1.611673636
1.667701617
1.719211946
1.77546836
1.823995825
1.871637044
1.933481052
1.993563668
1.97758073
1.966004511
1.946517037
1.889425345
1.826480829
1.778032793
1.733182886
1.677598785
1.61491813
1.556502308
1.50064856
1.446920777
1.432246161
1.41454942
1.399123935
1.381520164
1.369120903
1.348286848
1.327452793
1.306618739
1.285784684
1.264950629
];

vn=log(amp);   %using the log amplification-linear frequency interpolation 
vv=interp1(fn,vn,f,'linear','extrap');
Am=exp(vv);

end

function [Tpath]=FUNTpath(R)

if R<=50
    Tpath=0.2*R;
elseif R<=90
    Tpath=3*(R-50)/40+10;
elseif R<=250
    Tpath=-0.5*(R-90)/160+13;
elseif R<320
    Tpath=4.5*(R-250)/70+12.5;
else
    Tpath=0.145*(R-320)+17;
end

end
