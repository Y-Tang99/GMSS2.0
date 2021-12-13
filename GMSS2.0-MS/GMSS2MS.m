%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Ground Motion Simulation System (GMSS) Version 2.0 -- Stochastic Finite
% Fault Source Simulation
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

run GMSS2Input.m

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
    s1f=-FL/2;         % Along strike near edge
    s2f=FL/2;        % Along strike far edge 
    w1f=-FW/2;         % Down dip near edge
    w2f=FW/2;        % Down dip far edge
else
    s1f=0;         % Along strike near edge
    s2f=FL;        % Along strike far edge 
    w1f=0;         % Down dip near edge
    w2f=FW;        % Down dip far edge
end

nl=round(FL/dl);           % number of subfaults along strike
nw=round(FW/dw);           % number of subfaults along dip
if nl<=1
    nl=1;
    dl=FL;
end
if nw<=1
    nw=1;
    dw=FW;
end

NF=nl*nw;                  % total number of subfaults

fo=logspace(log10(fomin),log10(fomax),Ntn);   % Output frequency
Tn=1./fo;
Tn=Tn(end:-1:1);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       ^^ MAIN PROGRAM ^^
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--------------------------------------------------------------------------
% (1) Dtermination of Various Distances 
%--------------------------------------------------------------------------

% R: Distance with respect to the origin point
% Az: azimuth
% Rrup: Closest distance to fault rupture
% Rjb: Closet distance to the surface projection of the fault rupture
% Rseis: Closest distance to rupture portion in seismogenic crust
% Rx: Horizontal distance perpendicular to the rupture strike (site
% coordinate, positive on hanging wall side)
% Ry0: A distance along strike, used by Abrahamson et al. (2014)to taper
% their hanging wall effect.  It can only be greater than or equal to 0.0.

SiteLat(NSL)=zeros();
SiteLon(NSL)=zeros();
R(NSL)=zeros();
Az(NSL)=zeros();
Rrup(NSL)=zeros();
Rjb(NSL)=zeros();
Rseis(NSL)=zeros();
Rx(NSL)=zeros();
Ry0(NSL)=zeros();
Azrup(NSL)=zeros();
Azjb(NSL)=zeros();
Azseis(NSL)=zeros();

for i=1:1:NSL   
    [SiteLat(i),SiteLon(i),R(i),Az(i)]=FUNSL(SLIndex,SL1(i),SL2(i),FaultLat,FaultLon);
    [Rrup(i),Rjb(i),Rseis(i),Rx(i),Ry0(i),Azrup(i),Azjb(i),Azseis(i)]=FUNDist(SiteLat(i),SiteLon(i),FaultLat,FaultLon,Fstrike,Fdip,h_ref,s1f,s2f,w1f,w2f,h_min);
end

%%%-------------------------!!!Warning!!!-------------------------------%%%
NT=2^16;        % This is required for pre-determining the time-step number
% Please change this number if the Error is shown similar with 
% "Unable to perform assignment because the size of the left side is 1-by-16384 and the size of the right side is 1-by-32768."
%%%---------------------------------------------------------------------%%%
t=(1:1:NT)*dt;

At(NSL,NT)=zeros();
Vt(NSL,NT)=zeros();
Dt(NSL,NT)=zeros();
Husid(NSL,NT)=zeros();
PGA(NSL)=zeros();
PGV(NSL)=zeros();
AI(NSL)=zeros();
Td5_75(NSL)=zeros();
Td5_95(NSL)=zeros();
Tm(NSL)=zeros();
FAS(NSL,Ntn)=zeros();
SA(NSL,Ntn)=zeros();
SV(NSL,Ntn)=zeros();
SD(NSL,Ntn)=zeros();

parfor ii=1:NSL
    [At(ii,:),Vt(ii,:),Dt(ii,:),Husid(ii,:),PGA(ii),PGV(ii),AI(ii),Td5_75(ii),Td5_95(ii),Tm(ii),FAS(ii,:),SA(ii,:),SV(ii,:),SD(ii,:)]=GMSS2(beta0,dl,dt,dw,Fdip,Fstrike,fo,h_ref,kappa,M,Nhyp,npadl,npadt,NS,Ntn,pp,psi,roll0,stress,Tn,tpadl,tpadt,vrup,Vs30,y,z,nl,nw,R(ii),Az(ii),NF,f0factor,BLfactor);
end

figure(1)
subplot(311)
plot(t,At(1,:))
subplot(312)
plot(t,Vt(1,:))
subplot(313)
plot(t,Dt(1,:))

figure(2)
loglog(R,PGA,R,PGV,R,AI,'Linestyle','none','Marker','*')
legend('PGA (g)','PGV (cm/s)','AI (cm/s)')

figure(3)
loglog(fo,FAS(1,:))

figure(4)
loglog(Tn,SA(1,:))

% Save to Output

Output.Dist_Az.Repi=R;
Output.Dist_Az.Az=Az;
Output.Dist_Az.Rrup=Rrup;
Output.Dist_Az.Rjb=Rjb;
Output.Dist_Az.Rseis=Rseis;
Output.Dist_Az.Rx=Rx;
Output.Dist_Az.Ry0=Ry0;

Output.T=t';
Output.F=fo';
Output.Tn=Tn';
Output.At=At';
Output.Vt=Vt';
Output.Dt=Dt';
Output.Husid=Husid';
Output.PGA=PGA;
Output.PGV=PGV;
Output.AI=AI;
Output.Td5_75=Td5_75;
Output.Td5_95=Td5_95;
Output.Tm=Tm;
Output.FAS=FAS';
Output.SA=SA';
Output.SV=SV';
Output.SD=SD';

save Output.mat


function [At11,Vt11,Dt11,meanHusid,meanPGA,meanPGV,meanAI,D5_75,D5_95,Tm,FASo,meanSa,meanSv,meanSd]=GMSS2(beta0,dl,dt,dw,Fdip,Fstrike,fo,h_ref,kappa,M,Nhyp,npadl,npadt,NS,Ntn,pp,psi,roll0,stress,Tn,tpadl,tpadt,vrup,Vs30,y,z,nl,nw,R,Az,NF,f0factor,BLfactor)
%--------------------------------------------------------------------------
% (2) Dtermination of Number of Time Steps
%--------------------------------------------------------------------------
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
            if f0factor==1
                subTrise(nn,i,j)=subRadius/(2*vrup);    % rise time
            else
                subTrise(nn,i,j)=subRadius/(vrup);     % EXSIM original rise time
            end
            subTpath(nn,i,j)=FUNTpath(subR(nn,i,j));
            subSlip(nn,i,j)=10^(-22)*subM0(nn,i,j)/(roll0*(beta0^2)*dl*dw);
            subR(nn,i,j)=FUNsubR(R,h_ref,Fdip,Fstrike,Az,dl,dw,i,j);
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

%--------------------------------------------------------------------------
% (3) Main Procedures
%--------------------------------------------------------------------------

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
        PGA(nn,m) = max(abs(At(nn,m,:)))/980;        % unit "g"
        PGV(nn,m) = max(abs(Vt(nn,m,:)));     % unit "cm/s"
        PGD(nn,m) = max(abs(Dt(nn,m,:)));     % unit "cm"
        % Arias Intensity
        AI0(nn,m,:) = cumsum(At(nn,m,:).^2)*pi*dt/(2*980);  % unit "cm/s"
        AI(nn,m) = AI0(nn,m,end);
        Husid(nn,m,:) = AI0(nn,m,:)/AI(nn,m); 
%         td575_start(nn,m,:)=t(AI0(nn,m,:) >= 0.05*AI(nn,m));
%         td575_end(nn,m,:)=t(AI0(nn,m,:) <= 0.75*AI(nn,m));
        % Response Spectra
        [Sd(nn,m,:),Sv(nn,m,:),Sa(nn,m,:)]=CDM(At(nn,m,:),psi,dt,Tn);
        % Fourier Amplitude Spectra
        FAS0(nn,m,:) = abs(fft(At(nn,m,:),NT))*dt; % Fourier amplitudes 
        FAS(nn,m,:) = FAS0(nn,m,1:NT/2);                 % single sided FAS
    end
end

%--------------------------------------------------------------------------
% (4) Calculation of Averages
%--------------------------------------------------------------------------

% Arithmetic average over number of simulations (NS), 
% and geometric average over number of hypocenters (Nhyp)

[meanPGA,meanPGV,meanAI0,meanAI,meanHusid,meanFAS,meanSa,meanSv,meanSd]=FUNCalAve(Nhyp,NS,PGA,PGV,AI0,AI,Husid,FAS,Sa,Sv,Sd);

%--------------------------------------------------------------------------
% (5) Outputs
%--------------------------------------------------------------------------

% Duration
timed1 = t(meanAI0>=0.05*meanAI & meanAI0<=0.75*meanAI);
% Output.t5_75 = [timed1(1),timed1(end)];
% Output.D5_75 = timed1(end)-timed1(1);
D5_75=timed1(end)-timed1(1);
timed2 = t(meanAI0>=0.05*meanAI & meanAI0<=0.95*meanAI);
% Output.t5_95 = [timed2(1),timed2(end)];
% Output.D5_95 = timed2(end)-timed2(1);
D5_95=timed2(end)-timed2(1);

% Mean Period
fi = f(f>0.25 & f<20);
Ci = meanFAS(f>0.25 & f<20);
Tm = ((Ci(:)'.^2)*(1./fi(:)))/(Ci(:)'*Ci(:));
% Output.Tm = Tm;

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

end


%--------------------------------------------------------------------------
%------------This is the end of MAIN PROGRAM section-----------------------
%--------------------------------------------------------------------------

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       ^^ SUBROUTINE FUNCTIONS ^^
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [SiteLat,SiteLon,R,Az] = FUNSL(SLIndex,SL1,SL2,FaultLat,FaultLon)
% This function is used to determine the site locations

if SLIndex == 'LatLon'
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
    if SLIndex == 'DistAz'
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
        if f0factor==1
             subTrise(i,j)=subRadius/(2*vrup);     % rise time
        else
            subTrise(i,j)=subRadius/vrup;     % EXSIM original rise time
        end
        subSlip(i,j)=10^(-22)*subM0(i,j)/(roll0*(beta0^2)*dl*dw);
        subR(i,j)=FUNsubR(R,h_Ref,Fdip,Fstrike,Az,dl,dw,i,j);
        subTpath(i,j)=FUNTpath(subR(i,j));
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
        subAt(i,j,1:subN(i,j))=FUNsubAt(C,subN(i,j),dt,subdur(i,j),npadl,subR(i,j),subM0(i,j),f0,subf0(i,j),NF,Vs30,kappa);
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

function [subAt] = FUNsubAt(C,subN,dt,subdur,npadl,subR,subM0,f0,subf0,NF,Vs30,kappa)
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

subFAS=InputFAS(subf,subR,C,subM0,subf0,Vs30,kappa);

%--------------------------------------------------------------------------
% % Apply a low-cut filter (high-pass)
% fcut=0.05;
% norder=8;
% blcf=1./(1+(fcut./subf).^(2*norder));
% subFAS=blcf.*subFAS;
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

function [meanPGA,meanPGV,meanAI0,meanAI,meanHusid,meanFAS,meanSa,meanSv,meanSd] = FUNCalAve(Nhyp,NS,PGA,PGV,AI0,AI,Husid,FAS,Sa,Sv,Sd);

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

end
