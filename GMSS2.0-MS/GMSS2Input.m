%--------------------------------------------------------------------------
clear all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       ^^ INPUT PARAMETERS ^^
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------
% (1) General Input Parameters
%--------------------------------------------------------------------------

Nhyp=1;            % number of hypocenters
NS=1;             % number of simulations

M=6;               % Moment magnitude
stress=150;        % Stress drop, bars
Vs30=0.76;         % time-averaged shear wave velocity over top 30 m, km/s
kappa=0.04;       % high-frequency decay slope, s

beta0=3.7;         % Source shear wave velocity (SWV), km/s
stress_ref=70;     % reference stress drop, used for finding FL & FW
roll0=2.8;         % Source density, g/cm^3 
y=0.8;             % ratio of rupture velocity and Source SWV
vrup=y*beta0;      % rupture velocity
z=1.68;            % corresponding to x = 0.5
             
dt=0.002;          % Time step, no larger than 0.02 s 

%%%---------------------------------------------------------------------%%%
% Note: To be more strict, the time pad should be computed using the following 
% function: tpad=1.5*n*f0
% where n is the order of Butterworth filter, usually using 4; and f0 is the
% corner-frequency (Boore, BSSA, 2005).
%%%---------------------------------------------------------------------%%%
tpadl=40;          % Time pad before simulated series
tpadt=40;          % Time pad after simulated series
npadl=tpadl/dt;
npadt=tpadt/dt;


%--------------------------------------------------------------------------
% (2) Finite-fault Input Parameters
%--------------------------------------------------------------------------

FaultLat=23.61719;            % latitude of upper edge of fault
FaultLon=120.6885;            % longitude of upper edge of fault
Fstrike=3;                    % fault strike,degree (°)
Fdip=29;                      % fault dip, degree (°)
Rake=45;                      % rake angle, degree (°)

% Fault Dimensions

%FDfactor=1;  % use customdefined fault dimensions
% FL0=55;       % fault length
% FW0=40;       % fault width

% FDfactor=2; % use relation dveloped by Wells & Coppersmith (1994)
% FDfactor=3; % use relation dveloped by Leonard (2010)
 FDfactor=4; % use relation dveloped by Kumar et al. (2017)
% FDfactor=5; % use relation dveloped by Cheng et al. (2019)



% The following parameters are needed to locate the origin point


h_ref=0.94;    % fault depth to upper edge
h_min=3.0;     % Campbell depth to seismogenic region, usually set as 3.0
% Subfault Dimension
dl=2;                      % subfault length, no less than 1.5 km
dw=2;                      % subfault width, no less than 1.5 km
% pulsing percentage
pp=0.5;                    

%--------------------------------------------------------------------------
% (3) Site Input Parameters
%--------------------------------------------------------------------------

NSL=8;           % number of sites
% Site location

% Two options for determining site location: lattitude & longitude (LatLon),
% and distance & azimuth (DistAz). Users need to choose one for their purposes.

 SLIndex='LatLon';         % Latitude and longitude
% SLIndex='DistAz';          % Distance and Azimuth

%    SL1=5;                    % Input values to get site location
%    SL2=0;

SL1=[24 23.84 23.83333 24.4 24.2508 23.3927 21.9256 23.5894];
SL2=[120.8333 121.04 121.2 120.6 121.5196 122.0395 121.0255 117.7443];

%--------------------------------------------------------------------------
% (4) Input Parameters for Fourier domain and repsonse spectra Output
%--------------------------------------------------------------------------

fomin=0.1;        % lower output frequency
fomax=100;        % upper output frequency

Ntn=100;          % number of natural periods
psi=0.05;         % Viscous damping value


%--------------------------------------------------------------------------
% (5) Input Parameters for selecting the simulation methods
%--------------------------------------------------------------------------

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







