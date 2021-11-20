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
