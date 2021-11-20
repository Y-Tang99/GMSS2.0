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