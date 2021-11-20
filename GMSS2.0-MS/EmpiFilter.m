% The following is the crustal amplitude model
function [Am] = EmpiFilter(f)

fn=[0.01,0.02,0.03,0.04,0.05,0.1,0.2,0.5,1.0,1.5,2.0,3.0,5.0];
amp=[2.92,2.83,2.72,2.65,2.53,2.44,2.32,1.26,1.19,1.16,1.12,1.1,1.2];

vn=log(amp);   %using the log amplification-linear frequency interpolation 
vv=interp1(fn,vn,f,'linear','extrap');
Am=exp(vv);
end