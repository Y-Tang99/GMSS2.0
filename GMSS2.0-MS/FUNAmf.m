% The following is the crustal amplitude model
function [Am] = FUNAmf(f)

% H=[0,0.001,0.002,0.003,0.004,0.005,0.006,0.007,0.008,0.009,0.01,0.012,...
% 0.014,0.016,0.018,0.02,0.022,0.024,0.026,0.028,0.03,0.032,0.034,0.036,...
% 0.038,0.04,0.042,0.044,0.046,0.048,0.05,0.052,0.054,0.056,0.058,0.06,...
% 0.062,0.064,0.066,0.068,0.07,0.072,0.074,0.076,0.078,0.08,0.082,0.104,...
% 0.126,0.147,0.175,0.202,0.23,0.257,0.289,0.321,0.353,0.385,0.42,0.455,...
% 0.49,0.525,0.562,0.599,0.637,0.674,0.712,0.751,0.789,0.827,0.866,0.906,...
% 0.945,0.984,1.02,1.06,1.1,1.14,1.2,1.4,1.6,1.8,2,2.2,2.4,2.6,2.8,3,3.5,...
% 4,4.5,5,5.5,6,6.5,7,7.5,8,10,15,20,30,50];
% 
% % generic very hard rock
% V1=[2.768,2.7688,2.7696,2.7704,2.7712,2.772,2.7728,2.7736,2.7744,2.7752...
% 2.776,2.7776,2.7792,2.7808,2.7824,2.784,2.7856,2.7872,2.7888,2.7904,...
% 2.792,2.7936,2.7952,2.7968,2.7984,2.8,2.8016,2.8032,2.8048,2.8064,2.808,...
% 2.80956,2.81112,2.81268,2.81424,2.8158,2.81736,2.81892,2.82048,2.82204,...
% 2.8236,2.82516,2.82672,2.82828,2.82984,2.8314,2.83296,2.85004,2.86676,...
% 2.88272,2.9035,2.92344,2.9436,2.9629,2.9853,3.00686,3.06098,3.0821,3.0718,...
% 3.0941,3.1158,3.1365,3.15796,3.17942,3.20072,3.22048,3.24024,3.260835504,...
% 3.271637476,3.281964539,3.29211285,3.302087707,3.311425176,3.320409798,...
% 3.32841313,3.337002316,3.345294257,3.353309507,3.364853485,3.399786096,...
% 3.430339103,3.457516593,3.482010086,3.50431659,3.510651329,3.516529187,...
% 3.521980007,3.527062193,3.538443827,3.54833273,3.557078298,3.564919739,...
% 3.572028075,3.578529853,3.584521359,3.590077571,3.595258021,3.600110775,...
% 3.616939825,3.647720809,3.669718999,3.700949146,3.740673099];
% 
% 
% % generic rock
% V2=[0.245,0.245,0.406898105,0.454341586,0.491321564,0.522065927,0.548608647,...
% 0.572100299,0.593261253,0.612575285,0.630384474,0.662434293,0.690800008,...
% 0.716351448,0.739672767,0.761177017,0.781168059,0.799876552,0.817482113,...
% 0.834127602,0.84992869,0.872659052,0.894459078,0.915511227,0.935880677,...
% 0.955623835,0.974789894,0.993422052,1.011558478,1.02923309,1.046476183,...
% 1.063314944,1.079773879,1.095875162,1.111638937,1.127083562,1.142225823,...
% 1.157081117,1.171663602,1.185986333,1.200061379,1.21389992,1.22751234,...
% 1.240908299,1.254096801,1.26708626,1.279884545,1.409876676,1.524401504,...
% 1.623105361,1.742468929,1.822101865,1.869784562,1.91154454,1.956709713,...
% 1.998031044,2.036174042,2.071640608,2.10782397,2.141667263,2.173485515,...
% 2.203532354,2.23359926,2.262120148,2.289978882,2.315853324,2.341268645,...
% 2.366247007,2.389604645,2.412077883,2.434298272,2.456270778,2.4769581,...
% 2.496972474,2.514890987,2.534215818,2.55296508,2.571175941,2.597555276,...
% 2.678472608,2.750601065,2.815833418,2.875495547,2.930554778,2.981739972,...
% 3.029614889,3.074625172,3.117129605,3.214232374,3.300788279,3.33118689,...
% 3.361507951,3.389174372,3.414630616,3.438216903,3.460199618,3.48079135,...
% 3.500164545,3.567982556,3.694592725,3.787139494,3.921526466,4.097643229];
% 
% S1=1./V1;
% S2=1./V2;
% beta1=(1/Vs30-1/0.618)/(1/2.780-1/0.618);
% Nh=length(H);
% S(Nh)=zeros();
% for i=1:1:Nh
%     S(i)=(beta1)*S1(i)+(1-beta1)*S2(i);
% end
% V=1./S;
% 
% thick(Nh-1)=zeros();   % thickness of the layer
% wtt(Nh-1)=zeros();     % wave travelling time
% acct(Nh-1)=zeros();    % accumulated time
% period(Nh-1)=zeros();  % period
% avev(Nh-1)=zeros();    % average velocity
% fn(Nh-1)=zeros();      % accumulated time
% 
% for m=2:1:Nh
%     thick(m)=H(m)-H(m-1);
%     wtt(m)=thick(m)/((V(m)+V(m-1))/2);
%     acct(m)=sum(wtt(2:m));
%     period(m)=4*acct(m);
%     avev(m)=(H(m)-H(1))/acct(m);
%     fn(m)=1/period(m);
% end
% 
% Density(Nh)=zeros();
% Vp(Nh)=zeros();
% for i=1:1:Nh
%     if V(i)<0.3
%         Density(i)=1+(1.53*V(i)^0.85)/(0.35+1.889*V(i)^1.7);
%     else
%         if V(i)<3.55
%             Vp(i)=0.9409+V(i)*2.0947-0.8206*V(i)^2+0.2683*V(i)^3-0.0251*V(i)^4;
%             Density(i)=1.74*Vp(i)^0.25;
%         else
%             Vp(i)=0.9409+V(i)*2.0947-0.8206*V(i)^2+0.2683*V(i)^3-0.0251*V(i)^4;
%             Density(i)=1.6612*Vp(i)-0.4721*Vp(i)^2+0.0671*Vp(i)^3-0.0043*Vp(i)^4+0.000106*Vp(i)^5;
%         end
%     end
% end
% 
% fn=fn(2:Nh);
% avev=avev(2:Nh);   
% Density=Density(2:Nh);
%             
% amp=sqrt((beta0*roll0)./(avev.*Density.*1));
% 
% fn=fliplr(fn);
% amp=fliplr(amp);

fn=[0.00 0.10 0.24 0.45 0.79 1.38 1.93 2.85 4.03 6.34 12.54 21.23 33.39 82.0];
amp=[1.45 1.45 1.87 2.17 2.48 2.91 3.20 3.36 3.35 3.13 2.95 2.99 3.11 3.23];

% fn=[0.0001 0.1 0.24 0.45 0.79 1.38 1.93 2.85 4.03 6.34 12.5 21.2 33.4 82];
% amp=[1 1.07 1.15 1.24 1.39 1.67 1.88 2.08 2.2 2.31 2.41 2.45 2.47 2.5];

vn=log(amp);   %using the log amplification-linear frequency interpolation 
vv=interp1(fn,vn,f,'linear','extrap');
Am=exp(vv);
end