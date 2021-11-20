% The folowing is the seismological model defining the Fourier Amplitude
% Spectrum (FAS)

function [subFAS] = InputFAS(subf,subR,C,subM0,subf0,~,kappa)    
% This function gives the Fourier amplitude spectrum

% source spectrum
w=(2*subf.*pi).^2;
E1=C*w.*subM0;    
E2=1+(subf./subf0).^2;
E=E1./E2;
 
% geometric spreading
if subR<=40
    G=subR^(-1);
else
    G=40^(-1)*(subR/40)^(-0.5);
end

% if subR<=70
%     G=subR^(-1.3);
% else
%     if subR<=140
%         G=(70^(-0.2)/70^(1.3))*(subR^0.2);
%     else
%         G=(70^(-0.2)/70^(1.3))*(140^(0.5)/140^(-0.2))*(subR^(-0.5));
%     end
% end
   
% anelastic attenuation 
% Q0=180;
% nq=0.45;
% Qmin=60;

cq=3.7;        % note this is different from beta0!!
Q0=187;
nq=0.55;
Qmin=150;
Q=max(Qmin,Q0*(subf.^(nq)));

Ae1=-pi*subf.*subR;
Ae2=Q.*cq;
Ae=exp(Ae1./Ae2);

% crustal amplification

Am=FUNAmf(subf);
Filter=EmpiFilter(subf);

% Note: Users can add any amplification factors, including site
% amplification, empirical amplification, etc.

% High-frequency attenuation

% note kappa=kappa0+dk*(M-Mkappa), in this prpgram dk is set to be 0, thus
% kappa=kappa0. Users can change it if needed. 

An=exp(-pi*subf.*kappa);

% Butterworth low-cut (high-pass) filter, which allows the signals
% at high frequencies and remove the signals at low frequencies

% final subFAS

subFAS=E.*Ae.*Am.*An.*Filter.*G;     % unit cm/s;

end