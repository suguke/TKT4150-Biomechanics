function r = real_ifft(Z,w,t)
%%
% Inputs:
% Z the non-negative frequency components of the fft
% w the frequencies associated with these components
% t the times at which the interpolant is to be evaluated
%
% Outputs:
% r = the time domain Trigonometric Interpolant of the 
% frequency domain signal Z corresponding to the times 
% input, t.

%%
N=(length(w)-1)*2;
Z=Z./N;
r=zeros(size(t)) + Z(1);

for i=1:(N-1)/2+1
    r = r + real(2*Z(i+1)*exp(1j*w(i+1).*t));
end

end
    
