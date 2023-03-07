function y = genLinSweptSine(f1,f2,Nsamp,fsamp)
% generate Linear Swept Sine with given parameters
% f1 - starting freq in Hz
% f2 - final freq in Hz
% Nsamp - number of samples (choose to be in relation to sampl freq and f1
% and f2, e.g. 4800 samples which for 96 kHz sampling freq gives a chirp with 0.05 sec duration)
% fsamp - sampling freq (96 kHz)
% addZeros - add zero samples behind the chirp, the number of samples

T = Nsamp/fsamp;
%T_ = (Nsamp-1)/fsamp;
tx = 0:1/fsamp:(Nsamp-1)/fsamp;

A = 1;


y = zeros(length(tx),1);

for k=1:length(tx)
    y(k) = A*sin(2*pi*(f1+(f2-f1)/(2*(T))*tx(k))*tx(k));
end

%y = [y; zeros(addZeros,1)];
    