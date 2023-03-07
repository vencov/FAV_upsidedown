close all;
clear all;


[y,fs] = audioread('1z67931a_44kHz.wav');

sig = [];
for i=1:3
    sig = [sig; y];
    
    
end



dT = 0.3; % gat duration

NdT = round(dT*fs);
Nrep = NdT*4;

mask01 = [ones(NdT,1); zeros(NdT,1);];
mask02 = [zeros(NdT,1); tukeywin(NdT,0.1).*randn(NdT,1);];

Nwin = length(mask01);
sigGap = zeros(size(sig));
sigGapN = zeros(size(sig));
for i=1:floor(length(sig)/Nwin)
    sigGap(1+(i-1)*Nwin:Nwin*i) = sig(1+(i-1)*Nwin:Nwin*i).*mask01;
    sigGapN(1+(i-1)*Nwin:Nwin*i) = sig(1+(i-1)*Nwin:Nwin*i).*mask01+0.2*mask02;
end;

% sound(0.01*sigGap,fs)