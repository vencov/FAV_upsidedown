close all;
clear all;


[y,fs] = audioread('AE.wav');
sig = y(1:1190700,1);


dT = 0.1; % gat duration

NdT = round(dT*fs);
Nrep = NdT*4;

mask01 = [ones(3*NdT,1); zeros(NdT,1);];
mask02 = [zeros(3*NdT,1); tukeywin(NdT,0.1).*randn(NdT,1);];

Nwin = length(mask01);
sigGap = zeros(size(sig));
sigGapN = zeros(size(sig));
for i=1:floor(length(sig)/Nwin)
    sigGap(1+(i-1)*Nwin:Nwin*i) = sig(1+(i-1)*Nwin:Nwin*i).*mask01;
    sigGapN(1+(i-1)*Nwin:Nwin*i) = sig(1+(i-1)*Nwin:Nwin*i).*mask01+0.3*mask02;
end;

% sound(0.01*sigGap,fs)