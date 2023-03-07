


f1 = 200;
f2 = 1500;
dT = 1;
fs = 44.1e3;
Nsamp = round(dT*fs);

Nrep = 6;
dGapT = 0.1;
dGap = round(dGapT*fs);
zeroFrame = ones(Nsamp,1);
noiseFrame = zeroFrame;
zeroFrame(Nsamp/2-dGap:Nsamp/2+dGap) = 0;


noiseFrame = zeros(Nsamp,1);
overF = 0.01;
overFN = round(overF*fs);
% noiseFrame(Nsamp/2-dGap-overFN:Nsamp/2+overFN+dGap) = hann(2*overFN+dGap*2+1).*randn(2*overFN+dGap*2+1,1);
noiseFrame(Nsamp/2-dGap-overFN:Nsamp/2+overFN+dGap) = tukeywin(2*overFN+dGap*2+1,0.1).*randn(2*overFN+dGap*2+1,1);

sig = [];
sigGap = [];
sigGapN = [];
for i=1:Nrep
    
    sigup = genLinSweptSine(f1,f2,Nsamp,fs);
    sigdw = genLinSweptSine(f2,f1,Nsamp,fs);
    
    sig = [sig; sigup; sigdw];
    
    sigGap = [sigGap; zeroFrame.*sigup; zeroFrame.*sigdw];
    
    sigGapN = [sigGapN; (sigup+3*noiseFrame); (sigdw + 3*noiseFrame)];
    
    
    
end

%A = 0.05;
%sig = A*sig;

figure; 
plot(sigGapN)

%sound(0.01*sigGap,fs)

% sound(0.01*sigGapN,fs)

% sound(0.01*sig,fs)