clear; close all; clc;

addpath("bss_eval/");

% create source signals
T = 2;
fs = 8000;
timeAxis = (0 : 1/fs : T)';
srcSig1 = sin(2 * pi * 500 * timeAxis);
srcSig2 = sin(2 * pi * 750 * timeAxis);
srcSig3 = sin(2 * pi * 1000 * timeAxis);
srcSig = [srcSig1 srcSig2 srcSig3];
[signalLength, nSrc] = size(srcSig);

% create observed signals based on mixing model
nMic = nSrc;
A = randn(nSrc, nMic);
obsSig = (A * srcSig')';

% define some constants for FDSOBI
tauList = 2:2:20;
windowSize = 1024;
shiftSize = windowSize / 8;
fftNum = windowSize;
windowType = 'b';
nIter = 20;
permSolver = "IPS";
isFilt = false;
drawConv = true;

% apply FDSOBI (Frequency-Domain Second-Order Blind Source Identification)
[estSig, cost] = FDSOBI(obsSig, srcSig, nSrc, tauList, fs, windowSize, shiftSize, fftNum, windowType, nIter, permSolver, isFilt, drawConv);

% calcurate SDR
[inSdr, inSir, inSar] = bss_eval_sources(obsSig.', srcSig.');
[outSdr, outSir, outSar] = bss_eval_sources(estSig.', srcSig.');
inSdr1 = inSdr(1, 1); inSdr2 = inSdr(2, 1); inSdr3 = inSdr(3, 1);
inSir1 = inSir(1, 1); inSir2 = inSir(2, 1); inSir3 = inSir(3, 1);
inSar1 = inSar(1, 1); inSar2 = inSar(2, 1); inSar3 = inSar(3, 1);
outSdr1 = outSdr(1, 1); outSdr2 = outSdr(2, 1); outSdr3 = outSdr(3, 1);
outSir1 = outSir(1, 1); outSir2 = outSir(2, 1); outSir3 = outSir(3, 1); 
outSar1 = outSar(1, 1); outSar2 = outSar(2, 1); outSar3 = outSar(3, 1);
impSdr1 = outSdr1 - inSdr1; impSdr2 = outSdr2 - inSdr2; impSdr3 = outSdr3 - inSdr3;
impSir1 = outSir1 - inSir1; impSir2 = outSir2 - inSir2; impSir3 = outSir3 - inSir3;
impSar1 = outSar1 - inSar1; impSar2 = outSar2 - inSar2; impSar3 = outSar3 - inSar3;
disp("impSDR");
disp(impSdr1);
disp(impSdr2);
disp(impSdr3);

% plot spectrograms
F = DGTtool("windowLength", windowSize, "windowShift", shiftSize, "FFTnum", fftNum, "windowName", windowType);
srcSpecgram = F(srcSig);
obsSpecgram = F(obsSig);
estSpecgram = F(estSig);

freqAxis = (0:fftNum/2-1) * fs / fftNum;
timeAxis = (0:size(srcSpecgram,2)-1) * shiftSize / fs;

srcSpec_dB = 20*log10(abs(srcSpecgram(:,:,1:nSrc))/max(abs(srcSpecgram(:))) + 1e-12);
obsSpec_dB = 20*log10(abs(obsSpecgram(:,:,1:nSrc))/max(abs(obsSpecgram(:))) + 1e-12);
estSpec_dB = 20*log10(abs(estSpecgram(:,:,1:nSrc))/max(abs(estSpecgram(:))) + 1e-12);

figure('Name','Spectrogram Comparison','Position',[100 100 1200 800]);

for i = 1:nSrc
    subplot(nSrc,3,3*(i-1)+1);
    imagesc(timeAxis, freqAxis, srcSpec_dB(1:fftNum/2,:,i));
    axis xy; colormap jet; colorbar;
    title(['Source ', num2str(i)]);
    ylabel('Frequency [Hz]');
    if i == nSrc, xlabel('Time [s]'); end

    subplot(nSrc,3,3*(i-1)+2);
    imagesc(timeAxis, freqAxis, obsSpec_dB(1:fftNum/2,:,i));
    axis xy; colormap jet; colorbar;
    title(['Observed (Ch ', num2str(i), ')']);
    if i == nSrc, xlabel('Time [s]'); end

    subplot(nSrc,3,3*(i-1)+3);
    imagesc(timeAxis, freqAxis, estSpec_dB(1:fftNum/2,:,i));
    axis xy; colormap jet; colorbar;
    title(['Estimated ', num2str(i)]);
    if i == nSrc, xlabel('Time [s]'); end
end
sgtitle('Spectrograms of Sources, Observations, and Estimates');
