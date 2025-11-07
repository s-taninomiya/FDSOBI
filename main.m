clear; close all; clc;

% create source signals
T = 2;
fs = 8000;
timeAxis = (0 : 1/fs : T)';
srcSig1 = sin(2 * pi * 100 * timeAxis);
srcSig2 = sin(2 * pi * 500 * timeAxis);
srcSig3 = sin(2 * pi * 1000 * timeAxis);
srcSig = [srcSig1 srcSig2 srcSig3];
[signalLength, nSrc] = size(srcSig);

% create observed signals based on mixing model
nMic = nSrc;
A = randn(nSrc, nMic);
obsSig = (A * srcSig')';

% define some constants for FDSOBI
tauList = 2:2:20;
windowSize = 4096;
shiftSize = windowSize / 8;
fftNum = windowSize;
windowType = 'b';
nIter = 50;
drawConv = true;

% apply FDSOBI (Frequency-Domain Second-Order Blind Source Identification)
[estSig, cost] = FDSOBI(obsSig, nSrc, tauList, fs, windowSize, shiftSize, fftNum, windowType, nIter, drawConv);

F = DGTtool("windowLength", windowSize, "windowShift", shiftSize, "FFTnum", fftNum, "windowName", windowType);

srcSpec = F(srcSig);
obsSpec = F(obsSig);
estSpec = F(estSig);

freqAxis = (0:fftNum/2-1) * fs / fftNum;
timeAxis = (0:size(srcSpec,2)-1) * shiftSize / fs;

% normalize for visual comparison
srcSpec_dB = 20*log10(abs(srcSpec(:,:,1:nSrc))/max(abs(srcSpec(:))) + 1e-12);
obsSpec_dB = 20*log10(abs(obsSpec(:,:,1:nSrc))/max(abs(obsSpec(:))) + 1e-12);
estSpec_dB = 20*log10(abs(estSpec(:,:,1:nSrc))/max(abs(estSpec(:))) + 1e-12);

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
