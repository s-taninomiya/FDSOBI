clear; close all; clc;

% create source signals
fs = 8000;
timeAxis = (0 : 1/fs : 2)';
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
tauList = 10:10:100;
windowSize = 4096;
shiftSize = windowSize / 8;
fftNum = windowSize;
windowType = 'b';
nIter = 50;
drawConv = true;

% apply FDSOBI (Frequency-Domain Second-Order Blind Source Identification)
[estSig, cost] = FDSOBI(obsSig, nSrc, tauList, fs, windowSize, shiftSize, fftNum, windowType, nIter, drawConv);