function [estSig, cost] = FDSOBI(mixSig, nSrc, tauList, sampFreq, windowSize, shiftSize, fftNum, windowType, nIter, drawConv)
% Frequency-Domain Second-Order Blind Source Identification (FDSOBI)
%
% Coded by S. Taninomiya
%
% [syntax]
%   [estSig, cost] = FDSOBI(mixSig, nSrc, tauList, sampFreq, windowSize, shiftSize, fftNum, windowType, nIter, drawConv)
%
% [Inputs]
%       mixSig: observed Mixture (sigLen x nMic)
%         nSrc: number of sources in the mixture (scalar)
%      tauList: delay length [point] to obtain delayed SCM (nTau x 1)
%     sampFreq: sampling frequency [Hz] of mixSig (scalar)
%   windowSize: window length [points] in STFT (scalar, default: next higher power of 2 that exceeds 0.256*sampFreq)
%    shiftSize: shift length [points] in STFT (scalar, default: fftSize/2)
%       fftNum: fft length [points] in STFT (scalar, default: windowsize)
%   windowType: window function used in STFT (name of window function, default: 'blackman')
%        nIter: number of iterations in the joint diagonalization (scalar, default: 10)
%     drawConv: plot cost function values in each iteration or not (true or false, default: false)
%
% [outputs]
%       estSig: estimated signals (sigLen x nMic x nSrc)
%         cost: convergence behavior of cost function in ILRMA (nIter+1 x 1)
%

% Arguments check and set default values
arguments
    mixSig (:, :) double
    nSrc (1,1) double {mustBeInteger(nSrc)}
    tauList (:,1) double {mustBeInteger(tauList)}
    sampFreq (1,1) double
    windowSize (1,1) double {mustBeInteger(windowSize)} = 2^nextpow2(0.256*sampFreq)
    shiftSize (1,1) double {mustBeInteger(shiftSize)} = windowSize/2
    fftNum (1,1) double {mustBeInteger(fftNum)} = windowSize
    windowType char {mustBeMember(windowType,{'h','b'})} = 'b'
    nIter (1,1) double {mustBeInteger(nIter)} = 10
    drawConv (1,1) logical = false
end

% Error check
[sigLen, nMic] = size(mixSig); % sigLen: signal length, nMic: number of channels
if sigLen < nMic; error("The size of mixSig might be wrong.\n"); end
if nMic < nSrc || nSrc < 2; error("The number of channels must be equal to or grater than the number of sources in the mixture.\n"); end
if sampFreq <= 0; error("The sampling frequency (sampFreq) must be a positive value.\n"); end
if windowSize < 1; error("The FFT length in STFT (fftSize) must be a positive integer value.\n"); end
if shiftSize < 1; error("The shift length in STFT (shiftSize) must be a positive integer value.\n"); end
if nIter < 1; error("The number of iterations (nIter) must be a positive integer value.\n"); end

% Apply multichannel short-time Fourier transform
F = DGTtool("windowLength", sindowSize, "sindowShift", shiftSize, "FFTnum", fftNum, "windowName", windowType);
mixSpecgram = F(mixSig);

% Apply FDSOBI
[estSpecgram, demixMat, cost] = local_FDSOBI(maxSpecgram, nIter, tauList, drawConv);

% Apply back projection
[I, J, ~] = size(estSpecgram);
candEstSpecgram = zeros(I, J, nSrc, nMic);
candDemixMat = zeros(nSrc, nMic, I, nMic);
for iMic = 1 : nMic
    [candEstSpecgram(:,:,:,iMic), candDemixMat(:,:,:,iMic)] = local_backProjectionInit(estSpecgram, mixSpecgram(:,:,iMic), demixMat); % scale-fixed estimated signal
end
ind = (1:nMic).';
estSpecgramFix = zeros(I, J, nMic);
demixMatFix = zeros(nMic, nSrc, I);
for iMic = 1 : nMic
    estSpecgramFix(:, :, iMic) = candEstSpecgram(:, :, ind(iMic), iMic);
    demixMatFix(iMic, :, :) = candDemixMat(ind(iMic), :, :, iMic);
end

% Calculate estimated time-domain signal
estSig = F.pinv(estSpecgramFix);

end

%% Local function for FDSOBI
function [Y, W, cost] = local_FDSOBI(X, nIter, tauList, drawConv)
% [inputs]
%          X: observed multichannel spectrograms (I x J x M)
%      nIter: number of iterations of the parameter updates
%    tauList: delay length [point] to obtain delayed SCM (nTau x 1)
%   drawConv: plot cost function values in each iteration or not (true or false)
%
% [outputs]
%          Y: estimated spectrograms of sources (I x J x N)
%          W: demixing matrices (N x M x I)
%       cost: convergence behavior of cost function in ILRMA (nIter+1 x 1)
%
% [scalars]
%
%
% [matrices]
%
%

% Initialization
[I, J, M] = size(X); % I:frequency bins, J: time frames, M: channels
nTau = size(tauList, 1);
pX = permute(X, [3, 2, 1]); % permuted X whose dimensions are M x J x I
N = M; % N: number of sources, which equals to M in ILRMA
V = zeros(N, M, I); % frequency-wise whitening matrix
U = zeros(N, N, I); % frequency-wise Givens-rotation matrix
W = zeros(N, M, I); % frequency-wise demixing matrix
Y = zeros(I, J, N); % estimated spectrograms of sources (Y(i,:,n) =  W(n,:,i)*pX(:,:,i))
for i = 1:I
    U(:, :, i) = eye(N); % initial Givens-rotation matrices are set to identity matrices
end
cost = zeros(nIter+1, 1);

% whitening SCM
for i = 1:I
    RX0 = cov(pX(:, :, i)');
end

end