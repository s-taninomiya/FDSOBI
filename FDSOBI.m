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
F = DGTtool("windowLength", windowSize, "windowShift", shiftSize, "FFTnum", fftNum, "windowName", windowType);
mixSpecgram = F(mixSig);

% Apply FDSOBI
[estSpecgram, demixMat, cost] = local_FDSOBI(mixSpecgram, nIter, tauList, drawConv);

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
N = M; % N: number of sources, which equals to M in ILRMA
nTau = size(tauList, 1);
pX = permute(X, [3, 2, 1]); % permuted X whose dimensions are M x J x I
pZ = zeros(M, J, I);
pY = zeros(N, J, I);
RZ = zeros(M, M, I, nTau);
V = zeros(N, M, I); % frequency-wise whitening matrix
U = zeros(N, N, I); % frequency-wise Givens-rotation matrix
W = zeros(N, M, I); % frequency-wise demixing matrix

for i = 1:I
    U(:, :, i) = eye(N); % initial Givens-rotation matrices are set to identity matrices
end
cost = zeros(nIter+1, 1);

% Whitening X based on diagonalization of SCM
for i = 1:I
    RX0 = cov(pX(:, :, i)');
    [E, D] = eig(RX0);
    eigs = real(diag(D));
    eigs(eigs < eps) = eps;
    V(:, :, i) = E * diag(1 ./ sqrt(eigs)) * E';
    pZ(:, :, i) = V(:, :, i) * pX(:, :, i);
end

% Calcurate delayed SCMs
for i = 1:I
    for iTau = 1:nTau
        tau = tauList(iTau);
        RZTmp = (pZ(:, 1 + tau : end, i) * pZ(:, 1 : end - tau, i)') / (J - tau);
        RZ(:, :, i, iTau) = (RZTmp + RZTmp') / 2;
    end
end

% Calculate initial cost function value
if drawConv
    cost(1, 1) = local_calcCostFunction(RZ, U);
end

% Apply joint diadiagonalization to delayed SCMs
fprintf('Iteration:    ');
for iIter = 1:nIter
    fprintf('\b\b\b\b%4d', iIter);

    %%%%% Update parameters %%%%%
    for i = 1:I
        for p = 1:(N-1)
            for q = (p+1):N
                thetaNum = 0;
                thetaDen = 0;
                for iTau = 1:nTau
                    r = U(:, :, i)' * RZ(:, :, i, iTau) * U(:, :, i);
                    delta = 2 * real(r(p, q));
                    diff = r(p, p) - r(q, q);
                    thetaNum = thetaNum + delta;
                    thetaDen = thetaDen + diff;
                end
                theta = 0.5 * atan2(thetaNum, thetaDen);
                J = eye(N);
                c = cos(theta); s = sin(theta);
                J([p, q], [p, q]) = [c, -s; s, c];
                U(:, :, i) = U(:, :, i) * J;
            end
        end
    end

    %%%%% Calculate cost function value %%%%%
    if drawConv
        cost(iIter+1,1) = local_calcCostFunction(RZ, U);
    end
end

%%%%% Calculate demixing matrices W and estimated spectrograms Y %%%%%
for i = 1:I
    W(:,:,i) = U(:,:,i)' * V(:,:,i);
    pY(:,:,i) = W(:,:,i) * pX(:,:,i);
end
Y = permute(pY, [3, 2, 1]);

% Draw convergence behavior
if drawConv
    figure; plot((0:nIter), cost);
    set(gca, 'FontName', 'Times', 'FontSize', 16);
    xlabel('Number of iterations', 'FontName', 'Arial', 'FontSize', 16);
    ylabel('Value of cost function', 'FontName', 'Arial', 'FontSize', 16);
end

fprintf(' FDSOBI done.\n');
end

%% Local function for calculating cost function value in FDSOBI
function cost = local_calcCostFunction(RZ, U)
% [inputs]
%   RZ: delayed SCMs (M x M x I x nTau)
%   U : frequency-wise rotation matrices (N x N x I)
%
% [output]
%   cost: scalar cost value (smaller = more diagonalized)

[~, ~, I, nTau] = size(RZ);
cost = 0;

for i = 1:I
    Ui = U(:,:,i);
    for iTau = 1:nTau
        r = Ui' * RZ(:,:,i,iTau) * Ui;
        cost = cost + sum(sum(abs(r - diag(diag(r))).^2));
    end
end
end

%% Local function for applying initial back projection
function [fixY, fixW] = local_backProjectionInit(Y, S, W)
% Projection back technique to fix frequency-wise scales of estimated
% spectrogram obtained by FDICA
%
% [inputs]
%      Y: estimated spectrograms (I x J x N, nFreq x nTime x nSrc)
%      S: reference channel of observed spectrogram (I x J x 1)
%         or observed multichannel spectrogram (I x J x M, nFreq x nTime x nMic)
%      W: estimated emixing matrix (N x N x I, nSrc x nMic x nFreq)
%
% [outputs]
%   fixY: scale-fixed estimated spectrograms (I x J x N)
%         or scale-fitted estimated source images (I x J x N x M)
%   fixW: scale-fixed demixing matrix (N x N x I)
%         or scale-fitted demixing matrix for source images (N x N x I x M)
%

% Projection back
Yp = permute(Y, [3, 2, 1]); % N x J x I
Sp = permute(S, [3, 2, 1]); % 1 x J x 1 or M x J x I
Wp = permute(W, [4, 1, 2, 3]); % 1 x N x N x I
Yph = pagectranspose(Yp); % J x N x I, pagewise Hermitian transpose (Yp')
YpYph = pagemtimes(Yp, Yph); % N x N x I, pagewise matrix multiplication (Yp*Yp')
YphOnYpYph = pagemrdivide(Yph, YpYph); % J x N x I, pagewise matrix right-division (Yp'/(Yp*Yp'))
A = pagemtimes(Sp, YphOnYpYph); % 1 x N x I or M x N x I, pagewise matrix multiplication (Sp * Yp'/(Yp*Yp'))
Ap = permute(A, [1, 2, 4, 3]); % M x N x 1 x I
Ypp = permute(Yp, [4, 1, 2, 3]); % 1 x N x J x I
fixY = Ap .* Ypp; % M x N x J x I, using implicit expansion
fixY = permute(fixY, [4, 3, 2, 1]); % I x J x N x M
fixW = Ap .* Wp; % M x N x N x I, using implicit expansion
fixW = permute(fixW, [2, 3, 4, 1]); % N x N x I x M
end