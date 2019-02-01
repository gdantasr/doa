function S = steeredResponseDelayAndSum(R, e, w)
%steeredResponseDelayAndSum - calculate delay and sum in frequency domain
%
%Calculates the steered response from the delay-and-sum algorithm in the
%frequency domain based on sensor positions, input signal and scanning angles
%
%S = steeredResponseDelayAndSum(R, e, w)
%
%IN
%R - PxP correlation matrix / cross spectral matrix (CSM)
%e - NxMxP steering vector/matrix for a certain frequency
%w - 1xP vector of element weights
%
%OUT
%S - NxM matrix of delay-and-sum steered response power
%


[nPointsY, nPointsX, nMics] = size(e);

%Make the weighting vector a column vector instead of row vector
if isrow(w)
    w = w';
end

%Delay and sum steered response power
S = zeros(nPointsY, nPointsX);
for pointY = 1:nPointsY
    for pointX = 1:nPointsX
        ee = reshape(e(pointY, pointX, :), nMics, 1);
        %S(pointY, pointX) = (w.*ee)'*R*(ee.*w);
        S(pointY, pointX) = 1/((w.*ee)'*inv(R + 0.1*eye(length(R)))*(ee.*w));
    end
end



