function Px = minvar_method (x, p)
    % Estimates power spectrum Px based on Minimum Variance Method
%     if isempty(R)
%         R = covar(x,p);
%     end

    R = covar(x,p);
    [v ,d] = eig(R);
    U = diag(inv(abs(d) + eps));
    V = abs(fft(v, 1024)).^2;
    Px = 10*log10(p) - 10*log10(V*U);
end