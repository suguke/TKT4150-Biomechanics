function Z = ZWK(omega, R, C)
% calculate characteristic impedance of the two element WK model.
Z = R./(1 + 1i.*omega*R*C);