function [l1,l2,l3] = eig3x3symm(a,b,c,d,e,f)
%
% Compute eigenvalues of real, symmetric (or complex, Hermitian) 3x3 matrix:
% [a  b  c]
% [b* d  e]
% [c* e* f]
%
% See OK Smith (1961). Communications of the ACM 4(4): p. 168.
% Note this is similar but slightly different than the wikipedia code

% 3m = tr(A)
m = (a+d+f)./3;

% pB = A-qI = [a-q b   c  ]
%             [b*  d-q e  ]
%             [c*  e*  f-q]
% 2q = det(pB) = (a-q)[(d-q)(f-q)-|e|^2]-b[b*(f-q)-c*e]+c[b*e*-c*(d-q)]
% = (a-q)(d-q)(f-q)-(a-q)|e|^2-(f-q)|b|^2-(d-q)|c|^2+2Re{bc*e}
q = ((a-m).*(d-m).*(f-m) - (a-m).*abs(e).^2 - (f-m).*abs(b).^2 - (d-m).*abs(c).^2)./2 + real(b.*conj(c).*e);

% 6p = ||pB||_F^2 = tr(pB^2)
p = ((a-m).^2+(d-m).^2+(f-m).^2+2.*(abs(b).^2+abs(c).^2+abs(e).^2))./6;
sqrtp = sqrt(p);

% eigenvalues from hyperbolic formula
phi = max(0,p.^3 - q.^2); % should be nonnegative
phi = mod(atan2(sqrt(phi),q),2*pi)./3;

cosphi = cos(phi);
sinphi = sin(phi);

l1 = m + 2.*sqrtp.*cosphi;
l2 = m - sqrtp.*(cosphi-sqrt(3).*sinphi);
l3 = m - sqrtp.*(cosphi+sqrt(3).*sinphi);

end
