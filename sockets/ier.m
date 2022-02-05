function ier(portnumber)
%
%  Waits at socket portnumber (between 5000 and 5010) until a matrix arrives
% 
%  This file file is for the incompressible elastic rod (fic) problem
%

figure(1);
clg;
if ( nargin == 0 ) portnumber = 5001; end;
for i=0:10000
    m = receive(portnumber);
    if ( sum(size(m)) ~= 0 ) 
       n = size(m,1)
       nint = (n-3)/15
       h = 1/nint ;
       h1 = 0.25*h ;
       h2 = 0.50*h ;
       for i = 1: nint
           k = 3*(i-1) + 1 ;
           v1(k) = m(15*(i-1)+1,1) ;
           v1(k+1) = m(15*(i-1)+1,1) + h1*m(15*(i-1)+2,1) + ... 
                      0.5*h1^2*m(15*(i-1)+3,1) ;
           v1(k+2) = m(15*(i-1)+1,1) + h2*m(15*(i-1)+2,1) + ...
                      0.5*h2^2*m(15*(i-1)+3,1) ;
           v2(k) = m(15*(i-1)+6,1) ;
           v2(k+1) = m(15*(i-1)+6,1) + h1*m(15*(i-1)+7,1) + ... 
                      0.5*h1^2*m(15*(i-1)+8,1) ;
           v2(k+2) = m(15*(i-1)+6,1) + h2*m(15*(i-1)+7,1) + ...
                      0.5*h2^2*m(15*(i-1)+8,1) ;
       end
       k = 3*nint +1 ;
       h1 = 0.5*h ;
       h2 = h ;
       v1(k) = m(15*(i-1)+1,1) + h1*m(15*(i-1)+2,1) + ...
                           0.5*h1^2*m(15*(i-1)+3,1) ;
       v1(k+1) = m(15*(i-1)+1,1) + h2*m(15*(i-1)+2,1) + ...
                             0.5*h2^2*m(15*(i-1)+3,1) ;
       v2(k) = m(15*(i-1)+6,1) + h1*m(15*(i-1)+7,1) + ...
                           0.5*h1^2*m(15*(i-1)+8,1) ;
       v2(k+1) = m(15*(i-1)+6,1) + h2*m(15*(i-1)+7,1) + ...
                             0.5*h2^2*m(15*(i-1)+8,1) ;
       figure (1); plot(v1,v2); pause(2.0);
    else
       break;
    end;
end;
