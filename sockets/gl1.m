function gl1(portnumber)
%
%  Waits at socket portnumber (between 5000 and 5010) until a matrix arrives
% 
%  This file is for the Ginzburg-Landau 1-D (gl2) problem
%

figure(1);
clg;
if ( nargin == 0 ) portnumber = 5001; end;
for i=0:10000
    m = receive(portnumber);
    if ( sum(size(m)) ~= 0 ) 
       n = size(m,1)
       ds = 1.0 ;
       dn = 2.2 ;
       d = ds + dn ;
       n1 = n/4 ;
       n2 = n - 2*n1 ;
       h1 = dn/n1 ;
       h2 = 2*ds/n2 ;
       for i = 1: n1
           v1(i) = -d + (i-1)*h1 ;
       end
       for i = n1+1: n1+n2
           v1(i) = -ds + (i-n1-1)*h2 ;
       end
       for i = n1+n2+1: n+1
           v1(i) = ds + (i-n1-n2-1)*h1 ;
       end
       for i = 1: n
           v2(i) = m(i,1) ;
       end
       v2(n+1) = v2(1) ;
       figure (1); plot(v1,v2); pause(0.0);
    else
       break;
    end;
end;
