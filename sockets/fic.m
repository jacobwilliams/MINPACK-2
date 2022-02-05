function [x] = fic(portnumber)
%
%  Waits at socket portnumber (between 5000 and 5010) until a matrix arrives
% 
%  This file file is for the flow in a channel (fic) problem
%

figure(1);
clg;
if ( nargin == 0 ) portnumber = 5001; end;
for i=0:10000
    m = receive(portnumber);
    if ( sum(size(m)) ~= 0 ) 
       n = size(m,1)
       nint = n/8
       for i = 1: nint
           v1(i) = (i-1)/(nint+1) ;
           v2(i) = m(8*(i-1)+2,1) ;
       end
       v1(nint+1) = 1.0 ;
       v2(nint+1) = 0.0 ;
       figure (1); plot(v1,v2); pause(0.0);
    else
       break;
    end;
end;
