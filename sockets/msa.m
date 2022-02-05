function msa(portnumber)
%
%  Waits at socket portnumber (between 5000 and 5010) until a matrix arrives
% 
%  This file is for the minimal surface area (msa) problem
%

figure(1);
clg;
if ( nargin == 0 ) portnumber = 5001; end;
for i=0:10000
    m = receive(portnumber);
    if ( sum(size(m)) ~= 0 ) 
       n = size(m,1)
       nx = sqrt(n)
       ny = n/nx 
       for j = 1: nx
           for i = 1:ny
               v(i,j) = m((i-1)*nx+j,1) ;
           end
       end

%      Surface plots

       figure(1); surf(v);  pause (0.0);

   else
      break;
   end;
end;
