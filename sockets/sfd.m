function sfd(portnumber)
%
%  Waits at socket portnumber (between 5000 and 5010) until a matrix arrives
% 
%  This file is for the swirling flow between disks (sfd) problem
%

figure(1);
clg;
figure(2);
clg;
if ( nargin == 0 ) portnumber = 5001; end;
for i=0:10000
    m = receive(portnumber);
    if ( sum(size(m)) ~= 0 ) 
       n = size(m,1)
       nint = n/14
       for i = 1: nint
           v1(i) = -1.0 + 2*(i-1)/(nint+1) ;

%          For the radial velocity f'.

           v2(i) = m(14*(i-1)+2) ;

%          For the angular velocity g

           v3(i) = m(14*(i-1)+9) ;
       end
       v1(nint+1) = 1.0 ;

%      Plot of the radial velocity f'.

       v2(nint+1) = 0.0 ;
       figure (1); plot(v1,v2); pause(0.0);

%      Plot of the angular velocity g.

       v3(nint+1) = 1.0 ;
       figure (2); plot(v1,v3); pause(0.0);
    else
       break;
    end;
end;
