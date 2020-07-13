function [CONNECT_MLS,npoints]=FindNodeInteractingParticle...
    (p,spElems_corner,npoints,NN,corner,CONNECT_MLS)

 CONNECT_TEMP_N   = zeros(1,16);
%  CONNECT_TEMP     = zeros(1,4);
 
 % Node index for particle using cubic spline
 CONNECT_TEMP_N(6) = spElems_corner + ceil(spElems_corner/(NN(1)-1)) - 1;     
 npoints{CONNECT_TEMP_N(6)} = [npoints{CONNECT_TEMP_N(6)} [p;corner]];
 CONNECT_TEMP_N(7) = CONNECT_TEMP_N(6) + 1;                         
npoints{CONNECT_TEMP_N(7)} = [npoints{CONNECT_TEMP_N(7)} [p;corner]];
 CONNECT_TEMP_N(11)= CONNECT_TEMP_N(7) + NN(1);                     
npoints{CONNECT_TEMP_N(11)} = [npoints{CONNECT_TEMP_N(11)} [p;corner]];
 CONNECT_TEMP_N(10)= CONNECT_TEMP_N(6) + NN(1);                     
npoints{CONNECT_TEMP_N(10)} = [npoints{CONNECT_TEMP_N(10)} [p;corner]];

 CONNECT_TEMP_N(1) = CONNECT_TEMP_N(6) - NN(1) - 1;                 
npoints{CONNECT_TEMP_N(1)} = [npoints{CONNECT_TEMP_N(1)} [p;corner]];
 CONNECT_TEMP_N(2) = CONNECT_TEMP_N(1) + 1;                         
npoints{CONNECT_TEMP_N(2)} = [npoints{CONNECT_TEMP_N(2)} [p;corner]];
 CONNECT_TEMP_N(3) = CONNECT_TEMP_N(2) + 1;                         
npoints{CONNECT_TEMP_N(3)} = [npoints{CONNECT_TEMP_N(3)} [p;corner]];
 CONNECT_TEMP_N(4) = CONNECT_TEMP_N(3) + 1;                         
npoints{CONNECT_TEMP_N(4)} = [npoints{CONNECT_TEMP_N(4)} [p;corner]];

 CONNECT_TEMP_N(5) = CONNECT_TEMP_N(6) - 1;                         
npoints{CONNECT_TEMP_N(5)} = [npoints{CONNECT_TEMP_N(5)} [p;corner]];
 CONNECT_TEMP_N(8) = CONNECT_TEMP_N(7) + 1;                         
npoints{CONNECT_TEMP_N(8)} = [npoints{CONNECT_TEMP_N(8)} [p;corner]];

 CONNECT_TEMP_N(9) = CONNECT_TEMP_N(5) + NN(1);                     
npoints{CONNECT_TEMP_N(9)} = [npoints{CONNECT_TEMP_N(9)} [p;corner]];
 CONNECT_TEMP_N(12)= CONNECT_TEMP_N(11) + 1;                        
npoints{CONNECT_TEMP_N(12)} = [npoints{CONNECT_TEMP_N(12)} [p;corner]];

 CONNECT_TEMP_N(13)= CONNECT_TEMP_N(9) + NN(1);                     
npoints{CONNECT_TEMP_N(13)} = [npoints{CONNECT_TEMP_N(13)} [p;corner]];
 CONNECT_TEMP_N(14)= CONNECT_TEMP_N(13) + 1;                        
npoints{CONNECT_TEMP_N(14)} = [npoints{CONNECT_TEMP_N(14)} [p;corner]];
 CONNECT_TEMP_N(15)= CONNECT_TEMP_N(14) + 1;                        
npoints{CONNECT_TEMP_N(15)} = [npoints{CONNECT_TEMP_N(15)} [p;corner]];
 CONNECT_TEMP_N(16)= CONNECT_TEMP_N(15) + 1;                        
npoints{CONNECT_TEMP_N(16)} = [npoints{CONNECT_TEMP_N(16)} [p;corner]];

 CONNECT_MLS{p,corner}(1:16) = CONNECT_TEMP_N(1:16);
 
% %   CONNECT_TEMP(1)             = spElems_corner + floor(spElems_corner/(NN(1)-1));
%   CONNECT_TEMP(1)             = spElems_corner + ceil(spElems_corner/(NN(1)-1)) - 1;
%  npoints{CONNECT_TEMP(1)}    = [npoints{CONNECT_TEMP(1)} [p;corner]];
%  CONNECT_TEMP(2)             = CONNECT_TEMP(1) + 1;
%  npoints{CONNECT_TEMP(2)}    = [npoints{CONNECT_TEMP(2)} [p;corner]];
%  CONNECT_TEMP(3)             = CONNECT_TEMP(2) + NN(1);
%  npoints{CONNECT_TEMP(3)}    = [npoints{CONNECT_TEMP(3)} [p;corner]];
%  CONNECT_TEMP(4)             = CONNECT_TEMP(1) + NN(1);
%  npoints{CONNECT_TEMP(4)}    = [npoints{CONNECT_TEMP(4)} [p;corner]];
%  
%  CONNECT_MLS{p,corner}(1:4) = CONNECT_TEMP(1:4);