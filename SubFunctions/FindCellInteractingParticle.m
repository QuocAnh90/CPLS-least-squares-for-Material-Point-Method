function [CONNECT_MLS_G,cpoints]=FindCellInteractingParticle...
    (sp,spElems_corner,cpoints,NN,corner,CONNECT_MLS_G)

CONNECT_TEMP_G = zeros(1,9);

 % Corner1-cell interaction
 CONNECT_TEMP_G(1) = spElems_corner-(NN(1)-1)-1;     
 cpoints{CONNECT_TEMP_G(1)} = [cpoints{CONNECT_TEMP_G(1)} [sp;corner]];
 CONNECT_TEMP_G(2) = spElems_corner-(NN(1)-1);       
 cpoints{CONNECT_TEMP_G(2)} = [cpoints{CONNECT_TEMP_G(2)} [sp;corner]];
 CONNECT_TEMP_G(3) = spElems_corner-(NN(1)-1)+1;     
 cpoints{CONNECT_TEMP_G(3)} = [cpoints{CONNECT_TEMP_G(3)} [sp;corner]];
 CONNECT_TEMP_G(4) = spElems_corner-1;               
 cpoints{CONNECT_TEMP_G(4)} = [cpoints{CONNECT_TEMP_G(4)} [sp;corner]];
 CONNECT_TEMP_G(5) = spElems_corner;                 
 cpoints{CONNECT_TEMP_G(5)} = [cpoints{CONNECT_TEMP_G(5)} [sp;corner]];
 CONNECT_TEMP_G(6) = spElems_corner+1;               
 cpoints{CONNECT_TEMP_G(6)} = [cpoints{CONNECT_TEMP_G(6)} [sp;corner]];
 CONNECT_TEMP_G(7) = spElems_corner+(NN(1)-1)-1;     
 cpoints{CONNECT_TEMP_G(7)} = [cpoints{CONNECT_TEMP_G(7)} [sp;corner]];
 CONNECT_TEMP_G(8) = spElems_corner+(NN(1)-1);       
 cpoints{CONNECT_TEMP_G(8)} = [cpoints{CONNECT_TEMP_G(8)} [sp;corner]];
 CONNECT_TEMP_G(9) = spElems_corner+(NN(1)-1)+1;     
 cpoints{CONNECT_TEMP_G(9)} = [cpoints{CONNECT_TEMP_G(9)} [sp;corner]];
 
 CONNECT_MLS_G{sp,corner}(1:9) = CONNECT_TEMP_G(1:9);