function [funh,icase] = FUNh(rx,rz,w1,w2,s1,s2)

 if rx<=s1 && rz<=w1
     icase=7;
     funh = sqrt((s1-rx)^2+(w1-rz)^2);
 elseif rz<=w1 && rx>=s1 && rx<=s2
     icase=4;
     funh = w1-rz;
     elseif rx>=s2 && rz<=w1
         icase=1;
         funh = sqrt((rx-s2)^2+(w1-rz)^2);
         elseif rx>=s2 && rz>=w1 && rz<=w2
             icase=2;
             funh = rx-s2;
             elseif rx>=s2 && rz>=w2
                 icase=3;
                 funh = sqrt((rx-s2)^2+(rz-w2)^2);
                 elseif rz>=w2 && rx>=s1 && rx<=s2
                     icase=6;
                     funh = rz-w2;
                     elseif rz>=w2 && rx<=s1
                         icase=9;
                         funh=sqrt((s1-rx)^2+(rz-w2)^2);
                         elseif rx<=s1 && rz>=w1 && rz<=w2
                             icase=8;
                             funh = s1-rx;
                             elseif rx>=s1 && rx<=s2 && rz>=w1 && rz<=w2
                                 icase=5;
                                 funh = 0.0;
 end
                       
end