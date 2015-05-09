function testlab2(A,XA,XB)
countw1=size(XA,1);
countw2=size(XB,1);
W0=A(1);
W1=A(2);
W2=A(3);
correct=0;
for i = 1:countw1
   G=W0+W1*XA(i,2)+W2*XA(i,3);
   if G > 0
      correct=correct+1; 
   end
end
for i = 1:countw2
   G=W0+W1*XB(i,2)+W2*XB(i,3);
   if G < 0
      correct=correct+1; 
   end
end
percentage=correct*100/(countw1+countw2);
disp(percentage)
end