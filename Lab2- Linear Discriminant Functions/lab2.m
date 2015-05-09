function A=lab2(XA,XB,LR)
%XA = samples from W1
%XB = samples from W2
%LR = Learning Rate
A=[0,0,1];
countw1=size(XA,1);
countw2=size(XB,1);
totalcount=countw1+countw2;
Y=repmat(1,countw1+countw2,3);
Y(1:countw1,2:3)=XA(:,2:3);
Y(countw1+1:totalcount,2:3)=-1*XB(:,2:3);
Y(countw1+1:totalcount,1)=-1*Y(countw1+1:totalcount,1);
Y=Y.';
J=0;

for iteration = 1:300
AY=A*Y;

for i = 1:totalcount
    if AY(1,i) <= 0
        J=J-Y(:,i);
    end
end

if J == 0
    break
else
    A=A-LR*J.';
    J=0;
end
end
iteration
figure;

for i = 1:countw1
    plot(XA(i,2),XA(i,3),'ob')
    hold on;
end
for i = 1:countw2
    plot(XB(i,2),XB(i,3),'xr')
    hold on;
end
x=1:5;
plot(x,(-1*A(2)/A(3))*x+A(1)/A(3))

xlabel('sepal width')
ylabel('petal length')
hold off;
end