%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ELE 888/ EE 8209: LAB 1: Bayesian Decision Theory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [posteriors_x,g_x]=lab1(x,Training_Data,y)
%%
% x = individual sample to be tested (to identify its probable class label)
% featureOfInterest = index of relevant feature (column) in Training_Data 
% Train_Data = Matrix containing the training samples and numeric class labels
% posterior_x  = Posterior probabilities
% g_x = value of the discriminant function

D=Training_Data;

% D is MxN (M samples, N columns = N-1 features + 1 label)
[M,N]=size(D);    

feature_column = 2; %% use 1 when wanting to experiment with Sepal Length i.e question 6
f=D(:,feature_column);  % feature samples
la=D(:,N); % class labels


%% %%%%Prior Probabilities%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Hint: use the commands "find" and "length"

disp('Prior probabilities:');

Pr1 = sum(la==1)/length(la)
Pr2 = sum(la==2)/length(la)

%% %%%%%Class-conditional probabilities%%%%%%%%%%%%%%%%%%%%%%%

disp('Mean & Std for class 1 & 2');
m11 = mean(f(la==1)) % mean of the class conditional density p(x/w1)
std11 = std(f(la==1)) % Standard deviation of the class conditional density p(x/w1)

m12 = mean(f(la==2)) % mean of the class conditional density p(x/w2)
std12= std(f(la==2)) % Standard deviation of the class conditional density p(x/w2)

feature = x(feature_column);

disp(['Conditional probabilities for x=' num2str(x)]);
cp11 = 1/(sqrt(2*pi)*std11)*exp(-.5*((feature-m11)/std11)^2) % use the above mean, std and the test feature to calculate p(x/w1)

cp12 = 1/(sqrt(2*pi)*std12)*exp(-.5*((feature-m12)/std12)^2) % use the above mean, std and the test feature to calculate p(x/w2)

%% %%%%%%Compute the posterior probabilities%%%%%%%%%%%%%%%%%%%%

disp('Posterior prob. for the test feature');

Px = cp11*Pr1 + cp12*Pr2;

pos11 = cp11*Pr1/Px% p(w1/x) for the given test feature value

pos12= cp12*Pr2/Px% p(w2/x) for the given test feature value

posteriors_x = [pos11, pos12];

plot(f*pos11,f*Px);
figure;
plot(f*pos12,f*Px);
figure;

%% %%%%%%Discriminant function for min error rate classifier%%%

disp('Discriminant function for the test feature');

[C,I]= min(posteriors_x);

g_x = I  % compute the g(x) for min err rate classifier.

%% Finding Threshhold

% pos11 = 'cp11*Pr1/Px'
% pos12 = 'cp12*Pr2/Px'
% pos11 == pos12 @x=Th1

% =>
% cp11*Pr1/Px = cp12*Pr2/Px @x=Th1
% => 
% 1/(sqrt(2*pi)*std11)*exp(-.5*((feature-m11)/std11)^2) * Pr1
% ==
% 1/(sqrt(2*pi)*std12)*exp(-.5*((feature-m12)/std12)^2) * Pr2
syms th
eq = 1/(sqrt(2*pi)*std11)*exp(-.5*((th-m11)/std11)^2) * Pr1 == 1/(sqrt(2*pi)*std12)*exp(-.5*((th-m12)/std12)^2) * Pr2;
th = solve(eq, th);
th = double(th);
Th1 = th(th>0)

plot(th)
%% Unbalanced Loss
Theta = (0.8/0.9)*(Pr2/Pr1);

syms th
eq = 1/(sqrt(2*pi)*std11)*exp(-.5*((th-m11)/std11)^2) * Pr1 == 1/(sqrt(2*pi)*std12)*exp(-.5*((th-m12)/std12)^2) * Pr2;
th = solve(eq, th);
th = double(th);
Th1 = th(th>Theta)
%% Bonus

f=D(:,1:2);  % feature samples
la=D(:,N); % class labels
[dee,col]=size(la);
sigma=cov(f);

%% 

disp('Mean & Std for class 1 & 2');
mu1 = mean(f(la==1)) % mean of the class conditional density p(x/w1)
mu2 = mean(f(la==2)) % mean of the class conditional density p(x/w2)

lsigmal = det(sigma)
isigma=inv(sigma);

feature = y(1);
feature2 = y(2);

disp(['Conditional probabilities for x=' num2str(x)]);
cp11 = 1/((2*pi)^((length(sigma))/2)*(lsigmal^0.5))*exp((-1/2)*(transpose(feature-mu1))*isigma*(feature-mu1));% use the above mean, std and the test feature to calculate p(x/w1)
cp21 = 1/((2*pi)^((length(sigma))/2)*(lsigmal^0.5))*exp((-1/2)*(transpose(feature2-mu2))*isigma*(feature2-mu2));
cp12 = 1/((2*pi)^((length(sigma))/2)*(lsigmal^0.5))*exp((-1/2)*(transpose(feature-mu2))*isigma*(feature-mu2)); % use the above mean, std and the test feature to calculate p(x/w2)
cp22 = 1/((2*pi)^((length(sigma))/2)*(lsigmal^0.5))*exp((-1/2)*(transpose(feature2-mu1))*isigma*(feature2-mu1));

%% 

disp('Posterior prob. for the test feature');

Px = cp11*Pr1 + cp12*Pr2;
Px2 = cp21*Pr1 + cp22*Pr2;
pos11 = cp11*Pr1/Px% p(w1/x) for the given test feature value
pos21 = cp21*Pr1/Px2
pos12= cp12*Pr2/Px% p(w2/x) for the given test feature value
pos22 = cp22*Pr2/Px2
posteriors_x = [pos11, pos12, pos21, pos22];

plot3(pos11,pos12,pos21);hold on 
figure;
%plot3(pos21,pos22);

%% 

disp('Discriminant function for the test feature');

[C,I]= min(posteriors_x);

g_x = I  % compute the g(x) for min err rate classifier.