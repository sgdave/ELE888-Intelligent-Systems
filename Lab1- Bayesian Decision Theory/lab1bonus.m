function [posteriors_x,g_x]=lab1(x,Training_Data)
%% 

D=Training_Data;

% D is MxN (M samples, N columns = N-1 features + 1 label)
[M,N]=size(D);    

feature_column = 2; %% use 1 when wanting to experiment with Sepal Length i.e question 6
f=D(:,feature_column);  % feature samples
la=D(:,N); % class labels
%% 
