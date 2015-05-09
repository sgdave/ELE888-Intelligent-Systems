%% B
clc
clear all
runlab1

%% a
close all
x = [3.3, 4.4, 5.0, 5.7, 6.3];
y = [2 6;4.4 3;5 3.5;5.3 2;5.5 2.5;6.6 3.5;4.5 6.1]
[posteriors_x,g_x] = lab1(x,trainingSet,y)
