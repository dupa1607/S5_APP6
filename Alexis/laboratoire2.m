%% Question 11

clear all
close all
clc

e = [0.0843 0.02603 0.01048 0.00319 0.00040];
delta_t = [0.05 0.04 0.03 0.02 0.01];

p = [];
N = size(e);
N = N(2);
for i=2:1:N
    p_temp = log(e(i-1)/e(i))/log(delta_t(i-1)/delta_t(i));
    p = [p p_temp];
end

p = [];
for i=1:1:N-1
    p_temp = log(e(i)/e(5))/log(delta_t(i)/delta_t(5));
    p = [p p_temp];
end

%% Question 12