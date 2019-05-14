clear all;
close all;
clc;

dataset = DatasetInFolder('data/circle_dataset','*_circle_original.png','*_circle_noisy.png');

c = [];

r = 3.04:0.0001:3.07;

for i = r
    u = solve_ntr_lower_level(i,dataset.get_corrupt(1));
    cost = eval_ntr_upper_level(u,i,dataset);
    c = [c cost];
    fprintf("Evaluated sol=%d, cost=%f\n",i,cost);
end

plot(r,c);