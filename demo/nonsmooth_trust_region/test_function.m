clear all;
close all;
clc;

dataset = DatasetInFolder('data/circle_dataset','*_circle_original.png','*_circle_noisy.png');

c = [];

for i = 1:0.1:5
    fprintf("Evaluating %d\n",i);
    u = solve_ntr_lower_level(i,dataset.get_corrupt(1));
    cost = eval_ntr_upper_level(u,i,dataset);
    c = [c cost];
end