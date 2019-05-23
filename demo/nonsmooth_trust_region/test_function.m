clear all;
close all;
clc;

dataset = DatasetInFolder('data/circle_dataset','*_circle_original.png','*_circle_noisy.png');

c = [];

r = 0.5:0.5:90;

for i = r
    u = solve_ntr_lower_level(i,dataset.get_corrupt(1));
    cost = eval_ntr_upper_level(u,i,dataset);
    c = [c cost];
    fprintf("Evaluated sol=%d, cost=%f\n",i,cost);
end

plot(r,c);
matlab2tikz('/Users/dvillacis/OneDrive - Escuela Politécnica Nacional/DOCTORADO/ARTICLES/Bilevel_SD_ROF/plots/circle_1_cost.tex');