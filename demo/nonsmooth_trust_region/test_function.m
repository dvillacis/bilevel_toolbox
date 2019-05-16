clear all;
close all;
clc;

dataset = DatasetInFolder('data/circle_dataset','*_circle_original.png','*_circle_noisy.png');

c = [];

r = 2.91:1e-3:3.06;

for i = r
    u = solve_ntr_lower_level(i,dataset.get_corrupt(1));
    cost = eval_ntr_upper_level(u,i,dataset);
    gradient_parameters.regularized_model = false;
    grad = solve_ntr_gradient(u,i,dataset,gradient_parameters);
    c = [c cost];
    fprintf("Evaluated sol=%d, cost=%f, grad=%f\n",i,cost,grad);
end

plot(r,c);