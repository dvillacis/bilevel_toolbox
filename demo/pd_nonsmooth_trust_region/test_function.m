clear all;
clc;

dataset = DatasetInFolder('data/smiley','*_smiley_original.png','*_smiley_noisy.png');
noisy = dataset.get_corrupt(1);

Z = zeros(10,15);

for X = 1:5
    for Y = 10:15
        u = solve_pd_ntr_lower_level([X,Y],noisy);
        Z(X,Y) = eval_pd_ntr_upper_level(u,[X,Y],dataset);
        fprintf('Evaluated %d,%d\n',X,Y);
    end
end

[A,B] = meshgrid(1:5,10:15);