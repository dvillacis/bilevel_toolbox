function [cost] = eval_pd_ntr_upper_level(u,lambda,dataset)
    original = dataset.get_target(1);
    cost = 0.5*norm(u(:)-original(:)).^2 + 0.5*0.1*norm(lambda(:)).^2;
end

