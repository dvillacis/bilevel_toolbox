% Calculate the pointwise euclidean norm
function ngrepl=pointwise_norm_replicated(g, n, m)
    gg=reshape(g, [n, m, 2]);
    ng=sqrt(sum(gg.^2, 3));
    ngrepl=reshape(cat(3, ng, ng), size(g));
end