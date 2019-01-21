% Compute PSNR between iages im1 and im2, assumed to have dynamir range dr
function psnrval = psnr(im1, im2, dr)
%Compute MSE
sqr = (double(im1) - double(im2)).^2;
mse = 1/numel(im1) * sum(sqr(:));
%Compute PSNR
if mse == 0
    psnrval = Inf;
else
    psnrval = 10*log10(dr^2/mse);
end
end
