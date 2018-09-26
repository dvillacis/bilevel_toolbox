function y = norm_tv(I)
%NORM_TV 2 Dimentional TV norm

[dx, dy] = gradient_op(I);
temp = sqrt(abs(dx).^2 + abs(dy).^2);

y = sum(temp(:));
%y = reshape(sum(sum(temp,1),2),[],1);

end
