function outer = outer_product(p,q)
  a = p(:,1).*q(:,1);
  b = p(:,1).*q(:,2);
  c = p(:,2).*q(:,1);
  d = p(:,2).*q(:,2);
  outer = [a,b;c,d];
end
