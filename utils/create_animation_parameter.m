function [] = create_animation_parameter(input_tensor,pathname,op)
%CREATE_ANIMATION Create a GIF animation out of a tensor
%   Detailed explanation goes here
h = figure;
axis tight manual
filename = pathname;
sol_1 = op.val(input_tensor(:,:,1));
[~,~,NumSol] = size(input_tensor);
[Msol,Nsol] = size(sol_1);
[a,b] = meshgrid(1:Msol,1:Nsol);
for i = 1:NumSol
    sol_i = op.val(input_tensor(:,:,i));
    surf(a,b,sol_i);
    drawnow;
    % Capture the plot as an image 
    frame = getframe(h); 
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256); 
    % Write to the GIF File 
    if i == 1 
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
    else 
        imwrite(imind,cm,filename,'gif','WriteMode','append'); 
    end 
end
end

