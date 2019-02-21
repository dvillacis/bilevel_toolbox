function [] = create_animation_parameter(input_tensor,pathname)
%CREATE_ANIMATION Create a GIF animation out of a tensor
%   Detailed explanation goes here
h = figure;
axis tight manual
filename = pathname;
[Msol,Nsol,NumSol] = size(input_tensor);
[a,b] = meshgrid(1:Msol,1:Nsol);
for i = 1:NumSol
    surf(a,b,input_tensor(:,:,i));
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

