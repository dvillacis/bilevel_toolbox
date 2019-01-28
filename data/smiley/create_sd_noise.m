clear all;
close all;
clc;

smiley = im2double(imread('data/smiley/smiley.png'));
smiley = 1-imbinarize(smiley,1e-4);
[m,n] = size(smiley);
black_pixels = smiley == 1;
% Adding gaussian noise std 0.1 to black pixels
noise1 = 0.2*rand(m,0.5*n);
noise2 = 0.8*rand(m,0.5*n);
noise = horzcat(noise1,noise2);
noise(black_pixels) = 0;
smiley_noisy = smiley + noise;
    
%show image
figure(1)
imagesc_gray(smiley);
figure(2)
imagesc_gray(smiley_noisy);

% write image
imwrite(smiley,'data/smiley/1_smiley_original.png');
imwrite(smiley_noisy,'data/smiley/1_smiley_noisy.png');