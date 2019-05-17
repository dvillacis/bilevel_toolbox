clear all;
close all;
clc;

smiley = imread('data/smiley/1_smiley_original.png');
[m,n] = size(smiley);
black_pixels = smiley == 0;
% Adding gaussian noise std 0.1 to black pixels
smiley_noisy = smiley;
smiley_noisy(black_pixels) = imnoise(smiley(black_pixels),'gaussian',0,0.05);
smiley_noisy(black_pixels(:,1:n/2)) = imnoise(smiley(black_pixels(:,1:n/2)),'gaussian',0,0.5);

% smiley_noisy(:,1:n/2) = imnoise(smiley(:,1:n/2),'gaussian',0,0.1);
% smiley_noisy(:,n/2:n) = imnoise(smiley(:,n/2:n),'gaussian',0,0.6);

    
%show image
figure(1)
imagesc_gray(smiley);
figure(2)
imagesc_gray(smiley_noisy);

% write image
imwrite(im2double(smiley),'data/smiley/1_smiley_original.png');
imwrite(im2double(smiley_noisy),'data/smiley/1_smiley_noisy.png');