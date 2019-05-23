clear all;
close all;
clc;

smiley = imread('data/smiley/1_smiley_original.png');
[m,n] = size(smiley);
black_pixels = find(smiley == 0);
white_pixels = find(smiley ~= 0);
% Adding gaussian noise std 0.1 to black pixels
smiley_noisy = smiley;
smiley_noisy(black_pixels) = imnoise(smiley_noisy(black_pixels),'gaussian',0,0.5);
smiley_noisy(white_pixels) = imnoise(smiley_noisy(white_pixels),'gaussian',0,0.01);


% smiley_noisy(:,1:n/2) = imnoise(smiley(:,1:n/2),'gaussian',0,0.1);
% smiley_noisy(:,n/2:n) = imnoise(smiley(:,n/2:n),'gaussian',0,0.6);

    
%show image
figure(1)
imagesc_gray(smiley,1,'original','131');
imagesc_gray(smiley_noisy,1,'noisy','132');
imagesc_gray(smiley-smiley_noisy,1,'diff','133');

% write image
imwrite(im2double(smiley),'data/smiley/1_smiley_original.png');
imwrite(im2double(smiley_noisy),'data/smiley/1_smiley_noisy.png');