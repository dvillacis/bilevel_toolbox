clear all;
close all;
clc;

circle = imread('data/circle_dataset/circle.png');
[m,n] = size(circle);

% Adding gaussian noise std 0.1 to black pixels
%circle_noisy = imnoise(circle,'salt & pepper',0.2);
circle_noisy = imnoise(circle,'gaussian');
% circle_noisy = circle;
% circle_noisy(:,1:n/2) = imnoise(circle(:,1:n/2),'gaussian',0.1);
% circle_noisy(:,n/2:n) = imnoise(circle(:,n/2:n),'salt & pepper',0.2);
    
%show image
figure(1)
imagesc_gray(circle);
figure(2)
imagesc_gray(circle_noisy);

% write image
imwrite(im2double(circle),'data/circle_dataset/1_circle_original.png');
imwrite(im2double(circle_noisy),'data/circle_dataset/1_circle_noisy.png');