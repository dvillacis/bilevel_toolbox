clear all;
clc;

% Setup the size of the dataset
N = 1;

% Check if reduced playing dataset if available
dataset_name = 'playing_cards_reduced.mat';
if exist(dataset_name,'file') == 0
   % There is no file, so I download it from url
   url_path = 'https://storage.googleapis.com/pagina-personal.appspot.com/img_research/img_datasets/playing_cards_dataset_reduced.mat';
   websave(dataset_name,url_path);
end

% Load the mat dataset
load playing_cards_reduced.mat;

% Generate the patches for random rectangles of 128x128 pixels
[xmax,ymax] = size(im_clean);
for i=1:N
    rnd_x0 = randi(xmax-128);
    rnd_y0 = randi(ymax-128);
    clean_i = imcrop(im_clean,map,[rnd_x0 rnd_y0 rnd_x0+128 rnd_y0+128]);
    imagesc_gray(clean_i);
end
