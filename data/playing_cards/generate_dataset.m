clear all;
clc;

% Setup the size of the dataset
N = 10;

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
written = 0;
while written <= N
    rnd_x0 = randi([1500,xmax-1500],1,1);
    rnd_y0 = randi([1500,ymax-1500],1,1);
    clean_i = imcrop(im_clean,[rnd_x0 rnd_y0 511 511]);
    noise_i = imcrop(im_noise3,[rnd_x0 rnd_y0 511 511]);
    if ~isempty(clean_i)
        write_img(clean_i,sprintf('data/playing_cards/%d_playing_cards_original.tif',written));
        write_img(noise_i,sprintf('data/playing_cards/%d_playing_cards_noisy.tif',written));
        fprintf('Image %d wrote successfully\n',written);
        written = written + 1;
    end
end

function [] = write_img(img, img_file)
    img = uint16(img); % Convert data type
    BitsPerSample = 16;
    SampleFormat = 1;

    % Create a Tiff object
    t = Tiff(img_file,'w');

    tagstruct.ImageLength = size(img,1);
    tagstruct.ImageWidth = size(img,2);
    tagstruct.SamplesPerPixel = size(img,3);
    tagstruct.SampleFormat = SampleFormat;
    tagstruct.BitsPerSample = BitsPerSample;

    tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
    tagstruct.Compression = 1;
    tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
    tagstruct.Software = 'MATLAB';

    t.setTag(tagstruct);
    t.write(img);
    t.close()
end
