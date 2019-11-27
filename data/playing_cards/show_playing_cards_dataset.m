clear all
clc

%% Load dataset
dataset = DatasetInFolder('data/playing_cards','*_playing_cards_original.tif','*_playing_cards_noisy.tif');

%% Show Dataset
imagesc_gray(dataset.get_target(1),1,'Original 1','421');
imagesc_gray(dataset.get_corrupt(1),1,'Noisy 1','422');
imagesc_gray(dataset.get_target(2),1,'Original 2','423');
imagesc_gray(dataset.get_corrupt(2),1,'Noisy 2','424');
imagesc_gray(dataset.get_target(3),1,'Original 3','425');
imagesc_gray(dataset.get_corrupt(3),1,'Noisy 3','426');
imagesc_gray(dataset.get_target(6),1,'Original 4','427');
imagesc_gray(dataset.get_corrupt(6),1,'Noisy 4','428');