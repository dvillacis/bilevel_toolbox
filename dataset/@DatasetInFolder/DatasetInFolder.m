classdef DatasetInFolder < Dataset
    properties
        FolderPath
        OriginalRegex
        CorruptRegex
    end
    methods

        % Constructor
        function dataset = DatasetInFolder(FolderPath, OriginalRegex, CorruptRegex)
            dataset.FolderPath = FolderPath;
            dataset.OriginalRegex = OriginalRegex;
            dataset.CorruptRegex = CorruptRegex;

            % Get the number of training images
            training_pairs = dir(strcat(dataset.FolderPath,'/',dataset.OriginalRegex));
            dataset.NumPairs = length(training_pairs);
        end
        function [target, corrupt] = get_pair(dataset, id)
            % GET_PAIR Get image pair from dataset folder
            %   [TARGET, CORRUPT] = GET_PAIR(ID) gets a training sample tuple with id = ID
            %fprintf('Getting a pair of images for id:%d\n', id);
            target_path = strrep(dataset.OriginalRegex,'*',num2str(id));
            corrupt_path = strrep(dataset.CorruptRegex,'*',num2str(id));
            target = imadjust(im2double(imread(strcat(dataset.FolderPath,'/',target_path))));
            corrupt = imadjust(im2double(imread(strcat(dataset.FolderPath,'/',corrupt_path))));
        end
        function target = get_target(dataset, id)
            %fprintf('Getting target image for pair id:%d\n', id);
            target_path = strrep(dataset.OriginalRegex,'*',num2str(id));
            target = imadjust(im2double(imread(strcat(dataset.FolderPath,'/',target_path))));
        end
        function corrupt = get_corrupt(dataset, id)
            %fprintf('Getting corrupt image for pair id:%d\n', id);
            corrupt_path = strrep(dataset.CorruptRegex,'*',num2str(id));
            corrupt = imadjust(im2double(imread(strcat(dataset.FolderPath,'/',corrupt_path))));
        end

    end
end
