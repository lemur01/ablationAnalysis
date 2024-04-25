%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lam94@cam.ac.uk, 2013
% Read tif stack with name filename
% Input: filename
% Output:   A   - 3d matrix with image data
%           info - information on the image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [A]=  Read3d(filename, starti, endi)

info = imfinfo(filename,'tif');
    num_images = numel(info);   
    if ~exist('starti')
        starti = 1;
    end
    if ~exist('endi')
        endi = num_images;
    end
    num_images_read=endi-starti+1;
    A = zeros(info(1).Height(1),info(1).Width ,num_images_read);
    for k = starti:endi
        if rem(k,100) == 0
            disp(k)
        end
        A(:,:,k) = double(imread(filename, k));%, 'Info', info));
%         k
    end
end