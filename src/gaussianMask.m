% Name: gaussianMask.m
% Created 20/12/2011
%
% Description: Masks a matrix with a gaussian aperture
%               
% @author Kevin O'Holleran <ko311@cam.ac.uk>
%
% @param data mxn matrix
% @param sigma double - determines sigma for each dimension in pixels.
% @param invert - inverts mask by subtracting from array of ones
% @param quadswap - boolean that switches fftshift on or off before mask is
% applied.
function [ maskedData ] = gaussianMask(data,sigma,invert,quadswap, ind)
    if nargin < 3
        quadswap = false;
        invert = false;
        ind = [0,0];
    end
    if nargin < 4
        quadswap = false;
        ind = [0,0];
    end
    dims = size(data);
    N = dims(2);
    M = dims(1);
    if  mod(N,2)==0
        npdfx = exp(-(((-N/2-ind(1)):((N/2-1-ind(1)))).^2/(2*sigma^2)));      
    else
        npdfx = exp(-((-(N-1)/2-ind(1)):((N-1)/2-ind(1))).^2/(2*sigma^2));
    end 
    if  mod(M,2)==0       
        npdfy = exp(-(((-M/2-ind(2)):((M/2-1-ind(2)))).^2/(2*sigma^2)));
    else
        npdfy = exp(-((-(M-1)/2-ind(2)):((M-1)/2-ind(2))).^2/(2*sigma^2));
    end 
    a = repmat(npdfx,M,1);
    b =  repmat(transpose(npdfy),1,N);
    mask = a.*b;
    if quadswap
        mask = fftshift(mask);
    end
    if invert
        mask = -(mask-1);
    end
    maskedData = mask.*data;
end

