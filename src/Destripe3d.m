%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Destripe a 3d stack (plane by plane),
% lam94  23.04.2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function im_filtered = Destripe3d(A)

[w h] = size(A);
fA = fft(fft(A,[],1),[],2);

% figure; imagesc(log(1+(fftshift(sum(abs(fA),3)))))
% spotIndex = floor(ginput);

LL = log(1+(fftshift(sum(abs(fA),3))));
[pks,locs] = findpeaks(LL(:,floor(size(LL,2)/2+1)), 'SortStr' , 'descend');
N=4;
spotIndex = floor([ones(N,1)*floor(size(LL,2)/2+1) locs(2:5) ]);

mask = ones(size(A(:,:,1)));
for k = 1:size(spotIndex,1)    
    mask = min(mask, gaussianMask(mask,4,true,true, spotIndex(k,:) - w/2 ));
    mask = min(mask,gaussianMask(mask,4,true,true, w/2 - spotIndex(k,:)));
end

mask = repmat(mask, [1,1, size(fA,3)]);
im_filtered = real(ifft(ifft(fA.*mask,[], 2), [], 1));
