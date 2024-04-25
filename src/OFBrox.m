function [D of] = OFBrox(A, sigma, alpha, gamma)
addpath .\Brox\
if ~exist('sigma')
    sigma = 0.5;
end
if ~exist('alpha')
    alpha = 80;
end
if ~exist('gamma')
    gamma =5;
end
D = cell(size(A,3)-1,1);

for t = 2:size(A,3)
        t
        flow = mex_OF(double(repmat(A(:,:,t-1),[1 1 3])), double(repmat(A(:,:,t),[1 1 3])),sigma, alpha, gamma);
        
        im1 = A(:,:,t-1);
        imw2= imwarp( A(:,:,t), flow);
        figure(2);
        subplot(131)
        imagesc(im1-imw2); axis equal tight ; title(strcat('warp2:',num2str((mean2((im1-imw2).^2)))))
        dd = (im1-imw2).^2;
        quantile(dd(:), [0.5 0.75 0.9])
        c = caxis;
        subplot(132)
        %imagesc(imwarp( A(:,:,t-1), flow)-A(:,:,t)); axis equal tight ;  title(strcat('warp1:',num2str(mean2((imwarp( A(:,:,t-1), flow)-A(:,:,t)).^2))));
        
        g1 = gradient(im1);
        g2 = gradient(imw2);
        imagesc(sqrt(sum((g1-g2).^2,3))); axis equal tight ;  title(strcat('grad:',num2str(mean2(sqrt(sum((g1-g2).^2,3))))));
        caxis(c)
        subplot(133)
        imagesc(im1-A(:,:,t)); axis equal tight ; title(strcat('orig:',num2str(mean2((im1-A(:,:,t)).^2))))
        caxis(c)
        
        [x,y] = meshgrid(1:5:size(im1,2),1:5:size(im1,1));
        f1  = flow(:,:,1);
        f2  = flow(:,:,2);
        %quiver(x(:),y(:),f1(sub2ind(size(im1(:,:,1)), y(:), x(:))),f2(sub2ind(size(im1(:,:,1)), y(:), x(:))));
         
        h= figure('units','normalized','outerposition',[0 0 1 1]);
        clf
        imshowpair(A(:,:,t-1),A(:,:,t)); hold on  

        [x,y] = meshgrid(1:5:size(A,2),1:5:size(A,1));
        quiver(x(:),y(:),f1(sub2ind(size(im1(:,:,1)), y(:), x(:))),f2(sub2ind(size(im1(:,:,1)), y(:), x(:))),'Color', 'm');
        D{t-1} = cat(3, f1, f2);
        pause(0.01)

end

%%
mag = cellfun(@(x)(sqrt(x(:,:,1).^2+x(:,:,2).^2)), D,'UniformOutput',false);
for k = 1:length(mag)
    of(:,:,k) = mag{k};
end