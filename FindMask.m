function image_mask = FindMask(Image,pmask,threshold)
% FindMask Returns a binary mask of same size the input image
%   Use:
% image_mask = FindMask(Image,pmask,threshold)
%
% Radu Balan, 8 March 2024

if nargin < 2
    pmask = 3; % default number of voxels around a pixel detected in the brain region
    threshold = 0.1; % default threshold to define a brain pixel
elseif nargin < 3
    threshold = 0.1; % default threshold to define a brain pixel
end

% Compute the binary mask
image_mask = zeros(size(Image));
lmin = min(min(abs(Image)));
lmax = max(max(abs(Image)));
mask_threshold = lmin + threshold*(lmax-lmin); % 10% above min level
nx = size(Image,1);
ny = size(Image,2);
for ix =1:nx
    for iy = 1:ny
        if (abs(Image(ix,iy)) >= mask_threshold)
            % Fill out with 1 a small square of size pmask x pmask centered
            % around (ix,iy)
            for jx = ix-pmask:ix+pmask
                for jy = iy-pmask:iy+pmask
                    if (1 <= jx) & (jx<= nx) & (1 <= jy) & (jy <= ny)
                        image_mask(jx,jy) = 1;
                    end
                end
            end
        end
    end
end



end
