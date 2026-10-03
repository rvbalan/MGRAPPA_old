function res = ifft2c(x)

res = zeros(size(x));
parfor n=1:size(x,3)
	res(:,:,n) = ifft2(ifftshift(x(:,:,n)));
end
