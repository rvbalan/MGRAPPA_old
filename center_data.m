function A = center_data(A,N,channel_len)
size(A)
% switch to fft domain 
s = zeros(N,N,channel_len);
size(s)
for icha = 1:channel_len % parfor
    s(:,:,icha)=fftshift(fft2(A(:,:,icha)));
end

size(s)

s(:,:,1) = fftshift(s(:,:,1));
[m,x,y] = max2D(abs(s(:,:,1)));
s(:,:,1) = ifftshift(s(:,:,1));

size(s)
%center
for icha = 1:channel_len
    s(:,:,icha) = fftshift(s(:,:,icha));
    s(:,:,icha)=s([x:N 1:x-1],[y:N 1:y-1],icha);
    s(:,:,icha) = ifftshift(s(:,:,icha));
end
size(s)
parfor icha = 1:channel_len
    A(:,:,icha)=ifft2(ifftshift(s(:,:,icha)));
end
size(A)