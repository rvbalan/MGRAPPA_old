function [ChannelImage, channel_len] = getChannels(Image,theta_brain,theta_profiles,...
    x_offset,y_offset,profile_smooth_level)

    if theta_brain ~= 0
        Image = imrotate(Image,theta_brain,'bilinear','crop');
    end
    

    load('20210313_MPRAGE\ChannelImage_sos.mat','ChannelImage');
    profiles = ChannelImage;

    N = size(profiles,1);
    channel_len = size(profiles,3);
        
    %Center Data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %get signal data
    s = zeros(N,N,channel_len);
    parfor icha = 1:channel_len
        s(:,:,icha)=fftshift(fft2(profiles(:,:,icha)));
    end

    s(:,:,1) = fftshift(s(:,:,1));
    [m,x,y] = max2D(abs(s(:,:,1)));
    s(:,:,1) = ifftshift(s(:,:,1));

    %center
    for icha = 1:channel_len
        s(:,:,icha) = fftshift(s(:,:,icha));
        s(:,:,icha)=s([x:N 1:x-1],[y:N 1:y-1],icha);
        s(:,:,icha) = ifftshift(s(:,:,icha));
    end

    parfor icha = 1:channel_len
        profiles(:,:,icha)=ifft2(ifftshift(s(:,:,icha)));
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %DONE Center Data

    % Up sample
    profiles2 = zeros(2*N,2*N,channel_len);
    parfor icha = 1:channel_len
        profiles2(:,:,icha) = ...
            imresize(squeeze(profiles(:,:,icha)),2,'bilinear');
    end
    profiles = profiles2;
    N = size(profiles,1);
    
    % Crop
    profiles2 = zeros(N/2,N/2,channel_len);
    parfor icha = 1:channel_len
        profiles2(:,:,icha) = profiles((N/4):(3*N/4)-1,(N/4):(3*N/4)-1,icha);
    end
    profiles = profiles2;
    N = size(profiles,1);
    
    % Smooth profiles 
    parfor icha = 1:channel_len
        s = fftshift(fft2(squeeze(profiles(:,:,icha))));
        s(1:floor(N*profile_smooth_level/2),:)=0;
        s(ceil(N-N*profile_smooth_level/2):N,:)=0;
        s(:,1:floor(N*profile_smooth_level/2))=0;
        s(:,ceil(N-N*profile_smooth_level/2):N)=0;
        profiles(:,:,icha) = ifft2(ifftshift(s));
    end
    
    profiles_rot = zeros(N,N,channel_len);
    if theta_profiles ~= 0
        parfor icha = 1:channel_len
            profiles_rot(:,:,icha) = imrotate(squeeze(profiles(:,:,icha)),...
                                        theta_profiles,'bilinear','crop');
        end
    else
        profiles_rot = profiles;
    end

    profiles_rot_x_y = translate_x_y(profiles_rot, x_offset, y_offset);

    ChannelImage = zeros(N,N,channel_len);
    parfor icha = 1:channel_len % parfor
        ChannelImage(:,:,icha) = Image.*(profiles_rot_x_y(:,:,icha));
    end
