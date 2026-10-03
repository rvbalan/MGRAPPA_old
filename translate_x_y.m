function P = translate_x_y(P, x_offset, y_offset)

if length(size(P))==3
    Nx = size(P,1);
    if x_offset > 0
        P(1:x_offset,:,:)=0;
        P = P([1+x_offset:Nx, 1:x_offset],:,:);
    elseif x_offset < 0
        P(Nx+x_offset+1:Nx,:,:)=0;
        P = P([Nx+x_offset+1:Nx, 1:Nx+x_offset],:,:);
    end

    Ny = size(P,2);
    if y_offset > 0
        P(:,1:y_offset,:)=0;
        P = P(:,[1+y_offset:Ny, 1:y_offset],:);
    elseif y_offset < 0
        P(:,Ny+y_offset+1:Ny,:)=0;
        P = P(:,[Ny+y_offset+1:Ny, 1:Ny+y_offset],:);
    end
elseif  length(size(P))==2
    Nx = size(P,1);
    if x_offset > 0
        P(1:x_offset,:)=0;
        P = P([1+x_offset:Nx, 1:x_offset],:);
    elseif x_offset < 0
        P(Nx+x_offset+1:Nx,:)=0;
        P = P([Nx+x_offset+1:Nx, 1:Nx+x_offset],:);
    end

    Ny = size(P,2);
    if y_offset > 0
        P(:,1:y_offset)=0;
        P = P(:,[1+y_offset:Ny, 1:y_offset]);
    elseif y_offset < 0
        P(:,Ny+y_offset+1:Ny)=0;
        P = P(:,[Ny+y_offset+1:Ny, 1:Ny+y_offset]);
    end
else
    assert(false);
end
    
