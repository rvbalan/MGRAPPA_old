function f=lowpass2D(f,N,c_o)

    f(1:(N/2-c_o),:)=0;
    f(:,1:(N/2-c_o))=0;

    f((N/2+c_o):N,:)=0;
    f(:,(N/2+c_o):N)=0;
