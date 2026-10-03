function [m,x,y] = max2D(A)

[mx,ix]=max(A);
[m,iy]=max(mx);
x=ix(iy);
y=iy;
