function imshowSc(P)

P=squeeze(P);
imshow(abs(P),[0,max(max(abs(P)))]);
