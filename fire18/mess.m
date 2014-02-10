function [ img ] = mess( img )
%Complements the hue channel of an HSV image
img=rgb2hsv(img);
%img(:,:,1)=mod(img(:,:,1)+0.5,1);
img(:,:,3)=1-img(:,:,3);
img=hsv2rgb(img);
img=uint8(img*255);

end

