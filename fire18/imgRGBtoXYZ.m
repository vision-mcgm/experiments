function [ img2 ] = imgRGBtoXYZ( img )
%Converts an image from RGB to XYZ. Assumes it's an uint, so will have
%0-255 color range.

[nil,M]=XYZToSRGBPrimary([]); %Get conv matrix
img=double(img)/255;
[h w c]=size(img);
M=inv(M); %Invert matrix for conversion
it=1;
img2=zeros(h, w, c);
    for iw=1:w
        img2(:,iw,:)=squeeze(img(:,iw,:))*M;
       % img2(ih,iw,:)=M*squeeze(img(ih,iw,:));
    end


end

