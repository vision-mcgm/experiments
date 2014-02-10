function [ img2 ] = imgXYZtoRGB( img )
%Converts an image from RGB to XYZ. Assumes it's an uint, so will have
%0-255 color range.

[nil,M]=XYZToSRGBPrimary([]); %Get conv matrix
if strcmp(class(img),'uint8')
img=double(img)/255;
end
[h w c]=size(img);

it=1;
img2=zeros(h, w, c);
    for iw=1:w
        img2(:,iw,:)=squeeze(img(:,iw,:))*M;
       % img2(ih,iw,:)=M*squeeze(img(ih,iw,:));
    end


end

