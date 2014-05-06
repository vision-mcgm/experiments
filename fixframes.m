function [  ] = fixframes(  )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
f='D:\Fintan\faceframes\';

ld=dir(f);

nd=size(ld,1);

for i=1:nd
    nm=ld(i).name;
    if ~strcmp(nm,'.') & ~strcmp(nm,'..') & isdir([f nm])
   %Inner loop
   fullnm=[f nm];
   
   frames=dir([fullnm '\*.bmp']);
   [~,order]=sort_nat({frames.name});
   newframes=frames(order);
   nfr=size(newframes,1);
   for f=1:nfr
       % ['frame' num2str(f,'%05d') '.bmp']
       fprintf('%d\n',f);
       movefile([fullnm '\' newframes(f).name], [fullnm '\frame' num2str(f,'%05d') '.bmp']);
   end
end



end

