function EPVoxList=eteetscureEndpointsv2(NWsk)
% eteetscurEndpointsv2 scans each pixel neighborhood from a 5x5x5 matrix

sdw=size(NWsk);

NWsk([1:2 sdw(1)-2:sdw(1)],:,:)=false;
NWsk(:,[1:2 sdw(2)-2:sdw(2)],:)=false;
NWsk(:,:,[1:2 sdw(3)-2:sdw(3)])=false;

VoxList=find(NWsk);
cubcent=false(5,5,5);
cubcent(63)=true;
TF=false(numel(VoxList),1);
[x, y, z]=ind2sub(sdw,VoxList);

for i=1:numel(VoxList)


    
    a=x(i)-2;b=x(i)+2;
    c=y(i)-2;d=y(i)+2;
    e=z(i)-2;f=z(i)+2;
        cub=NWsk(a:b,c:d,e:f);
        cubj=imreconstruct(cubcent,cub);
        TF(i)=IsEndPoint(cubj);

end



EPVoxList=VoxList(TF);





