function tf=IsEndPoint(cubj)


tf=false;
nb=[32, 33, 34, 37, 38, 39, 42, 43, 44, 57,...
    58, 59, 62, 63, 64, 67, 68, 69, 82, 83,...
    84, 87, 88, 89, 92, 93, 94];
conn6=false(3,3,3);
conn6([5 11 13 15 17 23])=true;
cubjj=cubj(2:4,2:4,2:4);
cubjj(14)=false;
cmorph=bwmorph3(cubjj,'clean');

nbdiff=nb-63;

k=find(cubj);

C=intersect(k,nb);
C(C==63)=[];

if numel(C)==1
    tf=true;
elseif numel(C)==2
    for i=1:2
        testidx=C(i);
        nbtest=testidx-nbdiff;
        Ctest=intersect(k,nbtest);
        if numel(Ctest)==3 && numel(intersect(Ctest,nb))==3
            tf=true;
        end
    end
elseif numel(C)==3 && sum(cmorph(:))==3
    tf=true;
end