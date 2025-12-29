function [Icleantrcrop, NWsk, DWcured, DistT, FinalT] = eteetsNeurons(tiffile)
% eteetsNeurons computes metrics related to neurons-receptive endings of
% AFD neurons from 3D stacks collected on vtiSIM
%
% INPUT:
%       tiffile: full path of a '.tif' file name
%
% OUTPUT:
%       Icleantrcrop: processed, cropped image matrix
%       NWsk: binary stack of skeletonized object resulting from the adaptive resolution orientation space segmentation
%       DWcured: Binary stack containing NRE endpoints
%       DistT: table containing computed distances
%       FinalT: table containing computed properties of the NRE

%% get image info and load image
try %been used to prevent termination in case of a corrupted image during batch processing. Remove the try catch end operation if debugging
    [~, ~, sinfo] = bfGetInfo(tiffile); % reading metadata
    Xvol = sinfo.PhysicalXY;
    Zvol = sinfo.PhysicalZ;

    IS = bfopen(tiffile);
    Istack = cat(3, IS{1,1}{:,1});

    %% Preprocess image
    [~, ~, nZ] = size(Istack);
    Iclean = zeros(size(Istack), 'like', Istack);


    for i = 1:nZ
        I = Istack(:,:,i);
        Iclean(:,:,i) = I - imgaussfilt(I, 50);
    end


    Iclean(Iclean < 0) = 0;
    T = minminT(max(Iclean, [], 3));
    Icleant = Iclean - T;
    Icleant(Icleant < 0) = 0;
    Icleant = single(Icleant);
    Icleantr = Icleant / max(Icleant(:));

    Icleantr(:,:,1:3) = []; %remove the first 3 stacks that often contain data from another neuron
    nZ = size(Icleantr, 3);


    %% crop to ROI
    IW = imbinarize(Icleantr);
    NWD = imclearborder(imdilate(IW, strel('disk', 5)));
    if sum(NWD(:)) > sum(IW(:))
        S = regionprops3(double(NWD), 'BoundingBox');
    else
        S = regionprops3(double(IW), 'BoundingBox');
    end
    rect = S.BoundingBox;
    rect(3) = 1; rect(end) = nZ - 1;
    Icleantrcrop = imcrop3(Icleantr, rect);

    %% get weighted centroid of object
    Scent = regionprops3(double(imbinarize(Icleantrcrop)), Icleantrcrop, 'WeightedCentroid');
    centm = round(Scent.WeightedCentroid);
    M = false(size(Icleantrcrop));
    M(centm(2),centm(1),centm(3)) = true;
    M = imdilate(M, strel('disk', 10)); %creates a sphere to connect vilii when computing geodesic distances


    %% AROS computation and skeletonization
    NMS = zeros(size(Icleantrcrop));

    f = waitbar(0, 'Performing AROS segmentation...');
    tic
    for i = 1:nZ
        Iz = Icleantrcrop(:,:,i);
        if ~any(Iz,"all")
            nms = zeros(size(Iz));
        else
            [~, ~, nms] = steerableAdaptiveResolutionOrientationSpaceDetector...
                ( double(Iz));
        end
        NMS(:,:,i) = nms;
        waitbar(i / nZ);
    end
    toc
    close(f);
    NW = imbinarize(NMS);
    NWsk = bwskel(NW, 'MinBranchLength', 10);
    if sum(NWsk(:)) == 0
        Icleantrcrop = [];
        NWsk = [];
        DW = [];
        Volume = NaN;
        tvolume = table(Volume);
        Centroid = nan(1,3);
        tcentroid = table(Centroid);
        EquivDiameter = NaN;
        ted = table(EquivDiameter);
        Extent = NaN;
        textent = table(Extent);
        PrincipalAxisLength = [NaN NaN NaN];
        tpal = table(PrincipalAxisLength);
        ConvexVolume = NaN;
        tcv = table(ConvexVolume);
        Solidity = NaN;
        tsol = table(Solidity);
        Sample = cellstr(tiffile);
        tsample = table(Sample);
        Distance = NaN;
        tdistance = table(Distance);
        geoDistance = NaN;
        tgeodistance = table(geoDistance);
        nearestDistance = NaN;
        tnearestdistance = table(nearestDistance);

        DistT = [tdistance tgeodistance tnearestdistance];
        FinalT = [tsample tvolume tcentroid ted textent tpal tcv tsol];
    else
        % connecting vilii to spherical center and increas
        CW = NWsk & M;

        while ~any(CW,'all')
            M = imdilate(M, ones(3));
            CW = NWsk & M;
        end

        %% get endpoints and extract metrics
        EP = bwmorph3(NWsk, 'endpoints');

        Ssk = regionprops3(NWsk, 'Centroid', 'PrincipalAxisLength',...
            'Volume', 'ConvexVolume', 'Extent', 'Solidity', 'EquivDiameter');
        tf = Ssk.Volume == max(Ssk.Volume);

        Ssk(~tf,:) = [];

        Sample = cellstr(tiffile);
        tsample = table(Sample);
        FinalT = [tsample Ssk];
        FinalT.Volume = FinalT.Volume * Xvol^2 * Zvol;
        FinalT.ConvexVolume = FinalT.ConvexVolume * Xvol^2 * Zvol;
        FinalT.EquivDiameter = FinalT.EquivDiameter * Xvol;

        D = bwdistgeodesic(NWsk, M);
        Drmax = D;


        Drmax(isnan(Drmax)) = 0; Drmax(Drmax == Inf) = 0;
        DW = imregionalmax(Drmax);
        DW = DW | EP;


        EPVoxList = eteetscureEndpointsv2(NWsk); %screen each pixel for correct endoint configuration

        Dgeo = D(EPVoxList);

        Dgeo = Dgeo * Xvol;

        Dcart = bwdist(M);
        PD = Dcart(EPVoxList);
        PD = PD * Xvol;

        [r, c, z] = ind2sub(size(DW), EPVoxList);
        cent = [r c z];
        [~, DD] = knnsearch(cent,cent,'K',2);
        nearD = DD(:,2) * Xvol;

        Dist = [PD, Dgeo nearD];
        DistT = array2table(Dist, 'VariableNames',...
            {'Distance', 'geoDistance', 'nearestDistance'});
        DWcured = false(size(DW));
        DWcured(EPVoxList)=true;


        %% Plot


        figure; imshow(max(Icleantr, [], 3), []);

        rendvol = viewer3d;
        fig = rendvol.Parent;
        rendvol.BackgroundGradient = 'off';
        rendvol.BackgroundColor = [1 1 1];
        NWski = imdilate(NWsk, ones(3));
        vsk = volshow(NWski,Parent= rendvol);
        vsk.Transformation.A(3,3) = sinfo.PhysicalZ / sinfo.PhysicalXY;
        DWi = imdilate(DWcured,ones(3));
        vdw = volshow(DWi, Colormap = [1 0 0], Parent = rendvol);
        vdw.Transformation.A(3,3) = sinfo.PhysicalZ / sinfo.PhysicalXY;

        pause;
        delete(fig)

        close all;
    end

catch
    disp(['Unable to process file ' tiffile]);
    Icleantrcrop = [];
    NWsk = [];
    DW = [];
    Volume = NaN;
    tvolume = table(Volume);
    Centroid = nan(1,3);
    tcentroid = table(Centroid);
    EquivDiameter = NaN;
    ted = table(EquivDiameter);
    Extent = NaN;
    textent = table(Extent);
    PrincipalAxisLength = [NaN NaN NaN];
    tpal = table(PrincipalAxisLength);
    ConvexVolume = NaN;
    tcv = table(ConvexVolume);
    Solidity = NaN;
    tsol = table(Solidity);
    Sample = cellstr(tiffile);
    tsample = table(Sample);
    Distance = NaN(1);
    tdistance = table(Distance);
    geoDistance = NaN(1);
    tgeodistance = table(geoDistance);
    nearestDistance = NaN(1);
    tnearestdistance = table(nearestDistance);

    DistT = [tdistance tgeodistance tnearestdistance];
    FinalT = [tsample tvolume tcentroid ted textent tpal tcv tsol];

end
