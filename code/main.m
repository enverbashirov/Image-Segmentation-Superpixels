clear;
clc;
close all;

% Image folder location contains collection of .jpgs
location = 'E:\JrJKR\Courses\CS484 - Image Analysis\HW#3\';

% Superpixel Algorithm
sppAlgo = 'slic'; %slic / slic0

contents = dir(strcat(location, 'src\raw\', '*.jpg')); % or whatever the filename extension is
for k = 1:numel(contents)
    k
    filename = contents(k).name;
    folder = contents(k).folder;
    imgRaw = imread( strcat(folder, '\', filename));
    [~,name,~] = fileparts(filename);
    
    %GrayScale Conversion
    %--------------------
    imgGray = rgb2gray(imgRaw);
    imwrite(imgGray, fullfile(strcat(location, 'res\gs\', name, '_gs.jpg')));
    %--------------------
    
    %Reading .dat files
    %------------------
    [imgH,imgW,~] = size(imgRaw);
    SP = transpose(fread(fopen(strcat(location, 'src\', sppAlgo, '\', name, '.dat'), 'r'), [imgW,imgH], 'ulong')) + 1;
    %------------------
    
    %Gabor Convolution
    %----------------
    [gaborC, ~] = gaborconvolve(rgb2gray(imgRaw), 4, 4, 3, 1.7, 0.65, 1.3);

    f = figure('Visible', 'off');
    for t1 = 1:4
        for t2 = 1:4
            subaxis(4,4,((t1-1)*4)+t2, 'Spacing', 0.03, 'Padding', 0, 'Margin', 0);
            imshow(real(gaborC{t1,t2}))
            axis tight
            axis off
        end
    end
    saveas(f, fullfile(strcat(location, 'res\gabor\', name, '_gab.jpg')));
    close(f);
    %----------------
    
    %Color & Texture Feature Vector
    fvPCnt2 = double(zeros(255,3)); %pixel count for averaging for part 4
    fvCol2 = double(zeros(255,3)); %color vector for part 4
    
    maxSPP = max(max(SP)); %maximum superpixel value
    fvPCnt = double(zeros(maxSPP,1)); %pixel count for averaging
    fvCol = double(zeros(maxSPP,3)); %feature vector Color (R,G,B)
    fvTex = double(zeros(maxSPP,size(gaborC,1)*size(gaborC,2))); %feature vector Texture
    
    %Summing RGB channels & Texture Vectors
    for i = 1:imgH
        for j = 1:imgW
            if  SP(i,j) > 0
                %Part 3 - Color Vectors
                fvCol(SP(i,j),1) = fvCol(SP(i,j),1) + double(imgRaw(i,j,1));
                fvCol(SP(i,j),2) = fvCol(SP(i,j),2) + double(imgRaw(i,j,2));
                fvCol(SP(i,j),3) = fvCol(SP(i,j),3) + double(imgRaw(i,j,3));

                %Part 3 - Texture Vectors
                for t1 = 1:size(gaborC,1)
                    for t2 = 1:size(gaborC,2)
                        fvTex(SP(i,j),(((t1-1)*4)+t2)) = fvTex(SP(i,j),(((t1-1)*4)+t2)) + double(abs(gaborC{t1,t2}(i,j)));
                    end
                end
                
                fvPCnt(SP(i,j)) = fvPCnt(SP(i,j)) + 1;
            end
            %Part 4 - Color Vectors (not looking to Superpixels)
            if imgRaw(i,j,1) > 0
                fvCol2(imgRaw(i,j,1),1) = fvCol2(imgRaw(i,j,1),1) + double(imgRaw(i,j,1));
                fvPCnt2(imgRaw(i,j,1),1) = fvPCnt2(imgRaw(i,j,1),1) + 1;
            end
            if imgRaw(i,j,2) > 0
                fvCol2(imgRaw(i,j,2),2) = fvCol2(imgRaw(i,j,2),2) + double(imgRaw(i,j,2));
                fvPCnt2(imgRaw(i,j,2),2) = fvPCnt2(imgRaw(i,j,2),2) + 1;
            end
            if imgRaw(i,j,3) > 0
                fvCol2(imgRaw(i,j,3),3) = fvCol2(imgRaw(i,j,3),3) + double(imgRaw(i,j,3));
                fvPCnt2(imgRaw(i,j,3),3) = fvPCnt2(imgRaw(i,j,3),3) + 1;
            end
        end
    end
    
    %Averaging RGB channels & Texture Vectors
    for i = 1:maxSPP
        for j = 1:3
            fvCol(i,j) = fvCol(i,j) / fvPCnt(i);
        end

        for t1 = 1:size(gaborC,1)
            for t2 = 1:size(gaborC,2)
                fvTex(i,(((t1-1)*4)+t2)) = fvTex(i,(((t1-1)*4)+t2)) / fvPCnt(i);
            end
        end
    end
    for i = 1:255
        fvCol2(i,1) = fvCol2(i,1) / fvPCnt2(i,1);
        fvCol2(i,2) = fvCol2(i,2) / fvPCnt2(i,2);
        fvCol2(i,3) = fvCol2(i,3) / fvPCnt2(i,3);
    end
    %------------------------------
    
    %Forming part 3 a/b/c vectors
    if k == 1
        fvCombCol = fvCol;
        fvCombTex = fvTex;
        fvComb2 = fvCol2;
    else
        fvCombCol = vertcat(fvCombCol, fvCol);
        fvCombTex = vertcat(fvCombTex, fvTex);
        fvComb2 = vertcat(fvComb2, fvCol2);
    end
    
    %----------------------------
end

%-------STEP 3---------
%Normalizing Vectors
fvCombCol = fvCombCol./max(fvCombCol);
fvCombTex = fvCombTex./max(fvCombTex);
fvCombAll = horzcat(fvCombCol, fvCombTex);
fvComb2(isnan(fvComb2))=0;
fvComb2 = fvComb2./max(fvComb2);

%Normalizing Combined Vectors
fvCombCol = ceil(double(fvCombCol./max(fvCombCol)*10000));
fvCombTex = ceil(double(fvCombTex./max(fvCombTex)*10000));
fvCombAll = ceil(double(fvCombAll./max(fvCombAll)*10000));
fvComb2 = ceil(double(fvComb2./max(fvComb2)*10000));
%----------------------

%K-Means Clustering
kCoef = 8; %3 / 5 / 8
p3a = kmeans(fvCombCol, kCoef);
p3b = kmeans(fvCombTex, kCoef);
p3c = kmeans(fvCombAll, kCoef);
p4 = kmeans(fvComb2, kCoef);
%------------------

for k = 1:numel(contents)
    k
    filename = contents(k).name;
    folder = contents(k).folder;
    imgRaw = imread( strcat(folder, '\', filename));
    imgGray = rgb2gray(imgRaw);
    [~,name,~] = fileparts(filename);

    %Reading .dat files
    %------------------
    [imgH,imgW,~] = size(imgRaw);
    SP = transpose(fread(fopen(strcat(location, 'src\', sppAlgo, '\', name, '.dat'), 'r'), [imgW,imgH], 'ulong')) + 1;
    
    if k == 1
        indexSPP = double(0); 
        indexP4 = double(0);
    end
    %------------------
    
    %Exchanging feature vectors with superpixels
    clusterComba = uint8(zeros(imgH, imgW));
    clusterCombb = uint8(zeros(imgH, imgW));
    clusterCombc = uint8(zeros(imgH, imgW));
    clusterComb2 = uint8(zeros(imgH, imgW));
    
    for i = 1:imgH
        for j = 1:imgW
            if SP(i,j) > 0
                clusterComba(i,j) = uint8(p3a(indexSPP + SP(i,j)));
                clusterCombb(i,j) = uint8(p3b(indexSPP + SP(i,j)));
                clusterCombc(i,j) = uint8(p3c(indexSPP + SP(i,j)));
            end
            if imgGray(i,j) > 0
                clusterComb2(i,j) = uint8(p4(indexP4 + double(imgGray(i,j))));
            end
        end
    end
    
    if k == 1
        indexSPP = max(max(SP)); 
        indexP4 = double(max(max(imgGray)));
    else
        indexSPP = indexSPP + max(max(SP));
        indexP4 = double(indexP4 + double(255));
    end
    %--------------------------------------------
    
    f = figure('Visible', 'off');

    subaxis(1,4,1, 'Spacing', 0.03, 'Padding', 0, 'Margin', 0);
    imshow(label2rgb(clusterComba, 'jet'))
    title('Part 3A')
    axis tight
    axis off

    subaxis(1,4,2, 'Spacing', 0.03, 'Padding', 0, 'Margin', 0);
    imshow(label2rgb(clusterCombb, 'jet'))
    title('Part 3B')
    axis tight
    axis off

    subaxis(1,4,3, 'Spacing', 0.03, 'Padding', 0, 'Margin', 0);
    imshow(label2rgb(clusterCombc, 'jet'))
    title('Part 3C')
    axis tight
    axis off

    subaxis(1,4,4, 'Spacing', 0.03, 'Padding', 0, 'Margin', 0);
    imshow(label2rgb(clusterComb2, 'jet'))
    title('Part 4')
    axis tight
    axis off

    saveas(f, fullfile(strcat(location, 'res\', sppAlgo, '\k', int2str(kCoef), '\', name, '_k', int2str(kCoef), 'comb.jpg')));
    close(f);

%     imwrite(label2rgb(clusterComba, 'jet'), fullfile(strcat(location, 'res\', sppAlgo, '\k', int2str(kCoef), '\', name, '_k', int2str(kCoef), 'p3a.jpg')));
%     imwrite(label2rgb(clusterCombb, 'jet'), fullfile(strcat(location, 'res\', sppAlgo, '\k', int2str(kCoef), '\', name, '_k', int2str(kCoef), 'p3b.jpg')));
%     imwrite(label2rgb(clusterCombc, 'jet'), fullfile(strcat(location, 'res\', sppAlgo, '\k', int2str(kCoef), '\', name, '_k', int2str(kCoef), 'p3c.jpg')));
%     imwrite(label2rgb(clusterComb2, 'jet'), fullfile(strcat(location, 'res\', sppAlgo, '\k', int2str(kCoef), '\', name, '_k', int2str(kCoef), 'p4.jpg')));

    %Segment Boundaries as Overlay on Original Image
    mask = boundarymask(clusterComba);
    perim = bwperim(mask, 26); %// Get perimeter of the mask
    red = imgRaw(:,:,1); %// Extract the colour planes of the original image
    green = imgRaw(:,:,2);
    blue = imgRaw(:,:,3);
    red(perim) = 255;
    green(perim) = 255;
    blue(perim) = 255;
    out = cat(3, red, green, blue); %// Make an output image
    imwrite(out, fullfile(strcat(location, 'res\', sppAlgo, '\k', int2str(kCoef), '\', name, '_k', int2str(kCoef), 'p3a.jpg')));
    
    mask = boundarymask(clusterCombb);
    perim = bwperim(mask, 26); %// Get perimeter of the mask
    red = imgRaw(:,:,1); %// Extract the colour planes of the original image
    green = imgRaw(:,:,2);
    blue = imgRaw(:,:,3);
    red(perim) = 255;
    green(perim) = 255;
    blue(perim) = 255;
    out = cat(3, red, green, blue); %// Make an output image
    imwrite(out, fullfile(strcat(location, 'res\', sppAlgo, '\k', int2str(kCoef), '\', name, '_k', int2str(kCoef), 'p3b.jpg')));
    
    mask = boundarymask(clusterCombc);
    perim = bwperim(mask, 26); %// Get perimeter of the mask
    red = imgRaw(:,:,1); %// Extract the colour planes of the original image
    green = imgRaw(:,:,2);
    blue = imgRaw(:,:,3);
    red(perim) = 255;
    green(perim) = 255;
    blue(perim) = 255;
    out = cat(3, red, green, blue); %// Make an output image
    imwrite(out, fullfile(strcat(location, 'res\', sppAlgo, '\k', int2str(kCoef), '\', name, '_k', int2str(kCoef), 'p3c.jpg')));    
    
    mask = boundarymask(clusterComb2);
    perim = bwperim(mask, 26); %// Get perimeter of the mask
    red = imgRaw(:,:,1); %// Extract the colour planes of the original image
    green = imgRaw(:,:,2);
    blue = imgRaw(:,:,3);
    red(perim) = 255;
    green(perim) = 255;
    blue(perim) = 255;
    out = cat(3, red, green, blue); %// Make an output image
    imwrite(out, fullfile(strcat(location, 'res\', sppAlgo, '\k', int2str(kCoef), '\', name, '_k', int2str(kCoef), 'p4.jpg')));
    %-----------------------------------------------
end