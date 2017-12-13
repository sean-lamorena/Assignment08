clear;
clc;
close all;



addpath(genpath('MatlabFns'));

mkdir('_output');

StartingFrame = 1;
EndingFrame = 230;

directory = 'img1/img';

for k = StartingFrame : EndingFrame
    
Ipts1 = [];
ipts2 = [];

D1 = [];
D2 = [];

BaseIndex = [];
SubIndex = [];
index = [];

BaseLength = [];
SubLength = [];

SubValue = [];
    
    % Reading in images.. if past the first frame, use the Jregistered
    % image from the previous step as pos1, and load a new image into pos2
    if k == 1
        pos1 = imread([directory, sprintf('%3.3d',k),'.jpg']);
        pos2 = imread([directory, sprintf('%3.3d',k+1),'.jpg']);
    else
        pos1 = Jregistered;
        pos2 = imread([directory, sprintf('%3.3d',k+1),'.jpg']);
    end
    
    % Convert images to grayscale (necessary for Harris corner detection)
    Im1 = rgb2gray(pos1);
    Im2 = rgb2gray(pos2);  

    % SURF portion to help with Stabilization

    Ipts1 = OpenSurf(Im1);
    Ipts2 = OpenSurf(Im2);

    for k = 1:length(Ipts1)
    D1(:,k) = Ipts1(k).descriptor;
end

for k = 1:length(Ipts2)
    D2(:,k) = Ipts2(k).descriptor;
end

BaseLength = length(Ipts1);
SubLength = length(Ipts2);

for i = 1:BaseLength
    subtract = (repmat(D1(:,i), ...
        [1 SubLength]) - D2).^2;
    distance = sum(subtract);
    [SubValue(i) SubIndex(i)] = min(distance);
end

[value, index] = sort(SubValue);

index = index(1:100);

BaseIndex = index;
SubIndex = SubIndex(index);

m1 = floor([[Ipts1(BaseIndex).y]; ...
    [Ipts1(BaseIndex).x]]);

m2 = floor([[Ipts2(SubIndex).y]; ...
    [Ipts2(SubIndex).x]]);

    % RANSAC
    [H, inliers] = ransacfithomography(m1, m2, 0.001);
    
    % Moving = second frame, Fixed = first frame or Jregistered
    fixedPoints = [m1(2,inliers)' m1(1,inliers)'];
    movingPoints = [m2(2,inliers)' m2(1,inliers)'];
     
    % Determining the transform based on the relationship matrices between
    % the coordinates in the two images
    tform = fitgeotrans(movingPoints,fixedPoints,'NonreflectiveSimilarity');
    
    % Image registration (alignment)
    Jregistered = imwarp(pos2,tform,'OutputView',imref2d(size(pos1)));
    falsecolorOverlay = imfuse(pos1,Jregistered);

    % Putting everything together into a 2x2 image
    I1 = cat(2,pos2,Jregistered);
    im1rgb = cat(3,Im1,Im1,Im1);
    I2 = cat(2,im1rgb,falsecolorOverlay);
    I = cat(1,I1,I2);

    % Accounting for the concatenation
    shiftY = size(Im1,1);
    shiftX = size(Im1,2);

    % Displaying the four images
    imshow(I,'Border','tight'); hold on;

    % Plotting the relationships between POS1 and POS2 images (RANSAC
    % results only at this point)
    plot(m1(2,inliers),m1(1,inliers)+shiftY,'r+');
    plot(m2(2,inliers),m2(1,inliers)+shiftY,'b+');    

    for n = inliers
        line([m1(2,n) m2(2,n)], [m1(1,n)+shiftY m2(1,n)+shiftY],'color',[1 1 0])
    end

    
    % Text and other stuff
    
    white = [0.99,0.99,0.99];

    text(25,25,'ORIGINAL ','fontsize',16,'color',white,'fontweight','bold');
    text(25+shiftX,25,'STABILIZED using RANSAC ','fontsize',16,'color',white,'fontweight','bold');

    text(25,25+shiftY,'INLYING MATCHES ','fontsize',16,'color',white,'fontweight','bold');
    text(25+shiftX,25+shiftY,'REGISTRATION ERROR ','fontsize',16,'color',white,'fontweight','bold');

    frameNumber = sprintf('%3.3d',k);
    text(25,shiftY-50,'f-no:','fontsize',16,'color',white,'fontweight','bold');
    text(75,shiftY-50,frameNumber,'fontsize',16,'color',white,'fontweight','bold');

    inliersText = sprintf('Inliers: %d (%d%%) \n',length(inliers),round(100*length(inliers)/length(m1)));
    putativeText = sprintf('Putative matches: %d \n', length(m1));

    text(25,2*shiftY-50,inliersText,'fontsize',16,'color',white,'fontweight','bold');
    text(25,2*shiftY-25,putativeText,'fontsize',16,'color',white,'fontweight','bold');

    text(25,shiftY-25,'A. Almagambetov','fontsize',16,'color',white,'fontweight','bold');

    % Display everything to the figure window
    drawnow;

    % Store the figure window contents into a .jpg file (uncomment if
    % necessary)
    %screen2jpeg(['_output/img',frameNumber,'.jpg']);

    % Clear figure window of all content (prevents MATLAB from getting
    % progressively slower with each displayed image)
    clf;
end