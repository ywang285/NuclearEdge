clear all
close all
mkdir('AfterAnalysis')
filelist=dir('*.tif');
edgethreshold=1; % default
for i=1:size(filelist,1)
filename=filelist(i).name(1:end-4);
TIFfilename=strcat(filename,'.tif');
outputDlist=strcat(filename,'_Dlist.txt');
outputThreshold=strcat(filename,'_InputThreshold.csv');
outputEdgeDlist=strcat(filename,'_EdgeDlist.dat');
outputNucleusEdge=strcat(filename,'_NucleusEdge.png');
outputBinarization=strcat(filename,'_Binarization.png');
outputFISHspots=strcat(filename,'_FISHspots.png');
%% load image
Channel1(:,:) = imread(TIFfilename,1); % channel1: Cy5
Channel2(:,:) = imread(TIFfilename,2); % channel2: Cy3
Channel3(:,:) = imread(TIFfilename,3); % channel3: Hoechst
unit8channel3=uint8(double(Channel3)/max(double(Channel3(:)))*50);
%%find edges
loop=0;
[outputImage, number_of_nuclei, afterLog,binarization0] = nuclei_counter(unit8channel3);
prompt = strcat('Is the edge finding good (if yes type enter, if no input new threshold (start with 1, higher number corresponds higher threshold))?');
Good_edgethreshold = input(prompt,'s');
if isempty(Good_edgethreshold)
    loop=1;
else
    edgethreshold = str2double(Good_edgethreshold);
end
while loop==0
[outputImage, number_of_nuclei, afterLog,binarization0] = nuclei_counter(unit8channel3,edgethreshold/10);
prompt = strcat('Is the edgethreshold=', num2str(edgethreshold), ' good(if yes type enter, if no input new threshold)?');
Good_edgethreshold = input(prompt,'s');
if isempty(Good_edgethreshold)
    loop=1;
else
    edgethreshold = str2double(Good_edgethreshold);
end
end
s = regionprops(binarization0,'centroid');
centroids = cat(1,s.Centroid);
hold on
plot(centroids(:,1),centroids(:,2),'b*')
hold on
outputNucleusEdgefig=imshow(outputImage, [min(unit8channel3(:)) max(unit8channel3(:))]);
binarization1=imfill(binarization0,'holes');
binarization = bwareafilt(binarization1,1);
% figure;imshowpair(binarization0,binarization,'montage')
%% Find ellipse fit and major/minor axes
s0 = regionprops(binarization,{'Centroid','Orientation','MajorAxisLength','MinorAxisLength'});
s=s0(1);
% Calculate the ellipse line
theta = linspace(0,2*pi);
col = (s.MajorAxisLength/2)*cos(theta);
row = (s.MinorAxisLength/2)*sin(theta);
M = makehgtform('translate',[s.Centroid, 0],'zrotate',deg2rad(-1*s.Orientation));
D = M*[col;row;zeros(1,numel(row));ones(1,numel(row))];
%Visualize the result
figure
imshow(binarization)
hold on
plot(D(1,:),D(2,:),'r','LineWidth',2)
centroids = cat(1,s.Centroid);
hold on
plot(centroids(:,1),centroids(:,2),'b*')
hold on
deltaX1=[s.Centroid(1)-0.5*s.MajorAxisLength * cosd(-1*s.Orientation) s.Centroid(1)+0.5*s.MajorAxisLength * cosd(-1*s.Orientation)]; % x coordinates of two edges of the major axis
deltaY1=[s.Centroid(2)-0.5*s.MajorAxisLength * sind(-1*s.Orientation) s.Centroid(2)+0.5*s.MajorAxisLength * sind(-1*s.Orientation)]; % y coordinates of two edges of the major axis
MajorAxis1=[s.Centroid(1)-0.5*s.MajorAxisLength * cosd(-1*s.Orientation), s.Centroid(2)-0.5*s.MajorAxisLength * sind(-1*s.Orientation),0]; % (x,y) coordinates of the left edges of the major axis
MajorAxis2=[s.Centroid(1)+0.5*s.MajorAxisLength * cosd(-1*s.Orientation), s.Centroid(2)+0.5*s.MajorAxisLength * sind(-1*s.Orientation),0]; % (x,y) coordinates of the right edges of the major axis
deltaX2=[s.Centroid(1)-0.5*s.MinorAxisLength * sind(s.Orientation) s.Centroid(1)+0.5*s.MinorAxisLength * sind(s.Orientation)]; % x coordinates of two edges of the minor axis
deltaY2=[s.Centroid(2)-0.5*s.MinorAxisLength * cosd(s.Orientation) s.Centroid(2)+0.5*s.MinorAxisLength * cosd(s.Orientation)]; % y coordinates of two edges of the minor axis
MinorAxis1=[s.Centroid(1)-0.5*s.MinorAxisLength * sind(s.Orientation), s.Centroid(2)-0.5*s.MinorAxisLength * cosd(s.Orientation),0]; % (x,y) coordinates of the left edges of the minor axis
MinorAxis2=[s.Centroid(1)+0.5*s.MinorAxisLength * sind(s.Orientation), s.Centroid(2)+0.5*s.MinorAxisLength * cosd(s.Orientation),0]; % (x,y) coordinates of the right edges of the minor axis
%% find FISH spots center of mass
figure;imshow(Channel1, [min(Channel1(:)) max(Channel1(:))])
prompt = strcat('How many FISH spots in the cell? (default=1 if press enter)');
N_FISHspot = input(prompt);
if isempty(N_FISHspot)
    N_FISHspot=1; % default FISH spot number
end
xinput=[];
yinput=[];
TempImage=double(Channel1).*double(binarization);
for n_FISHspot=1:N_FISHspot
[y_max,x_max]=find(TempImage==max(TempImage(:)));
xinput=[xinput;x_max];
yinput=[yinput;y_max]; 
TempImage((y_max-3:y_max+3),(x_max-3:x_max+3))=zeros(7);
end
figure;imshow(Channel1, [min(Channel1(:)) max(Channel1(:))])
for n_FISHspot=1:N_FISHspot
hold on
plot(xinput(n_FISHspot),yinput(n_FISHspot),'r*')
end
prompt = strcat('FISH-spot found correctly? (if yes press any key, if no, press n, then manully locate the FISH spot)');
FISHspotFinding = input(prompt,'s');
if FISHspotFinding == 'n'
    xinput=[];
    yinput=[];
    figureginput=imshow(Channel1, [min(Channel1(:)) max(Channel1(:))]);
    %imcontrast
    [xinput,yinput] = ginput;
end
CoMlist=[];
for n=1:size(xinput)
size_of_fitting_window=10;
FittingArea=Channel1(round(yinput(n)-size_of_fitting_window/2):round(yinput(n)+size_of_fitting_window/2),...
    round(xinput(n)-size_of_fitting_window/2):round(xinput(n)+size_of_fitting_window/2));
[Y, X]=find(FittingArea>prctile(FittingArea(:), 100*(1-5/size_of_fitting_window^2)));
Y_original=Y+round(yinput(n)-size_of_fitting_window/2)-1; % row
X_original=X+round(xinput(n)-size_of_fitting_window/2)-1; % column
%org_index= range(index);
%imshow(FittingArea, [min(FittingArea(:)) max(FittingArea(:))])

sumX=0;
sumY=0;
sumI=0;
for i=1:size(Y_original)
    sumX=sumX+double(Channel1(Y_original(i),X_original(i)))*X_original(i);
    sumY=sumY+double(Channel1(Y_original(i),X_original(i)))*Y_original(i);
    sumI=sumI+double(Channel1(Y_original(i),X_original(i)));
end
CoM=[sumX/sumI sumY/sumI 0];
CoMlist=[CoMlist;CoM];
end

%% Plot
Binarizationfig=figure;
imshow(binarization)
hold on
plot(D(1,:),D(2,:),'r','LineWidth',2)
FISHspotsfig=figure;
imshow(Channel1, [min(Channel1(:)) max(Channel1(:))])
%imcontrast
DList=[];
EdgeDList=[];
hold on
plot(D(1,:),D(2,:),'r','LineWidth',2)
hold on
line(deltaX1,deltaY1,'color','b','LineWidth',1.5)
hold on
line(deltaX2,deltaY2,'color','b','LineWidth',1.5)
[y_row_edge,x_column_edge]=find(afterLog);
for k=1:size(xinput)
hold on
plot(CoMlist(k,1),CoMlist(k,2),'m+','MarkerSize',10,'LineWidth',1)
FISH_point_To_Major_Y=point_to_line_distance(CoMlist(k,:),MajorAxis1,MajorAxis2);
FISH_point_To_Minor_X=point_to_line_distance(CoMlist(k,:),MinorAxis2,MinorAxis1);
FISHspot_to_edge_distance=min(sqrt((x_column_edge-CoMlist(k,1)).^2 + (y_row_edge-CoMlist(k,2)).^2));
EdgeDList=[EdgeDList;FISHspot_to_edge_distance/sqrt(sum(binarization(:)))];
DList=[DList;FISH_point_To_Minor_X/(s.MajorAxisLength/2) FISH_point_To_Major_Y/(s.MinorAxisLength/2)]; % Coordinate
end
movefile(TIFfilename,'AfterAnalysis');
cd('AfterAnalysis')
%writematrix(DList,outputDlist,'Delimiter','tab')
writematrix(EdgeDList,outputEdgeDlist,'Delimiter','tab')
%writematrix(edgethreshold,outputThreshold,'Delimiter','tab')
%saveas(outputNucleusEdgefig,outputNucleusEdge);
%saveas(Binarizationfig,outputBinarization);
%saveas(FISHspotsfig,outputFISHspots);
cd ..
% close all
% clearvars -except filelist edgethreshold
Channel1=[];
Channel2=[];
Channel3=[];
end