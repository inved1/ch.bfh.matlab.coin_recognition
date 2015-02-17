% Matlab Project by Julien Villiger (villj2@bfh.ch) and Daniel Inversini
% (inved1@bfh.ch)
% 
% Target: Recognize coins (swiss francs) from a picture and count the money
% How: search with Hough
% Secondly: Pattern matching
%

% Steps:
% 1. Get Settings
% 2. Read Image and create Grayscale Image 
% 3. Create Histogram
% 4. Create Binary File
% 5. Holes Filling
% 6. Closing
% 7. Diameter
% 8. Edge detection
% 9. Hough Transformation
% 10. 

clc; clear; clear all; clear workspace;


tic % T T T TIMER

% 1.
% --- get Settings (csv, first column is Amount (CHF), second is diameter --- %
mySettings = coin_recognition_settings();
% this is the 0.01 CHF coin, reference coin
MassReference = 16;


% 2. 
% Read Image and create Grayscale Image %
X = imread('img_coins.jpg');
XPattern = X; %save
I = im2double(X);
IGray = rgb2gray(I);

%figure; imshow(I); title('Original');
%figure; imshow(IGray); title('Graustufenbild');

% 3.
% Create Histogram %
Hist = imhist(IGray);

%figure; imhist(IGray); title('Histogramm Graustufenbild');



% 4. 
% Create binary File %
IBinary = im2double(im2bw(I, 0.15));
IBinaryAllCoins = im2double(im2bw(I, 0.9));

%figure; imshow(IBinaryAllCoins); title('Binaeres Bild');

% 5. 
% Holes Filling %
IBinary = imcomplement(IBinary);
IBinary = imfill(IBinary,'holes');
IBinary = imcomplement(IBinary);

IBinaryAllCoins = imcomplement(IBinaryAllCoins);
IBinaryAllCoins = imfill(IBinaryAllCoins,'holes');
IBinaryAllCoins = imcomplement(IBinaryAllCoins);


%figure; imshow(IBinaryAllCoins); title('Holes Filling');

% 6. 
% Closing %
se = strel('disk',3);
IBinary = imclose(IBinary, se);
IBinaryAllCoins = imclose(IBinaryAllCoins, se);

%figure; imshow(IBinaryAllCoins); title('Closing');
IMGpattern = IBinaryAllCoins; % save for later


% 7.
% Get Diameter
% Count white Pixels %
PixelsWhite = sum(IBinary(:)==0);
Diameter = 2 * sqrt(PixelsWhite / pi);


% 8.
% Edge Detection %
IEdge = edge(IBinary,'canny');
IEdgeAllCoins = edge(IBinaryAllCoins, 'canny');

%figure; imshow(IEdge); title('Kantendetektion');
%figure; imshow(IBinaryAllCoins); title('Binaeres Bild aller Muenzen');
%figure; imshow(IEdgeAllCoins); title('Kantendetektion alle Muenzen');

% Calculate Searched Diameters
SearchedDiameters = {
    getPercentage(MassReference, mySettings(8)) * Diameter
    getPercentage(MassReference, mySettings(9)) * Diameter
    getPercentage(MassReference, mySettings(10)) * Diameter
    getPercentage(MassReference, mySettings(11)) * Diameter
    getPercentage(MassReference, mySettings(12)) * Diameter
    getPercentage(MassReference, mySettings(13)) * Diameter
    getPercentage(MassReference, mySettings(14)) * Diameter
};


% 9.
% Hough Transformation
% Now search Coins with class HoughCircle
money = 0;
for d = 1:length(SearchedDiameters)
    
    SearchedRadius = SearchedDiameters{d} / 2;
    
    foundPositionsX = [];
    
    for x = -2:2
        [y0detect,x0detect,Accumulator] = HoughCircle2(IEdgeAllCoins, int32(SearchedRadius) + x, (SearchedRadius + x) * pi);
        
        % Transpose Matrix to avoid errors in groupResults
        x0detect = x0detect';
        
        x0detectGrouped = groupResults(x0detect, 10);
        
        % TODO group array and check if previous value was similiar
        
        % 2 Fehlerquellen:
        % 
        % 1) Ein detect-array kann so aussehen: [412, 413, 415, 912] ->
        % zusammenfassen! also [413, 912]
        % -> BEHOBEN :D
        %
        % 2) bei einem neuen x kann dieselbe m?nze wieder erkannt werden.
        % darum ?hnliche resultate vergessen bei weiteren loops!
        % -> Also ein Array f?hren mit gefunden xDetect-Werten.
        % -> Falls bei einem weiteren Loop ein ?hnlicher Wert vorkommt,
        % also dieselbe M?nze nochmals gefunden wird, einfach vergessen!
        % -> BEHOBEN :D
        
        
        
        %average = x0detect;
        
        x0detectGroupedFinal = [];
        
        % Loop through all Groups
        for v = 1:size(x0detectGrouped',1)
            
            arr = x0detectGrouped{v};
            similiarCumulated = 0;
            
            % Loop through all Arrays in Groups
            for v1 = 1:size(arr',1)
                similiarCumulated = similiarCumulated + arr(v1);
            end
            
            average = similiarCumulated / size(arr',1);
            
            if(checkExisting(foundPositionsX, average, 5))
            else
                x0detectGroupedFinal = [x0detectGroupedFinal, average];
            end
        end
        
        
        
        % todo: check if current average is already (or nearly) in foundPositionsX
        % true: dont add average and skip this loop
        % false: add average and add money
        
        % ismember([15 17],X)
        
%         if(checkExisting(foundPositionsX, average, 5))
%            continue; 
%         end
        
        %for e = 1:size(x0detectGrouped',1)
        for e = 1:size(x0detectGroupedFinal',1)    
            
            %if(x0detectGrouped{e})
            %    if(x0detectGrouped{e})
            
            if(x0detectGroupedFinal(e))
                if(x0detectGroupedFinal(e))
                
                    % coin found!
                    % now count money

                    %coinAmount = size(x0detect,1)
                    %money = money + coinAmount * CoinValues{d};
                    
                    money = money + mySettings(d:d);
                    
                    % Add Found Value (average) to Array
                    foundPositionsX = [foundPositionsX, x0detectGroupedFinal(e)];
                end
            end
        end
        
        % Add Found Value (average) to Array
        % foundPositionsX = [foundPositionsX, average];
    end
end

disp('Final Money: CHF ')
disp(money)
toc

money = 0;

disp('Pattern Matchin Start')
tic
%
% 1. get all patterns
%
fnames = dir('ref_coins/*.jpg');  
for i = 1 : length(fnames)
    filename = strcat('ref_coins/',fnames(i).name);
    k = strfind(filename, '_');
    j = strfind(filename,'/');
    ref_coin_value = filename(j+1:k(end)-1); % get last index of _ with keyword end 
    %cannot im2double do here, because i need to stretch them first 
    %ref_coins{i} = {ref_coin_value,im2double(imread(filename))};
    
    ref_coins{i} = {ref_coin_value,imread(filename)};
    %figure, imshow(I);
end

cnt_ref_coins = numel(ref_coins);

% 
% 2. get all areas
% 
%IMGpattern = imcomplement(IMGpattern); %swap b&w without comment it is
%white backgrou nd
%figure; imshow(IMGpattern); title('test dani');

PatternAttributes = regionprops(IEdgeAllCoins,'centroid', 'area', 'BoundingBox');
cntPatternAttributes = size(PatternAttributes);

% 
% 3. Crop found coins/areas
% 
for i = 1 : cntPatternAttributes
    bboxArea = PatternAttributes(i).BoundingBox;  
    %increase bounding box by 100 px each side, so that template is always
    %smaller than this cropped area
    bboxArea(1) = bboxArea(1) - 150;
    bboxArea(2) = bboxArea(2) - 150;
    bboxArea(3) = bboxArea(3) + 300;
    bboxArea(4) = bboxArea(4) + 300;
	subImage = imcrop(XPattern, bboxArea);
    diameter = sqrt(PatternAttributes(i).Area / pi) *2 ;
    img_coins{i} = {diameter, im2double(rgb2gray(subImage))};
    %figure; imshow(subImage);
end

% 
% 4. now match the ref_coins with the img_coins
% 

cnt_img_coins = numel(img_coins);

%debug
%cnt_img_coins = 1;

iFoundCoins = 1;

for i = 1 : cnt_img_coins
    
    % now: 
    % a) calculate diameter from img_coin
    diameter = img_coins{1,i}(1);
    diameter = diameter{1};
    
    IMGCoin2Search = img_coins{1,i}(2);
    IMGCoin2Search = IMGCoin2Search{:};
    %disp(i)
    %disp(diameter)
    % b) loop trough ref_coins
   
    DIMGCoin2Search = im2double(IMGCoin2Search);
    
    %debug
    %cnt_ref_coins = 1 ;
    
    for j = 1 : cnt_ref_coins
        
        % c) stretch/adjust ref_coin(j) image to this diameter
        unstretchedRefIMG = ref_coins{1,j}(2);
        unstretchedRefIMG = rgb2gray(unstretchedRefIMG{1});
        
        ref_coin_name = ref_coins{1,j}(1);
        
        %create binary file for area later
        IBinaryUnstrechedRefIMG = im2double(im2bw(unstretchedRefIMG, 0.9));
        %filling
        IBinaryUnstrechedRefIMG = imcomplement(IBinaryUnstrechedRefIMG);
        IBinaryUnstrechedRefIMG = imfill(IBinaryUnstrechedRefIMG,'holes');
        IBinaryUnstrechedRefIMG = imcomplement(IBinaryUnstrechedRefIMG);
        %closing
        se = strel('disk',3);
        IBinaryUnstrechedRefIMG = imclose(IBinaryUnstrechedRefIMG, se);
        %edges
        IEdgeUnstrechedRefIMG = edge(IBinaryUnstrechedRefIMG, 'canny');
        
        
        
        RefCoinAttributes = regionprops(IEdgeUnstrechedRefIMG,'centroid', 'area', 'BoundingBox');
        % (1) should work here if CANNY and stuff worked correclty
        diameter_ref_coin = sqrt(RefCoinAttributes(1).Area / pi) *2 ;
        %disp(diameter_ref_coin)
        %disp(j)
        
        strech_factor = diameter/diameter_ref_coin;
        stretchedRefIMG = imresize(unstretchedRefIMG,strech_factor);
        %
        % do not uncomment this - explodes
        %figure; imshow(stretchedRefIMG); title('stretched Ref IMG');
        %figure; imshow(unstretchedRefIMG); title('UN stretched Ref IMG');
        %figure; imshow(IMGCoin2Search); title('Coin2search');
        
        
        % d) make loop for 360 degrees
        %debug is from -5 : 5 degree
        %for k = 0 : 359
        for k = 0 : 2
            % e) rotate adjusted ref_coin(j) for this degree
            strechedRefIMGRotated = imrotate(stretchedRefIMG,k);
            
            %dont...
            %figure; imshow(strechedRefIMGRotated); title('streched and rotated');
            %figure; imshow(DIMGCoin2Search);
            
            % f) convert this rotated ref_coin(j) to double
            DstrechedRefIMGRotated = im2double(strechedRefIMGRotated);
            
            %hitMatrix = normxcorr2(DstrechedRefIMGRotated,IMGCoin2Search) > 0.7;
            hitMatrix = normxcorr2(DstrechedRefIMGRotated,IGray) > 0.425;
            
            %draw
            figure; subplot(1,2,1), imshow(IGray), title('Hits'); hold on;
                    
            ps = size(DstrechedRefIMGRotated,1)/2;
            for a=1:size(hitMatrix,1)
               for b=1:size(hitMatrix,2)
                   if hitMatrix(a,b) == 1
                        plot(b-ps,a-ps,'x','LineWidth',3,'Color','red');
                        %save data
                        foundData{iFoundCoins} = {ref_coin_name,a,b};
                        iFoundCoins = iFoundCoins+1;
                   end
               end
            end
            subplot(1,2,2), imshow(DstrechedRefIMGRotated), title('rotated Pattern'); hold on;
            hold off;
            
            

            
            disp(i)
            disp(j)
            disp(k)
            
            %disp(i); disp(j);disp(k);disp(xxx);
%             found = 0;
%             if sum(sum(hitMatrix)) > 1
%                 figure; imshow(DstrechedRefIMGRotated); title(sprintf('i = %d , j = %d , k = %d ',i,j,k));
%                 disp(ref_coin_name)
%                 found = 1;
%             end
            
            
            % try match both matrixes ( DstrechedRefIMGRotated) and  DIMGCoin2Search
            %diff{i,j,k} = sum(setdiff(DstrechedRefIMGRotated,DIMGCoin2Search));
            
            
            
            clear strechedRefIMGRotated;
            clear DstrechedRefIMGRotated;
            
            
        end
        
        

        
       
        clear stretchedRefIMG;
        clear unstretchedRefIMG;
        clear IBinaryUnstrechedRefIMG;
        clear IEdgeUnstrechedRefIMG;
        
    end 
    
    clear IMGCoin2Search;
    clear DIMGCoin2Search;
    %clear diameter;
    
   
    
end




toc
disp(money)
