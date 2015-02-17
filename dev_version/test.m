
mySettings = coin_recognition_settings();
mySettings();


%some tests
% pattern matching

% 1. get all images

fnames = dir('ref_coins/*.jpg');  % the folder in which ur images exists
for i = 1 : length(fnames)
    filename = strcat('ref_coins/',fnames(i).name);
    I = im2double(imread(filename));
    %coins{i} = {filename,I};
    coins{i} = I;
    %figure, imshow(I);
end

coins(2);
%imrotate(coins(2:2),1);


%actually, not I need to check the sizes of the coins in the picture
%so that if I know, if I have to shrink my templates
% so here is what we do: 
% 1. get all circles (diameters) from the image
% 2. adjust my template to it, and try all 360


X = imread('img_coins.jpg');
I = im2double(X);
IGray = rgb2gray(I);

stats = regionprops(X,'Area','BoundingBox');

iCoins = size(stats, 1);

for k=1:iCoins
    a = stats(k).Area; %print me
    
    bb = stats(k).BoundingBox;  
	
    subImage = imcrop(X, bb);
    figure(subImage);
    
end



%nCoins = i;
nCoins = 2;

for i = 1: 360
   for j = 1 : nCoins
       if mod(i,36) == 0 % nur jeder 36te, testweise 
        coins_rotated{j,i} =  imrotate(coins{j},i);    
       end      
   end
end


