function [twoDHeader twoDResults] = Analysis_DICOM_Input_Sub_Auto_Thresh(startSlice,stopSlice,numbones,paths,bwa,outputfile,boneNumber,allSlices)

% 2D Analysis Program
% Version 4.0
% Ryan Tomlinson
% 8/13/2008
% upd 7/13/2009

% INPUT: 2D CT images (preferably .tif format)
% OUTPUT: XLS Spreadsheet (tab-delimited) of bone measurements

% CHANGE LOG
% 8/29/2008
% Added try/catch for filename issues
% Added error dialog box if output does not exist

% 9/12/2008
% Adding Tissue Area (T.Ar, area inside the periosteal perimeter)
% Medullary Area (T.Ar - B.Ar)
% Periosteal X-Width (x-max - x-min)
% Periosteal Y-Width (y-max - y min)
% as per email 9/11 MJS

% 9/18/2008
% Changed from rotate_image to imrotate. Much less noise in rotation.

% 5/28/2009 (Mark Willingham)
%  Added Ixy 

% 6/25/2009 (Mark Willingham)
% Added Imin and Imax 

% 7/09/09 (Mark Willingham)
% Added automatic rotation by minimizing product MOI
    % Prompt for auto-rotation before asking for angles
    % Setting to ensure bone lays on its long side
% Fixed issues with "find" that switched Ixx & Iyy
% Fixed error in "slicenum" indexing that led to incorrect results
% Added user-defined threshold minima to fix issues with some radius slices
% Added averages to output file after each bone

%7/13/09 (Mark Willingham)
% Clear Ixx_rt(:) etc. after every bone - correct results now
% Remove iterative angle finding now that moments are correct
% Store all images in 'a' to increase speed when processing multiple bones
% Output angle after each bone (no  longer automatically computed before headers are written)
% Improve rotation dialog box order - Auto? no longer in manual roation loop
% Change how phi calculates angle for more accurate rotations
% Removed user-defined threshold minima as it's not less necessary

% 10/02/09 Ryan Tomlinson
% Fixed bug in ymax calculation

% 12/30/13 Dan Leib
% Altered to accept DICOM input so metadata can be used instead of tube
% size inputs


numslices = abs(stopSlice - startSlice);
ext = 'dcm*';
thresh_cutoff = 0.4;

% PICK FIRST IMAGE
fileList = dir([paths '\*.dcm*']);
if strcmpi(allSlices,'y') == 1
    startSlice = 1;
    stopSlice = length(fileList);
end
file = fileList(1).name;
specname = paths;
files=dir([paths '\*.' ext]);
start = startSlice;

stop = start+numslices-1;
if stop > length(files)
    stop = length(files);
end

% FORMAT OUTPUT FILE
proceed = 0;
autorotate = 0;
while proceed == 0
    try
        output = fopen(outputfile,'wt');
        fprintf(output,[specname '\t']); 
        fprintf(output,'%s \n',outputfile);
        fclose(output);
        proceed = 1;
    catch
        for i=1:length(outputfile)
            if outputfile(i) == '\\'
                pos = i;
            end
        end
        [xx,xxx] = system(['mkdir "' outputfile(1:pos-1) '"']);
        if xx>0
            uiwait(errordlg('Output filepaths is invalid. Please restart.','Error'));
            proceed = 1;
        else
            uiwait(errordlg('Output filepaths did not exist. The directory has been created.','Error'));
        end
    end
end

for q=1:numbones
    % LOAD IMAGES
    img = dicomread([paths '\' file]);
    uniques = unique(img);
    if length(uniques) == 2
        thresh = 0.5;
    else
    end 
    info = dicominfo([paths '\' file]);
    thresh = graythresh(img);
    if thresh < thresh_cutoff
        img = img > thresh_cutoff*double(max(img(:)));
    else
        img = im2bw(img,graythresh(img));
    end
    
    pixelwidth = info.SliceThickness;
    pixelarea = pixelwidth^2;
    
    % NEGATIVE IMAGE
    if q==1
        whiteblack = 'White';%questdlg('Bone: White or Black?','Negative Function','White','Black','White');
        if whiteblack(1) == 'B'
            img = not(img);
        end
%         close all;
    end
    cc = bwconncomp(img);
    numPixels = cellfun(@numel,cc.PixelIdxList);
    [biggest,idx] = max(numPixels);
    img = false(size(img));
    img(cc.PixelIdxList{idx}) = true;
    proceed = 0;
        if q == 1
            quest = 'Yes, Horizontal Long Axis';%questdlg('Rotate All Bones Automatically?','Rotation Function','Yes, Horizontal Long Axis', 'Yes','No, Manually Set Angle','Yes Horizontal Long Axis');
            if quest(1) == 'Y'
                autorotate = 1;
                proceed = 1;
                degree(1:numbones) = 0;
                if length(quest) > 3
                    autohoriz = 1;
                else
                    autohoriz = 0;
                end
            end
        end
    while proceed == 0 && autorotate == 0
        quest = 'Yes';%questdlg('Finished rotating?','Rotation Function','Yes','No','Yes');
        if quest(1) == 'Y'
            proceed = 1;
        end
%         close all;
    end
end

%% NOW PROCESSING
tic;
%Read in images
for k = start:stop
	img = dicomread([paths '\' files(k).name]);
	if k == start
        a = uint16(img);
	else
        a = [a img];
	end
end
res = length(a(1,:)) / (stop-start+1);
% close('all');
for q=1:numbones
    % Calculate the correct rotation - computes phi to minimize the average Ixy
    if autorotate == 1
        slicenum = 1;
        for i = 1:res:length(a)
            img = a(:,i:i+res-1);
            thresh = graythresh(img);
            if thresh < thresh_cutoff
                img = img > thresh_cutoff*double(max(img(:)));
            else
                img = im2bw(img,graythresh(img));
            end
            if whiteblack(1) == 'B'
                img = not(img);
            end
            img = imrotate(img, degree(q));
            cc = bwconncomp(img);
            numPixels = cellfun(@numel,cc.PixelIdxList);
            [biggest,idx] = max(numPixels);
            img = false(size(img));
            img(cc.PixelIdxList{idx}) = true;
            % Find position of Centroid
            [I J] = find(img > 0);
            I = I .* pixelwidth; %mm (pixels * mm/pixels)
            J = J .* pixelwidth; %mm
            ycent = mean(I(:));
            xcent = mean(J(:));
        
            % Find MOIs
            for j=1:length(I)
                Iyy_rt(j) = double((J(j) - xcent)^2) * pixelarea;
                Ixx_rt(j) = double((I(j) - ycent)^2) * pixelarea;
                Ixy_rt(j) = (double((I(j) - ycent)*(J(j) - xcent)))*pixelarea;
            end
            
            Iyy(q,slicenum) = sum(Iyy_rt(:)); clear Iyy_rt;
            Ixx(q,slicenum) = sum(Ixx_rt(:)); clear Ixx_rt;
            Ixy(q,slicenum) = sum(Ixy_rt(:)); clear Ixy_rt;

            slicenum = slicenum + 1;    
        end
        start = 1;
        stop = start+(stopSlice-startSlice)-1;
        phi(q)  = pi/180*degree(q) - 0.5*atan(2*mean(Ixy(q,start:stop))/(mean(Ixx(q,start:stop))-mean(Iyy(q,start:stop))));
        
        % Ensure the long axis is horizontal (if that option was previously selected)
        if mean(Iyy(q,start:stop)) < mean(Ixx(q,start:stop)) && autohoriz == 1
            degree(q) = 90+180/pi*phi(q);
        else
            degree(q) = 180/pi*phi(q);
        end
    end
%     close('all');
    
	%Analyze all the slices with the given rotations, croppings, etc.
    slicenum = 1;
    for i = 1:res:length(a)
        img = a(:,i:i+res-1);
        thresh = graythresh(img);
        if thresh < thresh_cutoff
            img = img > thresh_cutoff*double(max(img(:)));
        else
            img = im2bw(img,graythresh(img));
        end
%         img = imcrop(img,rect{q});
        if whiteblack(1) == 'B'
            img = not(img);
        end
        img = imrotate(img, degree(q));
        cc = bwconncomp(img);
        numPixels = cellfun(@numel,cc.PixelIdxList);
        [biggest,idx] = max(numPixels);
        img = false(size(img));
        img(cc.PixelIdxList{idx}) = true;

        
        % Find position of Centroid
        [I J] = find(img > 0);      
        I = I .* pixelwidth; %mm (pixels * mm/pixels)
        J = J .* pixelwidth; %mm
        ycent = mean(I(:));
        xcent = mean(J(:));
  
        % Find Periosteal X width and Y width
        pywidth(q,slicenum) = max(I) - min(I);
        pxwidth(q,slicenum) = max(J) - min(J);
        
        % Find Bone Area
        bonearea(q,slicenum) = length(I)*pixelarea;
        
        % Find MOIs
        for j=1:length(I)
            Iyy_rt(j) = double((J(j) - xcent)^2) * pixelarea;
            Ixx_rt(j) = double((I(j) - ycent)^2) * pixelarea;
            Ixy_rt(j) = (double((I(j) - ycent)*(J(j) - xcent)))*pixelarea;
        end
        
        Iyy(q,slicenum) = sum(Iyy_rt(:)); clear Iyy_rt;
        Ixx(q,slicenum) = sum(Ixx_rt(:)); clear Ixx_rt;
        Ixy(q,slicenum) = sum(Ixy_rt(:)); clear Ixy_rt;
        
        %Debugging code to record & plot the images
        if slicenum == start 
            b = img;
		elseif slicenum == stop
            b = [b img];
%             subplot(numbones,1,q); imshow(b);
        else
            b = [b img];
        end

        % Find Ymax (distance from centroid to bottom of the bone)
        ymax(q,slicenum) = abs(max(I - ycent));
        
        % Find Tissue Area (T.Ar)
        tissuearea = imfill(img,[round(ycent/pixelwidth) round(xcent/pixelwidth)]);
        tarea(q,slicenum) = sum(tissuearea(:))*pixelarea;
        
        % Find Medullary Area (M.Ar)
        medarea = not(imfill(img,[1 1]));
        marea(q,slicenum) = sum(medarea(:))*pixelarea;
        
        % Find average Cortical Thickness
        ycentpixel = round(ycent/pixelwidth);
        xcentpixel = round(xcent/pixelwidth);
        
        horpixels = (img(ycentpixel,:));
        verpixels = (img(:,xcentpixel));
        
        %12:00
        first = 0;
        for j=1:ycentpixel
            if sum(img(j, xcentpixel-1:xcentpixel+1)) > 0
                if first == 0
                    noonfirst = j;
                    first=1;
                else
                    noonlast = j;
                end
            end
        end
        noonth(q,slicenum) = abs(noonlast-noonfirst) * pixelwidth;
        
        %3:00
        first = 0;
        for j=xcentpixel:length(horpixels)
            if sum(img(ycentpixel-1:ycentpixel+1,j)) > 0
                if first == 0
                    threefirst = j;
                    first=1;
                else
                    threelast = j;
                end
            end
        end
        threeth(q,slicenum) = abs(threelast-threefirst) * pixelwidth;
        
        %6:00
        first = 0;
        for j=ycentpixel:length(verpixels)
            if sum(img(j,xcentpixel-1:xcentpixel+1)) > 0
                if first == 0
                    sixfirst = j;
                    first=1;
                else
                    sixlast = j;
                end
            end
        end
        sixth(q,slicenum) = abs(sixlast-sixfirst) * pixelwidth;
        
        %9:00
        first = 0;
        for j=1:xcentpixel
            if sum(img(ycentpixel-1:ycentpixel+1,j)) > 0
                if first == 0
                    ninefirst = j;
                    first=1;
                else
                    ninelast = j;
                end
            end
        end
        nineth(q,slicenum) = abs(ninelast-ninefirst) * pixelwidth;
        avgth(q,slicenum) = (noonth(q,slicenum)+threeth(q,slicenum)+sixth(q,slicenum)+nineth(q,slicenum))/4;
        
        %Output rotations and headers:	
        if q == 1 && slicenum == start
            output = fopen(outputfile,'at+');
			fprintf(output,'Bone Number\t'); 
            fprintf(output,'Start and Stop Slices\t');
			fprintf(output,'Bone Area (mm^2)\t'); 
            fprintf(output,'Iyy (mm^4)\t'); 
			fprintf(output,'Ixx (mm^4)\t'); 
            fprintf(output,'Ixy (mm^4)\t'); 
			fprintf(output,'Ymax (mm)\t'); 
			fprintf(output,'T.ar (mm^2)\t'); 
            fprintf(output,'M.ar (mm^2)\t'); 
			fprintf(output,'Periosteal X width (mm)\t'); 
            fprintf(output,'Periosteal Y width (mm)\t'); 
			fprintf(output,'12:00 Thickness (mm)\t');
            fprintf(output,'3:00 Thickness (mm)\t'); 
			fprintf(output,'6:00 Thickness (mm)\t'); 
            fprintf(output,'9:00 Thickness (mm)\t'); 
			fprintf(output,'Average Thickness (mm)\t'); 
            fprintf(output,'Bone Rotation (°) \n'); fclose(output);
            twoDHeader = {'Bone Number' 'Bone Area' 'Iyy (mm^4' 'Ixx (mm^4)' 'Ixy(mm^4)' 'Ymax' 'T.ar' 'M.ar' 'Periosteal X Width'...
                'Periosteal Y Width' '12:00 Thickness' '3:00 Thickness' '6:00 Thickness' '9:00 Thickness' 'Average Thickness' 'Bone Rotation'};
        end
        slicenum = slicenum + 1;
    end
	%Output the averages after each bone's data
	output = fopen(outputfile,'at+');
    fprintf(output,'%s\t',num2str(boneNumber));
    fprintf(output,'%s\t',[num2str(startSlice) '-' num2str(stopSlice)]);
	fprintf(output,'%10.5f \t',mean(bonearea(q, start:stop)));
	fprintf(output,'%10.5f \t',mean(Iyy(q, start:stop)));
	fprintf(output,'%10.5f \t',mean(Ixx(q, start:stop)));
	fprintf(output,'%10.5f \t',mean(Ixy(q, start:stop)));
	fprintf(output,'%10.5f \t',mean(ymax(q, start:stop)));
	fprintf(output,'%10.5f \t',mean(tarea(q, start:stop)));
	fprintf(output,'%10.5f \t',mean(marea(q, start:stop)));
	fprintf(output,'%10.5f \t',mean(pxwidth(q, start:stop)));
	fprintf(output,'%10.5f \t',mean(pywidth(q, start:stop)));
	fprintf(output,'%10.5f \t',mean(noonth(q, start:stop)));
	fprintf(output,'%10.5f \t',mean(threeth(q, start:stop)));
	fprintf(output,'%10.5f \t',mean(sixth(q, start:stop)));
	fprintf(output,'%10.5f \t',mean(nineth(q, start:stop)));
	fprintf(output,'%10.5f \t',mean(avgth(q, start:stop)));
    fprintf(output, '%5.2f \n',degree(q));
	fclose(output);
    
    
end

% t = toc;
% disp(['Procesing ' num2str(numbones) ' bones of ' specname ' took ' num2str(round(t*10)/10) ' seconds.']);
% disp(['Output available in ' outputfile]);