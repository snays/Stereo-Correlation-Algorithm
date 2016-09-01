clear all;
%Read in the images 
disp('reading in images');
LI=imread('im6.png');

LI=LI(:,:,1);
[width,height]=size(LI); %image width n height
LI=im2double(LI);
LI(:,1:(height-width))=[];
%LI(1025:width,:)=[];

%LI=LI(1:2:end,1:2:end);
%LI=downsample(downsample(LI,4)',4);
%LI=LI';

%we WANT the abs value to be less than 20
RI=imread('im2.png');
RI=RI(:,:,1);
RI=im2double(RI);
RI(:,1:(height-width))=[];
%RI(1025:width,:)=[];


%RI=RI(1:2:end,1:2:end);
%RI=downsample(downsample(RI,4)',4);
%RI=RI';
%{
figure 
imshow(LI);
figure
imshow(RI);
%}
%% set parameters
hKern=9; %hKern width
maxDisp=120; %maximum disparity
minDisp=-5; %minimum disparity
[width,height]=size(LI); %image width n height
LR_thresh=10; %threshold for R->L check

xStart=hKern/2+0.5; %first valid xPosition (based on disparity of 0)
yStart=hKern/2+0.5; %first valid yPosition
xEnd=width-(hKern/2-0.5); %last valid xPosition 
yEnd=height-(hKern/2-0.5); %last valid yPosition
xMid1=xStart+maxDisp; 
xMid2=xEnd+minDisp;

count=0; %to count the number of pixels removed from the L<->R Consistency Check
minCount=0; %to count the number of pixels removed from minimum threshold
sum=0;
%% instantiate the buffers 
disp_map=ones([width,height])*3200; %disparity results
sub_disp_map=ones([width,height])*3200; %disparity results
cSum=zeros([(maxDisp-minDisp+1),width]); %column sum buffer
diff=zeros([(maxDisp-minDisp+1),width,hKern]); %3D matrix that stores absolute diff for all disparities, x values, and the hKern number of y values
SADbuffer=zeros([width,(maxDisp-minDisp+1)]); %Holds sum of absolute diff for the current x-d slice
next=32000; %this is equivalent to MISSING_PIXEL
tempResult=zeros([width,1]);

%initialize the first column sums (corresponding to the first row), for
%all disparities 
for d=minDisp:maxDisp
    if(d<0)
        xs=-d+1;
        xe=width;
    else
        xs=1;
        xe=width-d;
    end
    %for each (valid based on disparity) x position 
    for x=xs:xe
        %sum up the absolute difference
        for y=1:hKern
            %calculate and store absolute difference 
            diff(d-minDisp+1,x,y)=abs(LI(y,x)-RI(y,x+d));
            %add it to the column sum 
            cSum(d-minDisp+1,x)=cSum(d-minDisp+1,x)+diff(d-minDisp+1,x,y);
        end 
    end
   fprintf('Calculating column sums for disparity: %d\r',d); 
end

%find the SAD for the first point for each dimension
for d=minDisp:maxDisp
    if d<0
        %when d<0, the SAD values can only begin to be calculated at position xStart-d  
        for x=1-d:hKern-d
            SADbuffer(xStart-d,d-minDisp+1)=SADbuffer(xStart-d,d-minDisp+1)+cSum(d-minDisp+1,x);
        end
    else
        %otherwise, you can calculate the value of the SADbuffer starting from xStart 
        for x=1:hKern
            SADbuffer(xStart,d-minDisp+1)=SADbuffer(xStart,d-minDisp+1)+cSum(d-minDisp+1,x);
        end
    end
end


%Do the LEFT->RIGHT search (through the valid disparities) for every point on this row
for x=xStart+1:xEnd
        %reset the minimum values
        minval=3200000;
        minval2=3200000;
        minval3=3200000;
        %determine the valid disparity bounds 
        if(-minDisp-x>=minDisp)
            ds=-x-minDisp+1;
        else
            ds=minDisp;
        end
        if (xEnd-x>maxDisp)
            de=maxDisp;
        else
            de=xEnd-x;
        end
        %iterate through valid disparities
        for d=ds:de
            %calculate the SAD value(subtracting and adding column sums)
            SADbuffer(x,d-minDisp+1)=SADbuffer(x-1,d-minDisp+1)+cSum(d-minDisp+1,x+(hKern/2-0.5))-cSum(d-minDisp+1,x-(hKern/2+0.5));
            %update the minimum SAD value for this x position
            if(SADbuffer(x,d-minDisp+1)<minval)
                disp_map(y,x)=d;
                minval3=minval2;
                minval2=minval;
                minval=SADbuffer(x,d-minDisp+1);
            end
        end
        
        %check validitiy of minimum (make sure the minimum is sharp)
        if(minval3-minval<0.05*minval)
            disp_map(yStart,x)=3200;
        end
        %do the subpixel disparity if the disparity is not at the bounds
        if(disp_map(yStart,x)~=maxDisp && disp_map(yStart,x)~=minDisp && disp_map(yStart,x)~=3200)
            sub_disp_map(yStart,x)=subPixelApprox(disp_map(yStart,x),SADbuffer(x,disp_map(yStart,x)-minDisp+1-1),SADbuffer(x,disp_map(yStart,x)-minDisp+1),SADbuffer(x,disp_map(yStart,x)-minDisp+2));
        end    
    end

%do the RIGHT->LEFT search segmenting the x coordinates into three regions
%of valid disparity searching 

%REGION 1
for x=xStart:xMid1
    next=320000;
    for i=minDisp:x-xStart 
        %fprintf('Searching x:%d by looking at SAD(%d,%d)\r',x,x-i,i);
        if(SADbuffer(x-i,i-minDisp+1)<next)
            tempResult(x)=i;
            next=SADbuffer(x-i,i-minDisp+1);
        end
    end
end

%REGION 2
for x=xMid1+1:xMid2
    next=320000;
    for i=minDisp:maxDisp
           %fprintf('Searching x:%d by looking at SAD(%d,%d)\r',x,x-i,i);
            if(SADbuffer(x-i,i-minDisp+1)<next)
                tempResult(x)=i;
                next=SADbuffer(x-i,i-minDisp+1);
            end
    end
end

%REGION 3
for x=xMid2+1:xEnd
    next=320000;
    for i=minDisp+x-xMid2:maxDisp
        %fprintf('Searching x:%d by looking at SAD(%d,%d)\r',x,x-i,i);
        if (SADbuffer(x-i,i-minDisp+1)<next)
            tempResult(x)=i;
            next=SADbuffer(x-i,i-minDisp+1);
        end
    end
end

%CROSS CORRELATION : Compare and reset disparities (from the L->R and R->L) that are not within a threshold 
for x=xStart:xEnd
    if(disp_map(yStart,x)~=3200)
        if(abs(disp_map(yStart,x)-tempResult(x+disp_map(yStart,x))))>LR_thresh
            disp_map(yStart,x)=320000;
            count=count+1;
        end
    end
end

%% COMPUTE THE REST OF THE IMAGE 

for y=(yStart+1):yEnd
    fprintf('Computing SAD for row: %d\r',y);
    %reset the SADbuffer because we have now moved to a new y Slice
    SADbuffer=SADbuffer*0;

    %add a new slice to the bottom of the 3D difference matrix
    diff=cat(3,diff,zeros([(maxDisp-minDisp+1),width,1]));
    
    %calculate the difference buffer components
    for x=1:width
        %determine the bounds
        if(x+minDisp<=0)
            dStart=1-x;
        else
            dStart=minDisp;
        end   
        
        if(x+maxDisp<=width)
            dEnd=maxDisp;
        else
            dEnd=width-x;
        end
        %iterate through valid disparities
        for d=dStart:dEnd
                %compute diff for the new y-slice
                diff(d-minDisp+1,x,hKern+1)=abs(LI(y+hKern/2-0.5,x)-RI(y+hKern/2-0.5,x+d));
                %update the cSum by adding this new y-slice and subtracting the old
                cSum(d-minDisp+1,x)=cSum(d-minDisp+1,x)-diff(d-minDisp+1,x,1)+diff(d-minDisp+1,x,hKern+1);
        end
    end
    
    %free the top y-slice as we do not need it for any further computations
    diff(:,:,1)=[];

    
    %now that we have all the column sums..

    %find the SAD for the first point (for each disparity)
    for d=minDisp:maxDisp
        if d<0
            %when d<0, the SAD values can only begin to be calculated at position xStart-d
            for x=1-d:hKern-d
                SADbuffer(xStart-d,d-minDisp+1)=SADbuffer(xStart-d,d-minDisp+1)+cSum(d-minDisp+1,x);
            end
        else
    %otherwise, calculate SAD values for xStart 
            for x=1:hKern
                SADbuffer(xStart,d-minDisp+1)=SADbuffer(xStart,d-minDisp+1)+cSum(d-minDisp+1,x);
            end
        end
    end


    %Do the LEFT->RIGHT search (through the valid disparities) for every point on this row 
    for x=(xStart+1):xEnd
        %reset the minimum values 
        minval=3200000;
        minval2=3200000;
        minval3=3200000;
        %determine the valid disparity bounds
        if(-minDisp-x>=minDisp)
        ds=-x-minDisp+1;
        else
            ds=minDisp;
        end
        if (xEnd-x>maxDisp)
            de=maxDisp;
        else
            de=xEnd-x;
        end

        %fprintf('the bounds for x:%d %d %d\n',x,ds,de);
        %iterate through the valid disparities to search for minimum SAD
        for d=ds:de
            %calculate the next SAD value based on past value and column sums  
            SADbuffer(x,d-minDisp+1)=SADbuffer(x-1,d-minDisp+1)+cSum(d-minDisp+1,x+(hKern/2-0.5))-cSum(d-minDisp+1,x-(hKern/2+0.5));
            %update the minimum SAD   
            if(SADbuffer(x,d-minDisp+1)<minval)
                disp_map(y,x)=d;
                minval3=minval2;
                minval2=minval;
                minval=SADbuffer(x,d-minDisp+1);
            end               
        end
        %check validitiy of minimum (make sure minimum is sharp)
        if(minval3-minval<0.05*minval)
            disp_map(y,x)=3200;
        end
        %do the subpixel disparity if the disparity is not at the bounds
        if(disp_map(y,x)~=maxDisp && disp_map(y,x)~=minDisp && disp_map(y,x)~=3200)
            sub_disp_map(y,x)=subPixelApprox(disp_map(y,x),SADbuffer(x,disp_map(y,x)-minDisp+1-1),SADbuffer(x,disp_map(y,x)-minDisp+1),SADbuffer(x,disp_map(y,x)-minDisp+2));
        end    
    end
    
    
    %do the RIGHT->LEFT search (through the valid disparities)
    %by segmenting into the three valid disparity regions
    %REGION 1
    for x=xStart:xMid1
        next=320000;
        for i=minDisp:x-xStart 
            %fprintf('Searching x:%d by looking at SAD(%d,%d)\r',x,x-i,i);
            if(SADbuffer(x-i,i-minDisp+1)<next)
                tempResult(x)=i;
                next=SADbuffer(x-i,i-minDisp+1);
            end
        end
    end
    %REGION 2
    for x=xMid1+1:xMid2
        next=320000;
        for i=minDisp:maxDisp
            %fprintf('Searching x:%d by looking at SAD(%d,%d)\r',x,x-i,i);
            if(SADbuffer(x-i,i-minDisp+1)<next)
                tempResult(x)=i;
                next=SADbuffer(x-i,i-minDisp+1);
            end
        end
    end
    %REGION 3
    for x=xMid2+1:xEnd
        next=320000;
        for i=minDisp+x-xMid2:maxDisp
            %fprintf('Searching x:%d by looking at SAD(%d,%d)\r',x,x-i,i);
            if (SADbuffer(x-i,i-minDisp+1)<next && SADbuffer(x-i,i-minDisp+1)~=0 )
                tempResult(x)=i;
                next=SADbuffer(x-i,i-minDisp+1);
            end
        end
    end

    %CROSS CORRELATION : Compare and reset disparities (from the L->R and R->L) that are not within a threshold 
    for x=xStart:xEnd
        if(disp_map(y,x)~=3200)
            if(abs(disp_map(y,x)-tempResult(x+disp_map(y,x))))>LR_thresh
                disp_map(y,x)=320000;
                count=count+1;
            end
        end
    end
end

fprintf('There were %d / %d pixels removed from the R<->L Consistency Check', count, height*width);

%display results
figure
imshow(sub_disp_map/maxDisp)
figure
imshow(disp_map/maxDisp)