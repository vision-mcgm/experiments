function [  ] = LUTtest()
%Parameters
%FSN Jan 2013
%try
global g;



g.sim=0; %Simulate the run, not loading textures
g.fr=50; %Global maximum frame rate
g.int=1/g.fr;
g.pause=1;
g.smallPause=0.5;
g.videoFolder='frames_2\';
g.debugMode=1;
g.monitor=1; %1 for monitor with start bar
g.offset=100;
g.frames=200; %Number of frames to load in total
g.cropRect=[291   38  640  563]; %The rectangle of interest loaded from the images - x y length height
g.areaLength=10*50; %secs: areas we take both samples from
g.textColour=[0 0 0];

%We also have access to:
% g.centreX,centreY,black

%CONSTRAINT: total video length must be greater than three times
%sample length


%InitFramework must be called first as it sets groundwork variables

%InitialiseFramework();
%InitialiseExperiment();
%playClip(1,100,0,0,1,1,1);




screens=Screen('Screens');
sn=2;
[w1,windowRect]=Screen(sn,'OpenWindow',[127.5 127.5 127.5],[20 20 500 500],[],2);
[w2,windowRect]=Screen(sn,'OpenWindow',[127.5 127.5 127.5],[500 20 1000 500],[],2);


epar.GAMMA_TABLE = 'lut/nec-juin-2011';
epar.doGamma = 1;
[newGamma oldGamma]=exp_mon_init(w2,epar)




sca
ShowCursor;


%fclose all;

end

%END EXPERIMENT CODE

%----------------------------------------------------BEGIN PROJECT SPECIFIC HELPER FUNCTIONS



function []=InitialiseExperiment()
global g;

%We need to know the video dimensions

list=dir([g.videoFolder '*.bmp']);
alteredFolder=[g.videoFolder '\altered\'];
alteredList=dir([alteredFolder '*.bmp']);
img=imread([g.videoFolder list(1).name]);
vidHeight = size(img,1);
vidWidth = size(img,2);
g.h=g.cropRect(4);
g.w=g.cropRect(3);



% Preallocate movie structure.
mov(1:g.frames) = ...
    struct('cdata', zeros(g.h, g.w, 3, 'uint8'),...
    'colormap', []);

if ~g.sim
    for k = 1 :g.frames
        k
        img=imread([g.videoFolder list(k).name]);
        img=img(g.cropRect(2):g.cropRect(2)+g.cropRect(4),...
            g.cropRect(1):g.cropRect(1)+g.cropRect(3),:);
        %  altImg=imread([alteredFolder list(k).name]);
        %  altImg=altImg(g.cropRect(2):g.cropRect(2)+g.cropRect(4),...
        % g.cropRect(1):g.cropRect(1)+g.cropRect(3),:);
        
        % mov(k).cdata = read(xyloObj ,mod(k,nFrames)+1);
        
        g.texes(k)=Screen('MakeTexture',g.window,img);
        % g.texes(g.frames+k)=Screen('MakeTexture',g.window,altImg);
        
        % tex=g.texes(k);
        
    end
end
%This assumes that frames are numbered strictly alphabetically (leading
%zeros)
end

function []=firePractice()
global g;   
    reps=8;
    tLs=[0.2 0.5 1]*50;
    tLs=repmat(tLs,1,reps);
    n=reps*3;
    muls=[1.2 1.4 1.6 1.8 2];  
    for i=1:n        
     %  [num2str(g.targetLength) ' ' num2str(g.sampleLength) ]       
       correct=fireTrial(tLs(i),muls(randi(5))*tLs(i));
       corrects(i)=correct;
       %respLog([num2str(correct) ' ' num2str(g.targetLength) ' ' num2str(g.sampleLength)]);
       if correct, answer='Correct. '; else answer='Incorrect. '; end
       fireTextWait([ answer 'Accuracy: ' num2str(sum(corrects)/i)]);
       %fireText([ answer ]);      
    end
end


function [startA,startB]=pickSeparateSamples(areaLength,sampleLength)
global g;
%Returns two nonoverlapping sample offsets
%ASSUMPTION: areaL > 3* sampleL
%Ensure there is space left AFTER the first clip start, to play it
startA=randi(areaLength-sampleLength);
if startA>sampleLength && startA+sampleLength <areaLength-sampleLength
    %If there is space to pick clip B both before and after clip A
    choice=randi(2);
elseif startA<=sampleLength
    %If there is only space after
    choice=2;
else
    %Otherwise, there is only space before (provided that entire
    %video
    %length is greater than three times clip length)
    choice=1;
end

%Pick starting point of clip B
if choice ==1
    startB=randi(startA-sampleLength); %So that clip has time to play without overlapping
else
    Boffset=randi((areaLength-sampleLength)-(startA+sampleLength));
    %Boffset is the distance between the end of A and the start of B
    startB=startA+sampleLength+Boffset;
end

if startA > areaLength  || startB > areaLength
    ShowCursor
    keyboard
end
end



function [correct]=fireTrial(targetLength,sampleLength)
global g;

%Pick target area
areaStart=randi(g.frames-g.areaLength);

%Pick sample offsets
[offsetA offsetB]=pickSeparateSamples(g.areaLength,sampleLength);
startA=offsetA+areaStart; startB=offsetB+areaStart;

%Pick target source
targetSource=randi(2);

%Pick target offset within sample
targetOffset=randi(sampleLength-targetLength);
if targetSource==1, startTarget=targetOffset+startA;
else startTarget=targetOffset+startB; end

%Decide whether to swap order
swapOrder=randi(2)-1;
if swapOrder==0
    firstStart=startA;
    secondStart=startB;
else
    firstStart=startB;
    secondStart=startA;
end

%Play clips
inv=0;
neg=1;

pause(g.smallPause);
fixationSpot([255 0 0]);
fKbWait;
flip;
pause(g.smallPause);
playClip(startTarget,targetLength,inv,0,1,neg,3);
pause(g.smallPause);
fixationSpot([255 140 0]);
fKbWait;
flip;
pause(g.smallPause);
playClip(firstStart,sampleLength,inv,0,1,neg,3);
pause(g.smallPause);
fixationSpot([255 140 0]);
fKbWait;
flip;
pause(g.smallPause);
playClip(secondStart,sampleLength,inv,0,1,neg,3);
pause(g.smallPause);
fixationSpot([0 255 0]);

%Get response
if g.sim, r=randi(2)-1;
else r=getLR(); end
respClipNumber=xor(r,swapOrder)+1;
if respClipNumber==targetSource, correct=1; else correct=0; end
correct
respLog([num2str(targetSource) ' ' num2str(swapOrder) ' ' num2str(correct)...
    ' ' num2str(startA) ' ' num2str(startB) ' ' num2str(targetOffset)]);
end

function [  ] = playClip( start,length,negative,backwards,interval,regime,location )
global g;
%Plays a constant clip starting at start with a certain length
%Regimes:
% 1 Forwards
% 2 Negative
clear timestamps
fList=0:length-1;
fList=interval.*floor(fList/interval);

if backwards
    fList=-fList;
    fList=fList+(start+length-1);
else
    fList=fList+start;
end;

if negative, fList=fList+g.frames;end;
fList;

if regime==1
    angle=0;
elseif regime==2
    angle=180;
end

if location==1
    rect=[g.centreX-(g.w/2)-g.offset g.centreY-(g.h/2) g.centreX+(g.w/2)-g.offset g.centreY+(g.h/2)];
elseif location==2
    rect=[g.centreX-(g.w/2)+g.offset g.centreY-(g.h/2) g.centreX+(g.w/2)+g.offset g.centreY+(g.h/2)];
elseif location ==3
    rect=[g.centreX-(g.w/2) g.centreY-(g.h/2) g.centreX+(g.w/2) g.centreY+(g.h/2)];   
end

nFToPlay=size(fList,2);
intervals=0:nFToPlay;
intervals=intervals*g.int;
%Remember GetSecs, timestamps etc, all return the time in seconds since
%system startup. This is a high precision number.

if ~g.sim
    it=1; %Should always place the iterator declaration right before the loop starts
    for k=fList
        
        Screen('DrawTexture', g.window, g.texes(k),[],rect,angle);
        
        if 0 %Old method with incremental nextFlip calc
            if it==1
                
                % Update display:
                [VBLTimestamp StimulusOnsetTime FlipTimestamp Missed Beampos] = Screen('Flip', g.window);
                timestamps(it)=VBLTimestamp;
                nextStamp=VBLTimestamp+ g.int;
            else
                
                % [VBLTimestamp StimulusOnsetTime FlipTimestamp Missed Beampos] =
                % Screen('Flip', g.window,nextStamp); %Original
                %Flip immediately:
                [VBLTimestamp StimulusOnsetTime FlipTimestamp Missed Beampos] = Screen('Flip', g.window);
                % Release texture:Screen('Flip', window);
                timestamps(it)=VBLTimestamp;
                nextStamp=VBLTimestamp+ g.int;
            end
        else
            %New method with entire framerate calc at beginning
            if it==1
                
                % Update display:
                [VBLTimestamp StimulusOnsetTime FlipTimestamp Missed Beampos] = Screen('Flip', g.window);
                timestamps(it)=VBLTimestamp;
                nextStamps=VBLTimestamp + intervals;
            else
                
                % [VBLTimestamp StimulusOnsetTime FlipTimestamp Missed Beampos] =
                % Screen('Flip', g.window,nextStamp); %Original
                %Flip immediately:
                [VBLTimestamp StimulusOnsetTime FlipTimestamp Missed Beampos] = Screen('Flip', g.window,nextStamps(it));
                % Release texture:Screen('Flip', window);
                timestamps(it)=VBLTimestamp;
                nextStamp=VBLTimestamp+ g.int;
            end
        end
        it=it+1;
    end
    Screen('FillRect',g.window,g.grey, rect);
    Screen(g.window, 'Flip');
    
    for i=1:size(timestamps,2)-1
        diffs(i)=timestamps(i+1)-timestamps(i);
    end
   % 'Framerate'
   % 1/mean(diffs)
   % std(diffs)
   fprintf('Played %d %d\n',length,it-1);
end
end


%------------------------------------BEGIN COMMON EXPERIMENTAL HELPER FUNCTIONS

function []=flip()
global g;
Screen('Flip',g.window);
end

function[]=fKbWait()
global g;
if ~g.sim, KbWait; end
end

function [] = InitialiseFramework()
%Sets up PsychToolbox and other parts of the experimental framework
global g;

if g.sim, g.pause=0; end
if ~g.debugMode
    g.subjectName=input('Enter subject initials: ','s');
else
    g.subjectName='test';
end

g.logFileID=fopen([ 'output\' g.subjectName '_log.txt'],'w');
assert(g.logFileID ~= -1); %Check file is really open
fireLog(['Starting experiment on subject ' g.subjectName]);
g.responseFileID=fopen([ 'output\' g.subjectName '_results.txt'],'w');


if ~g.debugMode;  HideCursor(); end

Screen('Preference', 'SkipSyncTests', 1)
screens=Screen('Screens');
screenNumber=max(screens);


% Open a window.  Note the new argument to OpenWindow with value 2,
% specifying the number of buffers to the onscreen window.
if ~g.sim
    if g.debugMode
        %[g.window,windowRect]=Screen(screenNumber,'OpenWindow',0,[20 20 800 800],[],2);
        [g.window,windowRect]=Screen(screenNumber,'OpenWindow',[127.5 127.5 127.5],[20 20 700 700],[],2);
    else
        [g.window,windowRect]=Screen(screenNumber,'OpenWindow',[127.5 127.5 127.5],[],[],2);
    end
    
    %pixelSize=Screen('PixelSize', windowPtrOrScreenNumber);
    g.centreX=(windowRect(3)+windowRect(1))/2;
    g.centreY=(windowRect(4)+windowRect(2))/2;
    g.black=BlackIndex(g.window);
    white = WhiteIndex(g.window);  % Retrieves the CLUT color code for white.
    g.grey = (g.black + white) / 2;
    Screen('FillRect',g.window,g.grey);
end
end

function [y,i]=pick(x)

%Returns value and index randomly sampled
l=length(x);
if l
    i=randi(l);
    y=x(i);
    
else i=0;
    idx=0;
end
end

function [tf]=ItemsRemaining(list)
if find(list),tf=1, else tf=0;,end
end

function [index,list]=GetIndexAndDecrement(list)
%Returns the index to one of the nonzero slots in a list,
%and decrements that slot.
picked=0;
nonEmptySlots=find(list);
[index]=pick(nonEmptySlots);
list(index)=list(index)-1;
end

function fireLog(s)
%Log a string to the logfile for this experiment
global g;
logString=[datestr(now) ': ' s];
fprintf(g.logFileID,[logString '\n']);
fprintf([logString '\n']);

end

function respLog(s)
%Log a string to the logfile for this experiment
global g;

fprintf(g.responseFileID,[s '\n']);
end

function fireText(s)
global g;
if ~g.sim
    DrawFormattedText(g.window,s,'center','center',g.textColour)
    Screen('Flip',g.window);
    fKbWait;
end
end

function fireTextWait(s)
global g;
if ~g.sim
    DrawFormattedText(g.window,s,'center','center',g.textColour)
    Screen('Flip',g.window);
    pause(1);
end
end

function []=fixationSpot(varargin)
global g;
if nargin==0, colour=[0 0 0];
else colour=varargin{1}; end

if ~g.sim
    % Draw a fixation point
    Screen('DrawDots', g.window, [g.centreX g.centreY], 12, [255 255 255], [], 2);
    Screen('DrawDots', g.window, [g.centreX g.centreY], 8, colour, [], 2);
    % Show it at next retrace
    Screen('Flip',g.window);
end
end

function []=checkQuit()
[ keyIsDown, kt, keyCode ] = KbCheck;
response_key = KbName(keyCode);
%             a = strcmp(response_key, 'RightArrow');
%             b = strcmp(response_key, 'LeftArrow');
if strcmp(response_key, 'q')
    sca
end
end


function [r]=exitExperiment()
sca
ShowCursor;
keyboard
end

function [r]=getLR()
%Gets key response - 0 for left or 1 for right.
global g;
a = 0;
b=0;
while a == 0 && b ==0
    %WaitSecs(0.2);
    KbWait;
    [ keyIsDown, kt, keyCode ] = KbCheck;
    %keyCode is a 256-element list saying which keys are down, thus we must
    %not truncate it
    
    response_key = KbName(keyCode);
    %             a = strcmp(response_key, 'RightArrow');
    %             b = strcmp(response_key, 'LeftArrow');
    %If several keys are down, select only one
    if strcmp(class(response_key),'cell'), response_key=response_key{1}; end
    
    if 1
        if strcmp(response_key,'q')
            exitExperiment();
        end
    end
    a = strcmp(response_key, 'right');
    if a(1)==1; a=1; else a=0; end
    %Covers the case where a is a vector
    
    
    %Makes sure a is always 1 or 0 and not anything else
    b = strcmp(response_key, 'left');
    if b(1)==1; b=1; else b=0; end
    
    %
    if g.debugMode
        if strcmp(response_key, 'k')
            keyboard
        end
    end
end

if b
    r=0;
else
    r=1;
end
end
