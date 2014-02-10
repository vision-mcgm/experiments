function [  ] = fire18Run()
%Backwards detection over sampling rates
%Parameters
clear all;
global g;

g.sim=0; %Simulate the run, not loading textures
g.respSim=1;
g.debugMode=0; %Debug mode
g.fr=50; %Global maximum frame rate
g.int=1/g.fr;
g.pause=1;
g.smallPause=0.5;
g.videoFolder='../../cropped_20-24mins/'
g.monitor=1; %1 for monitor with start bar
g.offset=200;
g.frames=10000; %Number of frames to load in total
g.cropRect=[291   38  640  563]; %The rectangle of interest loaded from the images - x y length height
g.textColour=[0 0 0];
g.clusterOutputFolder='W:\Fintan\Experiments\autoOutput\';
g.sampleAreaPause=1;
g.decSpec='%05d';
g.komodo=0;
%Factor setup
g.inversion=1:2;
%Normal,inverted
g.direction=1:2;
g.resize=0; %Switch on and off scaling to save time
g.scaling=0.7; %Image scaling
%forwards, backward
g.feedback=1;
g.missed=0;
%Experiment parameters
g.areaLength=0.5*50; %secs: areas we take both samples from
g.sampleLength=0.4*50;
g.exName='fire17';
m=machineParams;
g.maxFramesInVM=500;
g.intercache=0;


InitialiseFramework();
InitialiseExperiment();
ShowCursor;

%RunTraining();

ScheduleExperiment(); %Get the trial scheduling
RunExperiment();

keyboard
exitExperiment();
end

%END EXPERIMENT CODE

%----------------------------------------------------BEGIN PROJECT SPECIFIC HELPER FUNCTIONS

function RunTraining()
global g;
fireTextConfirm('Ready to start training.\nHit any key to start.');
nTraining=30;
nWindow=10;
trainingParamsCell{1}=1;
trainingParamsCell{2}=0;
trainingParamsCell{3}=0;
for i=1:nTraining
    trainingParamsArray=fireScheduleTrial(trainingParamsCell);
    CacheN(500);
    tps=fireTrial(trainingParamsArray);
    resps(i)=tps.correct;
    fireTextWait(['Accuracy: ' num2str(mean(resps(max(1,end-nWindow):end)))]);
end

%VITAL: finish training
g.feedback=0;
g.training=0;
fireTextConfirm('Training over. WAIT FOR THE EXPERIMENTER.');
waitForX;
end


function ScheduleExperiment
global g;
%Init
g.nextFrI=1;
g.FramesInVM=0;

g.nQ=1;

[allTrialParams blockParamCells g.info g.design]=fireReadParams();
g.allTrialParamsArray(1400)=fireScheduleTrial(allTrialParams(1,:));

for i=1:g.info.nTrials
    i
    g.allTrialParamsArray(i)=fireScheduleTrial(allTrialParams(i,:));
    if i==1
         g.allTrialParamsArray(g.info.nTrials)= g.allTrialParamsArray(1);
    end
end
end

function RunExperiment
global g;
makeSureAtLeastN(500);
trial=1;
g.feedback=0;
g.training=0;

for b=1:g.info.nBlocks
    fireTextConfirm(['Block ' num2str(b) ' of ' num2str(g.info.nBlocks) ' is about to start.\nHit any key to begin.']);
    for tInBlock=1:g.info.nTrialsPerBlock
        fprintf('%d %d %d\n');
        tstart=now;
        thisTrialParams=fireTrial(g.allTrialParamsArray(trial));
        thisTrialParams.block=b;
        thisTrialParams.tend=now;
        thisTrialParams.tstart=tstart;
        tps(trial)=thisTrialParams;
        respSave(tps);
        trial=trial+1
        makeSureAtLeastN(500);
    end
    fireTextConfirm(['You have finished block ' num2str(b) ' of ' num2str(g.info.nBlocks) '.\nPlease have a short break, then hit any key to continue.']);
end


end


function trialParams=fireScheduleTrial(paramCell) 
%Matlab treats a passed cell array as multiple functions
%Loads frames into the video queue for a particular trial
%Remember we assume each frame will only be used once
global g;
%Leftover vars
neg=0;
chrom=0;
location=3;
direction=0;
sampleRate=1;
angle=1;
Lsample=50;
%Extract params
%CAREFUL WIRING THIS UP
inv=paramCell{1};
Lpre=paramCell{2};
Lpost=paramCell{3};

Ltest=Lsample+Lpre+Lpost;

%Pick target area
[Stest Sfalse]=pickSeparateSamples(g.frames, Ltest);
source=randi(2);
%Pick start point within sample
sampleOffset=Lpre;
YN=randi(2);

switch YN
    case 1 %Yes
        Ssample=sampleOffset+Stest;
    case 2 %No
        Ssample=sampleOffset+Sfalse;
end

%Copy params to tp
if ~g.sim
    trialParams.startAQ=EnqueueFrames(Ssample,Lsample,neg,direction,sampleRate,angle,location,chrom);
    trialParams.startBQ=EnqueueFrames(Stest,Ltest,neg,direction,sampleRate,angle,location,chrom);
end
%xpt-specific
trialParams.Lpre=Lpre;
trialParams.Lpost=Lpost;
trialParams.YN=YN;
trialParams.Stest=Stest;
trialParams.Sfalse=Sfalse;
trialParams.Ssample=Ssample;
trialParams.source=source;
trialParams.Lsample=Lsample;
trialParams.Ltest=Ltest;
%Common
trialParams;
end

function [tp]=fireTrial(tp)
global g;
neg=0;
chrom=0;
location=3;
length=100;
direction=0;
sampleRate=1;
angle=1;
%Play sample
firePause(1);
fixationSpot([0 255 0]);

if ~g.sim
    pause(1);
    playClip(tp.startAQ,tp.Lsample,neg,direction,sampleRate,angle,location,chrom);
    fireClear;
    firePause(1);
    fixationSpot([0 255 0]);
    firePause(1);
    
    playClip(tp.startBQ,tp.Ltest,neg,direction,sampleRate,angle,location,chrom);
        fixationSpot([0 0 255]);
end

%Get response
if g.sim
    response=makeResponse(direction,sampleRate);
elseif g.respSim
    response=randi(2)-1;
else
    response=getYN();
end
flip;
%Work out answer
tp.correct=0;
if tp.YN==1 %yes
    if response==0 tp.correct=0; else tp.correct=1;end
elseif tp.YN==2 %no
    if response==0,tp.correct=1;else tp.correct=0;end
end

if g.sim
    
        if rand() <0.6
            tp.correct=1;else tp.correct=0;end
   
end

%Remember to use tp.correct not correct

if g.feedback
    if tp.correct, answer='Correct';, else answer='Incorrect';,end
    fireTextWait(answer);
end
end

function []=CheckAssertions()
global g;
%Assertions
assert(g.frames > 3*(max(g.sampleLengths) ));
end

function []=InitialiseExperiment()
global g;
%We need alpha blending for this experiment
% Enable alpha-blending, set it to a blend equation useable for linear
% superposition with alpha-weighted source. This allows to linearly
% superimpose gabor patches in the mathematically correct manner, should
% they overlap. Alpha-weighted source means: The 'globalAlpha' parameter in
% the 'DrawTextures' can be used to modulate the intensity of each pixel of
% the drawn patch before it is superimposed to the framebuffer image, ie.,
% it allows to specify a global per-patch contrast value:
if ~g.sim
    %Things we don't want to do if it's a simulation
    % Screen('BlendFunction', g.window, GL_ONE, GL_ONE); %We don't want to
    % do this at all!
end
%We need to know the video dimensions

list=dir([g.videoFolder 'normal/*.bmp']);
%alteredFolder=[g.videoFolder '\altered\'];
%alteredList=dir([alteredFolder '*.bmp']);
img=imread([g.videoFolder 'normal/' list(1).name]);
img=imresize(img,g.scaling);
vidHeight = size(img,1);
vidWidth = size(img,2);
g.h=vidHeight;
g.w=vidWidth;
g.prevTextures=[];

if 0
    for k = 1 :g.frames
        k;
        fireTextNoWait(num2str(k));
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

function []=firePractice(n)
global g;

for i=1:n
    %  [num2str(g.targetLength) ' ' num2str(g.sampleLength) ]
    correct=fireTrial(randi(2),randi(2),30);
    corrects(i)=correct;
    %respLog([num2str(correct) ' ' num2str(g.targetLength) ' ' num2str(g.sampleLength)]);
    if correct, answer='Correct. '; else answer='Incorrect. '; end
    fireTextWait([ answer 'Accuracy: ' num2str(sum(corrects)/i)]);
    %fireText([ answer ]);
end
clear corrects

end


function [startA,startB]=pickSeparateSamples(areaLength,sampleLength)
global g;
%Returns two nonoverlapping sample offsets
%ASSUMPTION: areaL > 3* sampleL
assert(areaLength > 3* sampleLength);
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



function startQ=EnqueueFrames(start,length,negative,backwards,interval,angle,location,chromatic)
global g;
if ~g.sim
    %Plays a constant clip starting at start with a certain length
    %Regime controls INVERSION. 1=up 2=inverted
    fprintf('Enqueuing clip...\n');
    Tstart=GetSecs();
    if g.sim
        return
    end
    
    clear timestamps
    fList=0:length-1;
    fList=interval.*floor(fList/interval);
    
    if backwards
        fList=-fList;
        fList=fList+(start+length-1);
    else
        fList=fList+start;
    end;
    fList_handles=fList;
    if negative, fList_handles=fList+g.frames;end;
    if chromatic, fList_handles=fList+2*g.frames;end;
    fList;
    
    %Now we have the final frame list, so can cache.
    
    nTexes=size(fList,2);
    
    if ~g.sim
        for k=1:nTexes
            if k==1
                 startQ=EnqueueFrame(fList(k));
            else
                EnqueueFrame(fList(k));
            end
            %         k;
            %         if negative
            %             imPath=[g.videoFolder 'negative\frame' num2str(fList(k),g.decSpec) '.bmp'];
            %         end
            %         if chromatic
            %             imPath=[g.videoFolder 'hue\frame' num2str(fList(k),g.decSpec) '.bmp'];
            %             im=imread(imPath);
            %             im=uint8(im);
            %             im=im*255;
            %         else
            %             imPath=[g.videoFolder 'normal/frame' num2str(fList(k),g.decSpec) '.bmp'];
            %             im=imread(imPath);
            %         end
            %         if g.resize
            %             im=imresize(im,g.scaling);
            %         end
            %         tex=Screen('MakeTexture',g.window,im);
            %         g.texes(fList_handles(k))=tex;
            %         g.prevTextures=[g.prevTextures tex];
        end
    end
    time=GetSecs-Tstart;
    period=time/length;
    %fprintf('%d textures loaded in %d, %d per tex\n',length,time,period);
end
end

function nQ=EnqueueFrame(n)
global g;
g.nQ=g.nQ+1;
%g.FrameNumsQ=[g.FrameNumsQ n];
g.FrameNumsQ(g.nQ+1)=n;
 %Append produces a row vec
nQ=g.nQ;
end



function response=makeResponse(direction,sampleRate)
global g; %Comes up with a fake response that sets a trend

thresh=sampleRate/10;


if rand <thresh
    correct=1;
else
    correct=0;
end

correct=0;

if direction==0 %forwards
    if correct==0 response=0; else response=1;end
elseif direction==1 %backwards
    if correct==0,response=1;else response=0;end
end

end


function radians = deg2rad(degrees)
radians = degrees * pi / 180;
return
end

function []=testClip(start,length,negative,backwards,interval,regime,location,chromatic)
cacheClip(start,length,negative,backwards,interval,regime,location,chromatic);
playClip(start,length,negative,backwards,interval,regime,location,chromatic);
end

function loadLayout(layPath)

end

function makeSureAtLeastN(n)
global g;
if ~g.sim
while g.FramesInVM<g.maxFramesInVM
    
    LoadNextFrame();
end
end
end

function CacheForSeconds(n)
t=tic;
while toc(t)<n
    LoadNextFrame();
end
end

function CacheN(n)
t=tic;
for i=1:n
    LoadNextFrame();
end
end

function LoadNextFrame()
global g;

%Check whether video memory is full or not before here!
if ~g.sim
if g.nextFrI <= size(g.FrameNumsQ,2) &&  g.FramesInVM < g.maxFramesInVM%It's a vector, so must reduce to one first!
    LoadFrame(g.FrameNumsQ(g.nextFrI),0,0);
    g.FramesInVM=g.FramesInVM+1;
    g.nextFrI=g.nextFrI+1;
end
end
end

function LoadFrame(n,negative,chromatic)
global g;
t=GetSecs();
if negative
    imPath=[g.videoFolder 'negative\frame' num2str(n,g.decSpec) '.bmp'];
end
if chromatic
    imPath=[g.videoFolder 'hue\frame' num2str(n,g.decSpec) '.bmp'];
    im=imread(imPath);
    im=uint8(im);
    im=im*255;
else
    imPath=[g.videoFolder 'normal/frame' num2str(n,g.decSpec) '.bmp'];
    im=imread(imPath);
end
if g.resize
    im=imresize(im,g.scaling);
end
tex=Screen('MakeTexture',g.window,im);
g.texes(g.nextFrI)=tex;
g.numsByTex(g.nextFrI)=n;
% g.testTex=tex;
t=GetSecs()-t;
fprintf('Frame q %d loaded in %d to tex %d, %d in VM\n',n,t,tex,g.FramesInVM);
end

function [  ] = cacheClip( start,length,negative,backwards,interval,regime,location,chromatic )
global g;
%Plays a constant clip starting at start with a certain length
%Regime controls INVERSION. 1=up 2=inverted
fprintf('Caching clip...\n');
Tstart=GetSecs();
if g.sim
    return
end

clear timestamps
fList=0:length-1;
fList=interval.*floor(fList/interval);

if backwards
    fList=-fList;
    fList=fList+(start+length-1);
else
    fList=fList+start;
end;
fList_handles=fList;
if negative, fList_handles=fList+g.frames;end;
if chromatic, fList_handles=fList+2*g.frames;end;
fList;

%Now we have the final frame list, so can cache.

nTexes=size(fList,2);

if ~g.sim
    for k=1:nTexes
        LoadFrame(fList(k),fList_handles(k),negative,chromatic);
        %         k;
        %         if negative
        %             imPath=[g.videoFolder 'negative\frame' num2str(fList(k),g.decSpec) '.bmp'];
        %         end
        %         if chromatic
        %             imPath=[g.videoFolder 'hue\frame' num2str(fList(k),g.decSpec) '.bmp'];
        %             im=imread(imPath);
        %             im=uint8(im);
        %             im=im*255;
        %         else
        %             imPath=[g.videoFolder 'normal/frame' num2str(fList(k),g.decSpec) '.bmp'];
        %             im=imread(imPath);
        %         end
        %         if g.resize
        %             im=imresize(im,g.scaling);
        %         end
        %         tex=Screen('MakeTexture',g.window,im);
        %         g.texes(fList_handles(k))=tex;
        %         g.prevTextures=[g.prevTextures tex];
    end
end
time=GetSecs-Tstart;
period=time/length;
fprintf('%d textures loaded in %d, %d per tex\n',length,time,period);
end


function [  ] = playClip(start,length,negative,backwards,interval,angleCode,location,chromatic )
global g;
%The start here is a number in the Q!
fprintf('Playing clip at q %d length %d\n',start,length);
%Real time
%Plays a constant clip starting at start with a certain length
%Regime controls INVERSION - 1 upright 2 inverted
%Sampling rate is shown as the number of frames we sample every -
%1=fullrate, 10=every ten frames.
if g.sim
    return
end
Priority(2);
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
if chromatic, fList=fList+2*g.frames;end;
fList;

if angleCode==1
    angle=0;
elseif angleCode==2
    angle=90;
elseif angleCode==3
    angle=270;
elseif angleCode==4
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
fList
thisMissed=0;
thisNInQ=g.nextFrI;
if ~g.sim
    it=1; %Should always place the iterator declaration right before the loop starts
    for k=fList
        Screen('FillRect',g.window,[0 0 0],rect);
        % Screen('DrawTexture', g.window, g.texes(k),[],rect,angle);
        Screen('DrawTexture', g.window, g.texes(k),[],rect,angle);
        fprintf('f %d\n',g.numsByTex(k));
        %k is now the number of this frame in the queue
        
        %New method with entire framerate calc at beginning
        if it==1
            
            % Update display:
            %Flip is asynchronous unless we ask it not to be. Here it
            %doesn't return until the flip has happened.
            t1=GetSecs();
            [VBLTimestamp StimulusOnsetTime FlipTimestamp Missed Beampos] = Screen('Flip', g.window);
               Screen('Close',g.texes(k));
            g.FramesInVM=g.FramesInVM-1;
            t2=GetSecs();
            % fprintf('Free time was %d\n',t2-t1);
            timestamps(it)=VBLTimestamp;
            nextStamps=VBLTimestamp + intervals;
        else
            
            % [VBLTimestamp StimulusOnsetTime FlipTimestamp Missed Beampos] =
            % Screen('Flip', g.window,nextStamp); %Original
           if  g.intercache
               t1=GetSecs();
            expFree=nextStamps(it)-t1;
            if expFree>0.015;
                LoadNextFrame();
                %  fprintf('Doing work\n');
            end
            end
            [VBLTimestamp StimulusOnsetTime FlipTimestamp Missed Beampos] = Screen('Flip', g.window,nextStamps(it));
              Screen('Close',g.texes(k));
            g.FramesInVM=g.FramesInVM-1;
            t2=GetSecs();
            
            
            % fprintf('Free time exp. %d, was %d, missed %d\n',expFree,t2-t1,Missed);
            if Missed>0
                thisMissed=thisMissed+1;
            end
            % Release texture:Screen('Flip', window);
            timestamps(it)=VBLTimestamp;
            nextStamp=VBLTimestamp+ g.int;
        end
        
        it=it+1;
    end
    g.missed=g.missed+thisMissed;
    fprintf('Missed this time: %d, %d total. Loaded %d\n',thisMissed,g.missed,g.nextFrI-thisNInQ);
    Screen('FillRect',g.window,g.grey, rect);
    Screen(g.window, 'Flip');
    
    
    
    for i=1:size(timestamps,2)-1
        diffs(i)=timestamps(i+1)-timestamps(i);
    end
    %     'Framerate'
    %     1/mean(diffs)
    %     std(diffs)
    Priority(0);
end

end





%------------------------------------BEGIN COMMON EXPERIMENTAL HELPER FUNCTIONS

function []=flip()
global g;
if ~g.respSim
    Screen('Flip',g.window);
end
end

function[]=fKbWait()
global g;
if ~g.sim, KbWait; end
end

function [] = InitialiseFramework()
%Sets up PsychToolbox and other parts of the experimental framework
global g;

dbstop if error

g.training=1;
if g.sim, g.pause=0; g.smallPause=0; end
if ~g.debugMode
    g.subjectName=input('Enter subject initials: ','s');
else
    g.subjectName='test';
end
g.logFileName=[g.subjectName '_log.txt'];
g.logFileID=fopen(['output\' g.logFileName],'w');
assert(g.logFileID ~= -1); %Check file is really open
fireLog(['Starting experiment on subject ' g.subjectName]);
g.responseFileName=[g.subjectName '_results.txt'];
g.responseMatFileName=[ g.subjectName '_results.mat'];
g.responseFileID=fopen([ 'output\' g.responseFileName],'w');
g.FrameNumsQ=zeros(1,1400);
g.nextFrI=1;
g.FramesInVM=0;

%Check cluster connection
if g.komodo
    if ~g.debugMode
        try
            
            if ~exist(g.clusterOutputFolder)
                error('Cannot mount output folder!');
            end
        catch
            error('Something wrong with output folder - !');
        end
    end
end
if ~g.debugMode;  HideCursor(); end

if ~g.sim    
    AssertOpenGL;
    screenid = max(Screen('Screens'));
    PsychImaging('PrepareConfiguration');
    PsychImaging('AddTask', 'General', 'FloatingPoint32BitIfPossible');
    
    Screen('Preference', 'SkipSyncTests', 1)
    
    screens=Screen('Screens');
    screenNumber=max(screens);
    screenNumber=min(screenNumber,g.monitor);
    if g.debugMode
        %  [g.window,windowRect]=Screen(screenNumber,'OpenWindow',0,[20 20 800 800],[],2);
        % [g.window,g.windowRect]=Screen(screenNumber,'OpenWindow',[127.5 127.5 127.5],[20 20 800 800],[],2);
        [g.window,g.windowRect]=PsychImaging('OpenWindow',screenNumber,128,[20 20 800 800]);   
    else
        % [g.window,g.windowRect]=Screen(screenNumber,'OpenWindow',[127.5 127.5 127.5],[],[],2);
        [g.window,g.windowRect]=PsychImaging('OpenWindow',screenNumber,128);
    end
    
    %pixelSize=Screen('PixelSize', windowPtrOrScreenNumber);
    g.centreX=(g.windowRect(3)+g.windowRect(1))/2;
    g.centreY=(g.windowRect(4)+g.windowRect(2))/2;
    g.black=BlackIndex(g.window);
    white = WhiteIndex(g.window);  % Retrieves the CLUT color code for white.
    g.grey = (g.black + white) / 2;
    Screen('FillRect',g.window,g.grey);
end
%Screen('BlendFunction',g.window,GL_ONE,GL_ZERO);
HideCursor;
end

function fireClear()
global g;
Screen('FillRect',g.window,g.grey);
Screen('Flip',g.window);
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

function firePause(s)
%Log a string to the logfile for this experiment
global g;
if ~g.sim 
    pause(s);
end

end

function respLog(s)
%Log a string to the logfile for this experiment
global g;

fprintf(g.responseFileID,[s '\n']);
end

function respSave(trialParams)
%Log a string to the logfile for this experiment
global g;
design=g.design;
save([g.responseMatFileName],'trialParams','design');
end

function fireText(s)
global g;
if ~g.sim
    DrawFormattedText(g.window,s,'center','center',g.textColour)
    Screen('Flip',g.window);
    fKbWait;
    pause(1);
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

function fireTextConfirm(s)
global g;
if ~g.sim
    DrawFormattedText(g.window,s,'center','center',g.textColour)
    Screen('Flip',g.window);
    fKbWait;
    pause(1);
end
end

function fireTextNoWait(s)
global g;
if ~g.sim
    DrawFormattedText(g.window,s,'center','center',g.textColour)
    Screen('Flip',g.window);
end
end

function []=fixationSpotAdjustNoFlip(varargin)
global g;
if nargin==0, colour=[0 0 0];
else colour=varargin{1}; end

if ~g.sim
    % Draw a fixation point
    x=g.centreX-10;
    y=g.centreY-10;
    Screen('DrawDots', g.window, [x y], 12, [255 255 255], [], 2);
    Screen('DrawDots', g.window, [x y], 8, colour, [], 2);
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
global g;
fclose all
sca
ShowCursor;
copyfile(['output\' g.logFileName],[g.clusterOutputFolder g.logFileName]);
copyfile(['output\' g.responseFileName],[g.clusterOutputFolder g.responseFileName]);
keyboard
end

function y=trueWithProb(x)
%Returns 1 with prob x
if x==0, y=0;return; end
if x==1, y=1;return; end
n=rand;
if n<x, y=1;return; else y=0;return; end
end

function waitForX()
global g;
if g.sim
else
    a = 0;
    while a == 0
        %WaitSecs(0.2);
        KbWait;
        [ keyIsDown, kt, keyCode ] = KbCheck;
        %keyCode is a 256-element list saying which keys are down, thus we must
        %not truncate it
        response_key = KbName(keyCode);
        %             a = strcmp(response_key, 'RightArrow');
        %             b = strcmp(response_key, 'LeftArrow');
        %If several keys are down, select only one
        if strcmp(class(response_key),'cell'), response_key=response_key{1}; en
            if 1
                if strcmp(response_key,'q')
                    exitExperiment();
                end
            end
            a = strcmp(response_key, 'x') | strcmp(response_key,'X') ; %Mac
            if a(1)==1; a=1; else a=0; end
            %Covers the case where a is a vecto
            %
            if g.debugMode
                if strcmp(response_key, 'k')
                    keyboard
                end
            end
        end
    end
end
end

function [r]=getYN()
%Gets key response - 0 for left or 1 for right.
global g;
if g.sim
    if rand < 0.3
        r=1;
    else
        r=0;
    end
    
else
    
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
        a = strcmp(response_key, 'up') | strcmp(response_key,'UpArrow') ; %Mac
        if a(1)==1; a=1; else a=0; end
        %Covers the case where a is a vector
        
        
        %Makes sure a is always 1 or 0 and not anything else
        b = strcmp(response_key, 'down') | strcmp(response_key,'DownArrow');
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
end

function [ allTrialParams blockParamCells info params ] = fireReadParams(  )
%Generates block params for experiment, and runs
%ASSUMPTIONS: we have even block structure (same amount of trials per
%block)
%QUESTION: do we use matrices or cells for all this? Would be easier to
%work in numbers, then use cells if we need to...
%However, the advantage to cells is that they don't need to be
%square/rectangular. Although we can use zeros... hm...

%Dummy params

%Here we vary target length across blocks, sample length within.
%INBLOCK or ACROSSBLOCK

params(1).name='inv';
params(1).type='enum';
params(1).list=[1 2];
params(1).scheme='inblock';


params(2).name='Lpre';
params(2).type='enum';
params(2).list=[0 25 50 100];
params(2).scheme='inblock';

params(3).name='Lpost';
params(3).type='enum';
params(3).list=[0 25 50 100];
params(3).scheme='inblock';

blockReps=10; %Number of repetitions of each style of block
conditionReps=40; %Number of repetitions of each condition. Must be divisible by blockReps
conditionRepsPerBlock=conditionReps/blockReps;

if mod(conditionReps,blockReps)
    error('conditionReps is not divisble by blockReps');
end

%Get cardinalities
nFactors=size(params,2);
nConditions=1;
for i=1:nFactors
    params(i).n=size(params(i).list,2);
    nConditions=nConditions*params(i).n;
end
%We now know nConditions.
nBlockTypes=1;
for i=1:nFactors
    if strcmp(params(i).scheme,'acrossblock')
        nBlockTypes=nBlockTypes*params(i).n;
    end
end

%Work out block structure
nBlocks=nBlockTypes*blockReps;
%nTrialsPerBlock=
nTrials=nConditions*conditionReps;
nTrialsPerBlock=nTrials/nBlocks;
if mod(nTrials,nBlocks)
    error('trials do not divide into blocks!');
end
trials=zeros(nBlocks,nTrialsPerBlock,nFactors);
%Work out block params
for f=1:nFactors
    if strcmp(params(f).scheme,'acrossblock')
        blockParams{f}=randomiseWithReps(params(f).list,blockReps);
    end
end
%Of course we need to ensure nonindependence between block params...
%Work out params per trial


%Work out our groups of inblock params which must covary
levelLibrary=[];

for f=1:nFactors
    if strcmp(params(f).scheme,'inblock')
        %Criterion for joining the set
        levelLibrary=multiplyOutTransitive(levelLibrary,params(f).list);
    end
end
%Now we have levelLibrary - we still need to account for the condition reps
nLibraryConditions=size(levelLibrary,2);
conditionRepsPerBlock=nTrialsPerBlock/nLibraryConditions;
levelLibrary=repmat(levelLibrary,1,conditionRepsPerBlock);

%WIP - we know LevelLibrary is OK at this point.

it=1;
for b=1:nBlocks
    %Randomise indices for the group of block params - with group
    %validity!
    
    indices=1:nTrialsPerBlock;
    indices=Shuffle(indices);
    libraryRow=1; %We need this to index rows in LL, which are not the same as f
    for f=1:nFactors
        if strcmp(params(f).scheme,'inblock')
            values=levelLibrary(libraryRow,:);
            for i=1:nTrialsPerBlock
                thisPList(i)=values(indices(i));
            end
            blockParamLists{b,f}=thisPList;
            libraryRow=libraryRow+1;
        end
    end
    for f=1:nFactors
        if strcmp(params(f).scheme,'acrossblock')
            blockParamLists{b,f}=repmat(blockParams{f}(b),1,nTrialsPerBlock);
        end
    end
end

%Work out all trials
it=1;
for b=1:nBlocks
    for t=1:nTrialsPerBlock
        for f=1:nFactors
            allTrialParams{it,f}=blockParamLists{b,f}(t);
            thisBlockParams{t,f}=blockParamLists{b,f}(t);
        end
        it=it+1;
    end
    blockParamCells{b}=thisBlockParams;
end


%Sanity checking
assert(nTrials==nBlockTypes*blockReps*nTrialsPerBlock);

fprintf('--Experiment summary:\n');
fprintf('%d factors\n%d conditions\n%d reps per condition\n',nFactors,nConditions,conditionReps);
fprintf('%d block types\n',nBlockTypes);
fprintf('%d block reps\n',blockReps);
fprintf('%d blocks\n',nBlocks);
fprintf('%d trials per block\n',nTrialsPerBlock);
fprintf('%d trials\n',nTrials);

for s=[5 10 15]
    fprintf('%d sec trial -> %d min (%2.2f hr) experiment\n',s,(nTrials*s)/60,(nTrials*s)/3600);
end

%Load info
info.nTrials=nTrials;
info.nBlocks=nBlocks;
info.nTrialsPerBlock=nTrialsPerBlock;
end

function v=randomiseWithReps(in,reps)
v=repmat(in,1,reps);%Assume row vecs
v=Shuffle(v);
end

function x=crossExpandSeveral(varargin)
%Shuffle several variables so that they covary properly
x=varargin{1};
for i=2:nargin
    x=multiplyOutTransitive(x,varargin{i});
end
end


function z=multiplyOut2(x,y)
fassert(rowvec(x));
fassert(rowvec(y));
sx=size(x,2);
sy=size(y,2);
%Row vectors
z=repmat(x,1,size(y,2));
for i=1:size(z,2)
    z(2,i)=y(ceil(i/sx)); %Cunning
end
end

function z=multiplyOutTransitive(x,y)
%x: EXISTING MATRIX
%y: NEW MATRIX
%Repeat x n(y) times. Each x has y(i) repeated under it.
%fassert(rowvec(y));
%Make suitable for recursion:
if isempty(x)
    z=y;
    return
end
sx=size(x,2);
rowsx=size(x,1);
sy=size(y,2);
%Row vectors
z=repmat(x,1,size(y,2));
for i=1:size(z,2)
    z(rowsx+1,i)=y(ceil(i/sx)); %Cunning
end
end


