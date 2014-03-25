% Calibration of CRS colorCAL photometer
% based on rgbInandStimuli cogent program

% Usage: doColorCalCalib( file_name )
% records to the current directory
% - possibility to quit by pressing 'q'

% The psychophysics toolbox needs to be in the path, 
% as well as the visage toolbox (initColorCAL, colorCALread ...)

% April 2011: David Souto @ Dirk Kerzel's lab., Uni Geneva
function doColorCalCalib(file_name)

    DummyMode = 0;   
    if nargin ==0,
       error('doColorCalCalib(file_name)');
    end

    gray              = 192;
    number_of_bits    = 32;  % number of bits per pixel 

    % run on the second screen
    screenID = max(Screen('Screens'));
    [win winRect] = Screen('OpenWindow', screenID, gray, [], number_of_bits, 2, [], [], []);

    [center(1) center(2)] = RectCenter(winRect);
    resx = center(1)*2;
    resy = center(2)*2;
    
    rect = 400; % displayed rectangle pixels
    rect_coords = [center(1)-rect/2 center(2)-rect/2 center(1)+rect/2 center(2)+rect/2];
    
    Screen('TextFont',win,'Ariel');
    Screen('TextSize',win,20);
    Screen('TextStyle',win,1);
    labels  = ['r: '; 'g: '; 'b: ';'l: '];  

    s = sprintf('the resolution is %d x %d\n and the refresh frequency is %d (any key to resume)\n', resx, resy, ceil(1/Screen('GetFlipInterval', win)));   
    Screen('DrawText',win, s, center(1)-500, center(2)-rect+100);         
    Screen('Flip', win); 
    KbWait;
    
	if ~DummyMode,
		if exist('colorCALinit','file'), % check that VSG toolbox is in the path
			%vsgInit; % try this if init does not work
			colorCALinit;
			[ErrorCode]= colorCALautocalibrate;
			CIEcolour = colorCALread;
		else
			   error('you need a pointer in the matlab path directory to the VSG toolbox')
		end
	end
    HideCursor;

    RGBvalue = 0:255; % [0 4 8 16 24 32 64 96 128 140 150 160 170 180 190 200 210 220 230 240 250 255];
    channel  = [1 0 0; 0 1 0; 0 0 1;1 1 1];

    LUMvalue = [];
    CIEcolor = [];
   % VOLTvalue   = [];
    
    isi = 0.2;  
   
    Screen('FillRect', win, [255 255 255],rect_coords);
    Screen('DrawText',win,'press any key to start (press ''q'' to quit)', center(1)-200, center(2)-rect+100);         
    Screen('Flip', win);       
    
    kbWait;
    for j=1:4,
        for i=1:numel(RGBvalue);
            rect_color = [RGBvalue(i) * channel(j,1), RGBvalue(i) * channel(j,2), RGBvalue(i) * channel(j,3)];
            Screen('FillRect', win, rect_color, rect_coords);
            Screen('DrawText',win,[labels(j,:) int2str(RGBvalue(i))], center(1)-100, center(2)-rect+100);         
            
            Screen('Flip',win);
                        
            WaitSecs(isi);
            
            [d d kbCode] = KbCheck();
            if kbCode(KbName('q'));
                sca
                keyboard;
            end
            if ~DummyMode,
                disp('taking measurements');
                CIEcolor = colorCALread;
                LUMvalue(i,j) = CIEcolor(3);
            end
        end
    end   
    
    dlmwrite(file_name,[RGBvalue' LUMvalue],'delimiter',' ');

    Screen('CloseAll');
end
