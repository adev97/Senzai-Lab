
function SquareMoveContinuous(trials,s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function generates black or white square sequences with parameters
%defined by the array of trials structures (see trialStruct.m). Trials
%structures are automatically generated from the table values in the gui
%by trialsStruct.m so your stimulus should take only one input namely
%trials. You can access parameters of a structure in the trials structure
%array using dynamic field referencing (e.g. trials(1).Orientation ...
%returns the orientaiton of trial 1). As you write your stimulus you can
%test it by creating a Default trials structure as done below so you can
%see if it is behaving as expected before adding it to the stimGen gui.
%
% INPUTS:  TRIALSSTRUCT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by MSC 4-23-12 (Modified from DriftDemo2 in PTB)
% Modified by: MSC/2012-4-27,
% Modifid by: YS 2018-03-23,
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%% DEFAULTS FOR TESTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%UNCOMMENT THIS SECTION FOR RUNNING STIMULUS AS STAND ALONE; COMMENT ABOVE
%CONFLICTING FUNCTION FULLFIELDGRATING(TRIALS)
% function [trials] = FullFieldGrating(stimType,table)
% if nargin<1
%     table = {'Square Size (deg)', 5, 1, 5;...
%               'Square PositionX (deg)', -22.5, 5, 22.5;...
%               'Square PositionY (deg)', -22.5, 5, 22.5;...
%               'Square Luminance (binary)', 0, 1, 1;...
%               'Timing (delay,duration,wait) (s)', 0.1, 0.1, 0.1;...
%               'Blank', 0, [], [];
%               'Randomize', 1, [], [];...
%               'Interleave', 0, [], [];...
%               'Repeats', 1, [], [];...
%               'Initialization Screen (s)', 5, [],[]};
%    stimType = 'Receptive Field Mapping';
%
% end
% trials = trialStruct_RFmapS_Yuta(stimType, table);

%%%%%USE WITH CAUTION%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Screen('Preference', 'SkipSyncTests', 1);%% better to be commented
%%Screen('Preference', 'SkipSyncTests', 0);%% use this after above
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%


% Get monitor info from monitorInformation located in RigSpecificInfo dir.
% This structure contains all the pertinent monitor information we will
% need such as screen size and appropriate conversions from pixels to
% visual degrees
monitorInformation;


%%%%%%%%%%%%%%%%%%%%% TURN OFF PTB SYSTEM CHECK REPORT %%%%%%%%%%%%%%%%%%%%
Screen('Preference', 'Verbosity',1);
% This will suppress all but critical warning messages
% At the end of the code we will return the verbosity back to norm level 3

% please see the following page for an explanation of this function
% http://psychtoolbox.org/FaqWarningPrefs
% NOTE: as you debug your code comment this line because PTB will return
% back useful info about memory usage that will tell you about leaks that
% may casue problems

% When Screen('OpenWindow',w,color) is called, PTB performs many checks of
% your system. The time it takes to perform these checks depends on the
% noisiness of your system (up to two seconds on 2-photon rig). During this
% time it displays a white screen which is obviously not good for visual
% stimulation. We can disable the startup screen using the following. The
% sreen will now be black before visual stimulus
Screen('Preference', 'VisualDebuglevel', 3);
% see http://psychtoolbox.org/FaqBlueScreen for a reference


    
%%%%%%%%%%%%%%%%%%%%% OPEN A SCREEN & DETERMINE PARAMETERS %%%%%%%%%%%%%%%%
% Use a try except block to prevent the screen from hanging. During testing
% if the screen does hang press cntrl C or cntrl-alt del to bring up the
% task manager to stop PTB execution
try
    
    % Require OPENGL becasue some of the functions used here need the
    % OPENGL version of PTB
    AssertOpenGL;
    
    %%%%%%%%%%%%%%%%%%%%%% GET SPECIFIC MONITOR INFORMATION %%%%%%%%%%%%%%%%%%%
    
    % SCREEN WE WILL DISPLAY ON
    %Query monitorInformation for screenNumber
    screenNumber = monitorInfo.screenNumber;
    
    % COLOR INFORMATION OF SCREEN
    % Get black, white and gray color values for the current monitor
    whitePix = WhiteIndex(screenNumber);
    blackPix = BlackIndex(screenNumber);
    
    %Convert balck and white to luminance values to determine gray
    %luminance
    whiteLum = PixToLum(whitePix);
    blackLum = PixToLum(blackPix);
    grayLum = (whiteLum + blackLum)/2;
    
    % Now determine the pixel value of gray from the gray luminance
    grayPix = GammaCorrect(grayLum);
    
    
    % CONVERSION FROM DEGS TO PX AND SIZING INFO FOR SCREEN
    %conversion factor specific to monitor
    degPerPix = monitorInfo.degPerPix;
    
    RightBoxBG = grayPix;
    LeftBoxBG = grayPix;
    
    
    
    % CONVERSION FROM DEGS TO PX AND SIZING INFO FOR SCREEN

    % Size of the grating (in pix) that we will draw (1.5 times 
    % monitor width)
    visibleSize = 1*monitorInfo.screenSizePixX;

    X_shiftperframe = round(1/degPerPix *40/60); 
    Y_shiftperframe = round(1/degPerPix *10/60); 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%% INITIAL SCREEN DRAW %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % We start with a gray screen before generating our stimulus and displaying
    % our stimulus.
    
    % HIDE CURSOR FROM SCREEN
    HideCursor;
    % OPEN A SCREEN WITH A BG COLOR OF GRAY (RETURN POINTER W)
    [w, screenRect]=Screen('OpenWindow',screenNumber, grayPix);
    LeftBoxStim(w, screenRect, LeftBoxBG);
    RightBoxStim(w, screenRect, RightBoxBG);  
    
    % CREATE A DESTINATION RECTANGLE where the stimulus will be drawn to
    
    % dstRect=[0 0 visibleSize visibleSize];
    
    dstRect=[0 0 monitorInfo.screenSizePixX monitorInfo.screenSizePixY];
    %center the rectangle to the screen
    dstRect=CenterRect(dstRect, screenRect);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%% PREP SCREEN FOR DRAWING %%%%%%%%%%%%%%%%%%%%%%%%%
    
    % SCRIPT PRIORITY LEVEL
    % Query for the maximum priority level availbale on this system. This
    % determines the priority level of the matlab thread (0= normal,
    % 1=high, 2=realTime priority) note that a setting of 2 may cause the
    % keyboard to be unresponsive. You may want to play with this number if
    % you have trouble recovering the screen back
    
    priorityLevel=MaxPriority(w);
    Priority(priorityLevel);
    
    % INTERFRAME INTERVAL INFO
    % Get the montior inter-frame-interval
    ifi = Screen('GetFlipInterval',w);
    
    %on old slow machines we may not be able to update every ifi. If your
    %graphics processor is too slow you can buy a better one or adjust the
    %number of frames to wait between flips below
    
    waitframes = 1; %I expect most new computers can handle updates at ifi
    ifiDuration = waitframes*ifi;
    %
    % % CREATE A DESTINATION RECTANGLE where the stimulus will be drawn to
    %     dstRect=[0 0 monitorInfo.screenSizePixX monitorInfo.screenSizePixY];
    %     %center the rectangle to the screen
    %     dstRect=CenterRect(dstRect, screenRect);
    
    %%%%%%%%%%%%%%%%%%%%%% DRAW PRESTIM GRAY SCREEN %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % We call the function stimInitScreen to draw a screen to the window before
    % the stimulus appears to allow for any adaptation that is need to a change
    % in luminance
    stimInitScreen(w,trials(1).Initialization_Screen,grayPix,ifiDuration)
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%% CONSTRUCT AND DRAW TEXTURES %%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This is the main body of the code. We will loop through our trials array
    % structure, construct a grating texture based on the values for each trial
    % and then execute the drawing in a while loop. All of this must be done in
    % a single loop becasue we need to close the textures in the trial loop
    % after using each texture becasue otherwise they will hang around in
    % memory and cause the familiar Java runtime error: Out of memory.
    
    % Exit Codes and initialization
    
    % This is a flag indicating we need to break from the trials
    % structure loop below. The flag becomes true (=1) if the
    % user presses any key
    exitLoop=0;
    
    
    SquarePosXPool = [trials(1).Square_PositionX_vector; 100];
    SquarePosYPool = [trials(1).Square_PositionY_vector; 100];
    % MAIN LOOP OVER TRIALS TO CONSTRUCT TEXTURES AND DRAW THEM
    for trial=1:numel(trials)
        if exitLoop==1
            break;
        end
        
        %%%%%%%%%%%%%%%%%%%% GET STIMULUS TIMING INFORMATION %%%%%%%%%%%%%%%%%%%%%%
        % The wait, duration, and delay are stored in trials structure. They
        % may vary over the trials if an LED was shown so get them for each
        % trial
        
        if trials(trial).Square_PositionX==trials(trial).Square_PositionX_vector(1)...
                && trials(trial).Square_PositionY==trials(trial).Square_PositionY_vector(1)
            delay(trial) = trials(trial).Timing_End(1);
            duration(trial) = trials(trial).Timing_End(2);
            wait(trial) = trials(trial).Timing_End(3);
        else
            delay(trial) = trials(trial).Timing(1);
            duration(trial) = trials(trial).Timing(2);
            wait(trial) = trials(trial).Timing(3);
        end
        
        TexIdx(trial) = find(SquarePosXPool==trials(trial).Square_PositionX ...
            & SquarePosYPool==trials(trial).Square_PositionY);
    end
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%% CONSTRUCT STIMULUS TEXTURES %%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % To make a static grating texture of a drifting grating we need the
    % contrast and the spatial frequency. For each trial in our structure we
    % will get these two variables and convert them to appropraite units and
    % then make our texture.
    
    flagLuminance = trials(1).Square_Luminance;
    if flagLuminance==0
        SquareLum = blackPix;
    elseif flagLuminance==1
        SquareLum = whitePix;
    end
    degSize = trials(1).Square_Size;
    
%         for textureIdx=1:numel(SquarePosXPool)
    
            % Get the contrast and spatial frequency of the trial
    
    %         degPositionX = SquarePosXPool(textureIdx);
    %         degPositionY = SquarePosYPool(textureIdx);
            degPositionX = -30; %  0
            degPositionY = -30; % 0
    
            % convert to pixel units
            pxSize = ceil(degSize /degPerPix);
            pxPositionX = ceil(degPositionX /degPerPix);
            pxPositionY = ceil(degPositionY /degPerPix);
            pxHalfSize = ceil(pxSize/2);
    
            % compute square pixXlim and pixYlim in pix
            pixXlim(1) = ceil(monitorInfo.screenSizePixX/2) + pxPositionX - pxHalfSize;
            pixXlim(2) = pixXlim(1) + pxSize;
            pixYlim(1) = ceil(monitorInfo.screenSizePixY/2) + pxPositionY - pxHalfSize;
            pixYlim(2) = pixYlim(1) + pxSize;
    
            imMtrx = grayPix*ones(monitorInfo.screenSizePixY,monitorInfo.screenSizePixX);
    
            for i=1:monitorInfo.screenSizePixX
                for j=1:monitorInfo.screenSizePixY
                    if i>=pixXlim(1) && i<=pixXlim(2) ...
                            && j>=pixYlim(1) && j<=pixYlim(2)
    
                        imMtrx(j,i)=SquareLum;
                    end
                end
            end
%             squaretex{textureIdx}=Screen('MakeTexture', w, imMtrx);
    
%         end
    
    squaretex = Screen('MakeTexture', w, imMtrx);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DRAW TEXTURES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % In DRAW TEXTURES, we will obtain specific parameters such as the
    % orientation etc for each trial in the trials struct. We will then draw an
    % initial gray screen persisting for a time called delay. Then we will draw
    % our grating using the parameters we pulled from the trials structure.
    % Lastly we will draw another gray screen persisting for a time called
    % wait. We repeat until the end of trials.
    for trial=1:numel(trials)
        if exitLoop==1
            break;
        end
        
        %%%%%%%%%%%%%%%%%%%% DRAW DELAY GRAY SCREEN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % DEVELOPER NOTE: Prior versions of stimuli used the func WaitSecs to
        % draw gray screens. This is a bad practice because the function sleeps
        % the matlab thread making the computer unresponsive to KbCheck clicks.
        % In addition PTB only guarantees the accuracy of WaitSecs to the
        % millisecond scale whereas VBL timestamps described below uses
        % GetSecs() a highly accurate submillisecond estimate of the system
        % time. All times should be referenced to this estimate for better
        % accuracy.
        
        % We start by performing an initial screen flip using Screen, we return
        % back a time called vbl. This value is a high precision time estimate
        % of when the graphics card performed a buffer swap. This time is what
        % all of our times will be referenced to. More details at
        % http://psychtoolbox.org/FaqFlipTimestamps
        vbl=Screen('Flip', w);
        %LeftBoxStim(w, screenRect, LeftBoxBG);
        %RightBoxStim(w, screenRect, RightBoxBG);
        % The first time element of the stimulus is the delay from trigger
        % onset to stimulus onset
        delayTime = vbl + delay(trial);
        
        % Display a gray screen while the vbl is less than delay time. NOTE
        % we are going to add 0.5*ifi to the vbl to give us some headroom
        % to take possible timing jitter or roundoff-errors into account.
        while (vbl < delayTime)
            % Draw a gray screen
            Screen('FillRect', w,grayPix);
            
            %             if trials(trial).Audio_Cue==1
            %                 TopLeftBoxStim(w, screenRect, blackPix);
            %             end
            % update the vbl timestamp and provide headroom for jitters
            vbl = Screen('Flip', w, vbl + (waitframes - 0.5) * ifi);
            
            % exit the while loop and flag to one if user presses any key
            if KbCheck
                exitLoop=1;
                break;
            end
        end
        n=0; % This is a counter to shift our grating on each redraw
        
%%%%%%%%%%%%%%%%%%%%%%%%%% DRAW GRATING TEXTURE %%%%%%%%%%%%%%%%%%%%%%%%%%%         
        % If the trial is a blank then we do not need to set src and dst
        % rect and calculate grating shifts etc
   
            vbl=Screen('Flip', w);
            
            % Set the runtime of each trial by adding duration to vbl time
            runtime = vbl + duration;
            
            while (vbl < runtime)
                % calculate the offset of the grating and use the mod func
                % to ensure the grating snaps back once the border is
                % reached
                xoffset = -n*X_shiftperframe;
                yoffset = -n*Y_shiftperframe;
                
                n = n+1;
                
                % Set the source rect to excise the grating from
                % srcRect = [xoffset yoffset xoffset+visibleSize yoffset+visibleSize];
                
                srcRect = [xoffset yoffset xoffset+monitorInfo.screenSizePixX yoffset+monitorInfo.screenSizePixY];
                
                % Draw the grating texture for this trial to the dst
                % rectangle
                Screen('DrawTextures', w, squaretex, srcRect,dstRect);
                
%                 % Draw a box at the bottom right of the screen to record
%                 % all screen flips using a photodiode. Please see the file
%                 % FlipCheck.m in the stimulus directory for further
%                 % explanation
%                 FlipCheck(w, screenRect, [whitePix, blackPix], n)
                LeftBoxStim(w, screenRect, blackPix);
                FlipCheck(w, screenRect, [blackPix, whitePix], n);
                
                % update the vbl timestamp and provide headroom for jitters
                vbl = Screen('Flip', w, vbl + (waitframes - 0.5) * ifi);
                
                % exit the while loop and flag to one if user presses any
                % key
                if KbCheck
                    exitLoop=1;
                    break;
                end
            end
  
        
        
        %%%%%%%%%%%%%%%%%%%%% DRAW INTERSTIMULUS GRAY SCREEN %%%%%%%%%%%%%%%%%%%%%%
        % Between trials we want to draw a gray screen for a time of wait
        
        % Flip the screen and collect the time of the flip
        vbl=Screen('Flip', w);
        
    end
    
    % close textures

        Screen('Close',squaretex)

    
    % If keyboard was pressed to escape attempt close of last texture
    if exitLoop
        try
           Screen('Close',squaretex)
        catch
        end
    end
    
    
    
    
    %%%% blank screen after exp %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % outputSingleScan(s,[0 0 0]);
    
    % Flip the screen and collect the time of the flip
    vbl=Screen('Flip', w);
%     LeftBoxStim(w, screenRect, LeftBoxBG);
%     RightBoxStim(w, screenRect, RightBoxBG);
    % We will loop until delay time referenced to the flip time
    waitTime = vbl + 2; % wait 2 [sec]
    %
    while (vbl < waitTime)
        % Draw a gray screen
        Screen('FillRect', w,grayPix);
%         LeftBoxStim(w, screenRect, LeftBoxBG);
%         RightBoxStim(w, screenRect, RightBoxBG);
        % update the vbl timestamp and provide headroom for jitters
        vbl = Screen('Flip', w, vbl + (waitframes - 0.5) * ifi);
        
        %         % exit the while loop and flag to one if user presses any key
        %         if KbCheck
        %             exitLoop=1;
        %             break;
        %         end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Restore normal priority scheduling in case something else was set
    % before:
    Priority(0);
    
    %The same commands wich close onscreen and offscreen windows also close
    %textures. We still need to close any screens opened prior to the trial
    %loop ( the prep screen for example)
    Screen('CloseAll');
    
catch
    %this "catch" section executes in case of an error in the "try" section
    %above.  Importantly, it closes the onscreen window if its open.
    outputSingleScan(s,[0 0 0]);
    
    Screen('CloseAll');
    Priority(0);
    psychrethrow(psychlasterror);
    
end

%%%%%%%%%%%%%%%%%%%%%%%% Turn On PTB verbose warnings %%%%%%%%%%%%%%%%%%%%
Screen('Preference', 'Verbosity',3);
% please see the following page for an explanation of this function
%  http://psychtoolbox.org/FaqWarningPrefs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
java.lang.Runtime.getRuntime().gc % call garbage collect (likely useless)


return

