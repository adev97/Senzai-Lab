% Sync RF mapping data & construct maps from spike rate. Requires manual 
% correction of photodiode (PD) edges

% This script uses an example session

% Please copy the RFMapping_Analysis folder into your folder, and copy
% Example_Session_Mouse08_20251007_810to2250_RFMapping folder and
% Example_Kilosort_Data folder to your local drive
% EB 02.11.2026

%% Input: File_Info_RFmapping row number
rn = 1; % row number in tb
mode_sel = 'mean'; % 'sum' is the default. 'mean' calculates spike rate.
nbins = 25;
save_figs = true;
save_gauss = true;
use_events = true; % Use samplenumbers.npy from events folder instead of continuous data
%%
% Add RF mapping data information to tb. An example session is used here.
p1 = 'R:\Basic_Sciences\Phys\SenzaiLab\Shared\RFMapping_Analysis\';
oe_p = 'D:\OpenEphys_Data\'; %  Example_Session_Mouse08_20251007_810to2250_RFMapping
ks_p = 'D:\Kilosort\'; % Example_Kilosort_Data
tb = readtable(fullfile(p1,'Info_RFmapping.xlsx'));   
session_folder = tb.session_folder{rn};
stim_path = fullfile(p1,'Example_Psychtoolbox_Data\',tb.stim_path{rn}); 
ks_path = fullfile(ks_p,tb.ks_path{rn}); 
gdp = fullfile(p1,'GaussFit_Data\');
gdf = fullfile(gdp,[extractAfter(session_folder,'\'),'.mat']);
gdi = [extractBefore(gdf,'.mat'),'_idx','.mat'];

% Specify probe channel configuration for recording session 
% 1: 4shank_single_columns_0to1440.imro
% 2: 4shank_single_columns_1440to2880.imro
% 3: 4shank_single_columns_2880to4320.imro
% 4: 1shank_single_column_shank4_0to5760.imro 
% 5: Shank1_BankA
% 6: 4shank_150_to_1590um.imro 
% 7: 4shank_single_columns_810to2250.imro
chan_config = tb.chan_config(rn);

% If Kilosort was run manually:
if ~contains(ks_path,'ProbeA')
    fn = extractBefore(extractAfter(ks_path,'Kilosort\'),'\kilosort4');
    channel_map = readNPY(fullfile(ks_path,'channel_map.npy'));   % [nChannels x 1], 0-based
% If Kilosort was run through pipeline:
else
    fn = extractBefore(extractAfter(ks_path,'Kilosort\'),'\ProbeA');
    imro_p = 'R:\Basic_Sciences\Phys\SenzaiLab\Elissa_Belluccini\OpenEphys_Data\imro_files\';
    imro = [imro_p,tb.imro{rn}];
    txt = fileread(imro);
    tok = regexp(txt, '\((\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\)', 'tokens');    
    n = numel(tok);
    data = zeros(n,5);
    for i = 1:n
        data(i,:) = str2double(tok{i});
    end    
    channel_map = data(:,1);  
end

fp = ['R:\Basic_Sciences\Phys\SenzaiLab\Shared\RFMapping_Analysis\Example_RFMapping_Sync_Pulses\',...
    fn,'.mat'];
% prec_files: List of session folder names preceding the recording in stim_path 
prec_files = tb.prec_files(rn);  
protocol = tb.Protocol{rn};

if exist(fp, 'file')==2
    load(fp)
    load(stim_path)
else
    % Insert break point at line 446 (case 'RFMapping') to manually correct
    % photodiode edge detection. Current correction is for the example
    % session.
    [trials, Vq] = Sync_Signals(oe_p,session_folder,stim_path,prec_files,...
        protocol,[],use_events);
     save(fp,'Vq')
end

if exist(gdf,'file')==2
    load(gdf)
else
    gauss_2d_data = [];
end

RFmapping_EB(ks_path,fn,trials,Vq,chan_config,channel_map,mode_sel,nbins,...
    save_figs,save_gauss,gauss_2d_data,gdi)


% Functions
function [stimInfo, Vq] = Sync_Signals(oe_p,cf,stim_path,prec_files,protocol,...
    stim_name,use_events)
    close all
    addpath(genpath("D:\buzcode-master"))
    
    % Synchronizing data using function generator
    % ADC channel 0: Function generator signal
    % ADC channel 1: Photodiode signal
    % Probe channel 385: Function generator signal
    % FG: Function generator, PD: Photodiode
    
    sr = 30000; 
    sr_adc = 30300.5;
    nchan_probe = 385;
    nchan_adc = 12;
    
    % Stimulus file:
    load(stim_path);
    if exist('trials','var')
        stimInfo = trials;
    end
    
    % Load data streams first to find the shortest recording duration (probably 
    % Probe data) over which to align the signals
    fp = fullfile(oe_p,cf,'/Record Node 102');
    xml = xmlread([fp,'/settings.xml']);
    processors = xml.getElementsByTagName('PROCESSOR');
    for i = 0:processors.getLength-1
        p = processors.item(i);
        if strcmp(char(p.getAttribute('name')), 'OneBox')
            nodeId = char(p.getAttribute('nodeId'));
        end
    end

    fp_d0 = sprintf('%s/experiment1/recording1/continuous/OneBox-%s.OneBox-ADC',...
        fp,nodeId);
    cd(fp_d0)

    d0 = LoadBinary('continuous.dat','frequency',sr_adc,'nChannels',12,'channels',1);
    d0_dur = length(d0)/sr_adc;

    snd0_cont = readNPY('sample_numbers.npy');
    snd0_cont = double(snd0_cont);

    fts_d0 = sprintf('%s/experiment1/recording1/events/OneBox-%s.OneBox-ADC/TTL',...
         fp,nodeId);
    tsd0 = readNPY(fullfile(fts_d0,'timestamps.npy'));
    snd0 = readNPY(fullfile(fts_d0,'sample_numbers.npy'));
    snd0 = double(snd0);
    snd0 = snd0-snd0_cont(1);
    snd0t = snd0/sr_adc;

    if ~any(strcmp(protocol,{'RFmapping360','Grating360'}))
        d1_chan = 2; 
    else
        d1_chan = 4;
    end

    d1 = LoadBinary('continuous.dat','frequency',sr_adc,'nChannels',12,'channels',d1_chan);
    d1_dur = length(d1)/sr_adc; % recording duration (sec)
    
    fp_d385 = sprintf('%s/experiment1/recording1/continuous/OneBox-%s.ProbeA',...
     fp,nodeId); 
    cd(fp_d385)
    if ~strcmp(protocol,'RFmapping360')
        d385 = LoadBinary('continuous.dat','frequency',sr,'nChannels',nchan_probe,...
        'channels',385);
         d385_dur = length(d385)/sr;
    end
    snd385_cont = readNPY('sample_numbers.npy');
    snd385_cont = double(snd385_cont);

    fts_d385 = sprintf('%s/experiment1/recording1/events/OneBox-%s.ProbeA/TTL',...
        fp,nodeId); 
    tsd385 = readNPY(fullfile(fts_d385,'timestamps.npy'));
    snd385 = readNPY(fullfile(fts_d385,'sample_numbers.npy'));
    snd385 = double(snd385);
    snd385 = snd385-snd385_cont(1);
    snd385t = snd385/sr;

    % Get the duration (in seconds) of ADC0 and Probe channel 385 recordings 
    % in concatenated file, preceding the current file, in order of recording.
    % The sum of these durations gives the start time of the current
    % recording, relative to the concatenated file.    
           
    pfiles = split(prec_files,',');
    nfiles = length(pfiles);
    if strcmp(prec_files,'NaN')
        nfiles = 0;
    end

    d0_ar = cell(nfiles,1);
    d385_ar = cell(nfiles,1);
    d0_dur_ar = zeros(nfiles,1);
    d385_dur_ar = zeros(nfiles,1);
    
    d0_dur_ar = zeros(nfiles,1);
    d385_dur_ar = zeros(nfiles,1);
    bytesPerSample = 2; % int16 = 2 bytes
 
    for ff = 1:nfiles
        fn = pfiles{ff};
        cf = [cf(1:8),fn];
        fp = fullfile(oe_p, cf,'/Record Node 101');

        xml = xmlread([fp,'/settings.xml']);
        processors = xml.getElementsByTagName('PROCESSOR');
        for i = 0:processors.getLength-1
            p = processors.item(i);
            if strcmp(char(p.getAttribute('name')), 'OneBox')
                nodeId = char(p.getAttribute('nodeId'));
            end
        end

        fp_d0_prec = sprintf('%s/experiment1/recording1/continuous/OneBox-%s.OneBox-ADC/continuous.dat',...
            fp,nodeId); 
        d0_info = dir(fp_d0_prec);
        d0_nSamples = d0_info.bytes / (nchan_adc * bytesPerSample);    
        d0_dur_ar(ff) = d0_nSamples / sr;
        
        fp_d385_prec = sprintf('%s/experiment1/recording1/continuous/OneBox-%s.ProbeA/continuous.dat',...
            fp,nodeId); 
        d385_info = dir(fp_d385_prec);
        d385_nSamples = d385_info.bytes / (nchan_probe * bytesPerSample);
        d385_dur_ar(ff) = d385_nSamples / sr;
    end
    d0_start_time = sum(d0_dur_ar);
    d385_start_time = sum(d385_dur_ar);
    
    %% Function generator signal sent to ADC channel 0 %%
    if ~use_events
        td0 = (1:length(d0))/sr_adc;
        td0 = td0';
        
        figure
        i0 = td0<=d1_dur;
        plot(td0(i0),d0(i0)), hold on
        % plot(td0,d0), hold on
        title('FG signal, ADC0')
        
        % periods: Timepoints of the rising & falling edge of each pulse
        [periods0,in0] = Threshold([td0,double(d0)],'>',10000,'min',0.2);
        d0_midpts = mean(periods0,2); 
        periods0_1col = reshape(periods0.', [], 1);
        
        sweep_time = 20; 
        ns = d1_dur/sweep_time;
        blocks1 = 0:sweep_time:d1_dur-sweep_time;
        blocks2 = sweep_time:sweep_time:d1_dur;
        
        % Function generator: Changes from 1-2 Hz over 20 sec
        nperiods = arrayfun(@(n) sum(periods0(:,1)>blocks1(n) & periods0(:,1)<=blocks2(n)),...
            1:ns, 'UniformOutput', false);
        
        p0_tb = table(periods0(:,1),periods0(:,2),periods0(:,2)-periods0(:,1),...
            'VariableNames',{'RisingEdge','FallingEdge','Difference'});
        
        nperiods = cell2mat(nperiods);
        bf = cumsum(nperiods);
        bs = [1,bf(1:end-1)+1];
        
        % New strategy: Leading edge is the first rising edge following the 
        % smallest width pulse in each 20s sweep. 
        
        [~,min_pts] = arrayfun(@(n) min(p0_tb.Difference(bs(n):bf(n))), 1:length(bs));
        min_pts = min_pts+bs-1;
        start_pts = min_pts+1; % Add 1: The next edge should be the largest pulse
        
        % For each start_pt, keep selecting the next pulse if it is larger
        for i = 1:length(start_pts)
            idx = start_pts(i);
            w = p0_tb.Difference(idx);
            while p0_tb.Difference(idx+1)>w || p0_tb.Difference(idx)<0.4 
                start_pts(i) = idx+1;
                idx = idx+1;
                w = p0_tb.Difference(idx);
            end
        end
    
        % Remove duplicate points:
        [~, ~, ic] = unique(start_pts, 'stable');
        [~, first_occurrence] = unique(ic, 'first');
        all_indices = 1:numel(start_pts);
        dup = setdiff(all_indices, first_occurrence);
        start_pts(dup) = [];
               
        st_adc0 = p0_tb.RisingEdge(start_pts); % Leading pulse rising edge start times
        et_adc0 = p0_tb.RisingEdge(start_pts(2:end)); % Leading pulse rising edge start times
        
        y1 = ones(length(st_adc0),1)*11000;
        plot(st_adc0,y1,'gx','MarkerSize',12)
    
        % Break here & uncomment code to fix edge detection errors:
        % Manually add last edge if missed
        % [x,y] = ginput(1); 
        % x = 1601.695;
        % st_adc0 = [st_adc0;x];
        % et_adc0 = [et_adc0;x];
    else
        st_adc0 = snd0t;
    end    
    
    %% Photodiode signal sent to ADC channel 1 %%
    td1 = (1:length(d1))/sr_adc;
    td1 = td1';
    
    figure
    plot(td1,d1), hold on
    title(sprintf('PD signal, ADC%d',d1_chan-1))
    
    % Threshold function output is: [rising edges, falling edges]
    % Grating: Appears that PD signal is 30 Hz. Width between starting and 
    % falling edge expected to be: 1/30/2 = 0.0167
    % 'min': minimum interval between rising and falling edge for pulse to
    % be included
    % 'max': Intervals between pulses < 'max' will be excluded

    switch protocol
        case 'Grating'
            % PD flickers at 30 Hz for 2 sec grating presentation
            py = 14000;
            [periods1,in1] = Threshold([td1,double(d1)],'>',py,'min',0.01);    
            % Remove any pulses with width < 1.5 s or > 2.5
            too_narrow = diff(periods1,1,2) < 1.5;
            periods1(too_narrow,:) = [];
            too_wide = diff(periods1,1,2) > 2.5;
            periods1(too_wide,:) = [];

            plot(periods1(:,1),ones(length(periods1),1)*py,'gx') 
            plot(periods1(:,2),ones(length(periods1),1)*py,'rx')

            periods1 = fliplr(periods1); % ???

            % Manually add the first falling edge
            % [x1,y] = ginput(1); 
            x1 = 24.0036;

            % Manually add the last rising edge 
            % [x2,y] = ginput(1); 
            x2 = 1000.3832;
            
            periods1 = [[x1;periods1(:,1)],[periods1(:,2);x2]];
            plot(periods1(:,1),ones(length(periods1),1)*py,'bo') 
            plot(periods1(:,2),ones(length(periods1),1)*py,'bo')

            nstim_frames = height(trials);

         case 'Grating360'
            % Sync patch flickers On & Off during 2 sec grating, then Off
            % for 2 sec gray interval. PD signal is overlaid on 60 Hz
            % refresh rate

            threshold = 10000;           
            isAboveThresholdMask = (d1 < threshold); 
  
            figure
            plot(isAboveThresholdMask)
            ylim([-0.1,1.1])
            title('isAboveThresholdMask')

            figure
            plot(td1,isAboveThresholdMask)
            ylim([-0.1,1.1])
            title('isAboveThresholdMask')

            maxGapSamples = 400;
            pd_noGap = fill_short_gaps(isAboveThresholdMask, maxGapSamples);

            figure
            plot(pd_noGap)
            ylim([-0.1,1.1])
            title('pd_noGap')

            figure
            plot(td1,pd_noGap)
            ylim([-0.1,1.1])
            title('pd_noGap')

            pd_noGap_new = fill_short_gaps(~pd_noGap, 400);
            pd_noGap_new = ~pd_noGap_new;

            figure
            plot(td1, ~pd_noGap_new), hold on
            ylim([-0.1 1.1])
            title('PD, processed signal')

            py = 0.9;
            [periods1,in1] = Threshold([td1,double(~pd_noGap_new)],'>',py,'min',0.01);    
            % Remove any pulses with width < 1.5 s or > 2.5
            too_narrow = diff(periods1,1,2) < 1.3;
            periods1(too_narrow,:) = [];
            too_wide = diff(periods1,1,2) > 2.7;
            periods1(too_wide,:) = [];

            plot(periods1(:,1),ones(length(periods1),1)*py,'gx') 
            plot(periods1(:,2),ones(length(periods1),1)*py,'rx')

            % Break here and do manual corrections
            p1 = periods1;
            periods1 = p1(:,2);

            % Manually add the first 2 falling edges
            % [x1,y] = ginput(1); 
            xf1 = 41.8813;
            xf2 = 45.9642;

            periods1 = [xf1;xf2;periods1];

            % Manually add the first and last rising edge 
            % [xr,y] = ginput(1); 
            xr1 = 43.8810;
            xr2 = 1023.230;

            periods1 = [periods1,[xr1;p1(:,1);xr2]];

            % Correct the 4th falling edge:
            % [xf4,y] = ginput(1);
            xf4 =  54.1627;
            periods1(4,1) = xf4;

            plot(periods1(:,1),ones(length(periods1),1)*py,'bo') 
            plot(periods1(:,2),ones(length(periods1),1)*py,'bo')

            nstim_frames = height(trials);
        
        case 'RFmapping'
            tf = 14000;
            tr = 2000;
            [pf,~] = Threshold([td1,double(d1)],'<',tf,'min',0.07);
            [pr,~] = Threshold([td1,double(d1)],'>',tr,'min',0.07);

            plot(pf,ones(length(pf),1).*tf,'gx') 
            plot(pr,ones(length(pr),1).*tr,'rx')
            
            wf = diff(pf,[],2);
            pf(wf>0.15,:) = [];
            wr = diff(pr,[],2);
            pr(wr>0.15,:) = [];

            plot(pf(:,1),ones(length(pf),1).*tf,'go') 
            plot(pr(:,1),ones(length(pr),1).*tr,'ro')  

            %% Break here & uncomment code to fix edge detection errors:

            % Manually add the last rising edge: 
            % (For example session):
            % [x,y] = ginput(1);
            % sprintf('%.3f',x)
            x = 3262.996; 
            pr = [pr;[x,NaN]];
            %%
           
            nstim_frames = height(trials);
            periods1 = [pf(:,1),[pr(:,1)]];

            % Final check of detected PD edges:
            plot(periods1(:,1),ones(length(periods1),1).*tf,'bo') 
            plot(periods1(:,2),ones(length(periods1),1).*tr,'bo') 

        case 'RFmapping360'          
            threshold = 12000;           
            isAboveThresholdMask = (d1 < threshold);

            figure
            plot(isAboveThresholdMask), hold on
            ylim([-0.1 1.1])
            title('isAboveThresholdMask')

            figure
            plot(td1,isAboveThresholdMask), hold on
            ylim([-0.1 1.1])
            title('isAboveThresholdMask')

            maxGapSamples = 400;
            pd_noGap = fill_short_gaps(isAboveThresholdMask, maxGapSamples);            

            figure
            plot(td1,pd_noGap), hold on
            ylim([-0.1 1.1])
            title('pd_noGap')

            pd_noGap_inverted = ~pd_noGap;
            pd_noGap_new = fill_short_gaps(pd_noGap_inverted, 1000);
            pd_noGap_new = ~pd_noGap_new;

            figure
            plot(td1, pd_noGap_new), hold on
            ylim([-0.1 1.1])
            title('PD, processed signal')

            tf = 0.9;
            tr = 0.1;
            [pf,~] = Threshold([td1,double(pd_noGap_new)],'<',tf,'min',0.06);
            [pr,~] = Threshold([td1,double(pd_noGap_new)],'>',tr,'min',0.06);

            plot(pf,ones(length(pf),1).*tf,'gx') 
            plot(pr,ones(length(pr),1).*tr,'rx')

            wf = diff(pf,[],2);
            pf(wf>0.15,:) = [];
            wr = diff(pr,[],2);
            pr(wr>0.15,:) = [];

            %% Manual correction:
            % Make cutoff to remove edges corresponding to next stimulus (grating)
            % cutoff = 2180;
            % pf_after_stim = pf(:,1) > cutoff;
            % pf(pf_after_stim,:) = [];
            % pr_after_stim = pr(:,1) > cutoff;
            % pr(pr_after_stim,:) = [];

            % Add last rising edge
            % % [x,y] = ginput(1); 
            % x = 2176.678;
            % pr = [pr;[x,NaN]];
            %%

            plot(pf(:,1),ones(length(pf),1).*tf,'bo') 
            plot(pr(:,1),ones(length(pr),1).*tr,'bo')              
            wf = diff(pf,[],2);
            wr = diff(pr,[],2);

            periods1 = [pf(:,1),[pr(:,1)]];

            dif_periods1 = diff(periods1,1,2);
            dif_periods1 = [periods1,dif_periods1];


        case 'NaturalScenes' % MNIST or AllenScenes
            % Image display time: 250 ms
            tf = 14000;
            tr = 2000;
            [pf,~] = Threshold([td1,double(d1)],'<',tf,'min',0.2);
            [pr,~] = Threshold([td1,double(d1)],'>',tr,'min',0.2);

            plot(pf,ones(length(pf),1).*tf,'gx') 
            plot(pr,ones(length(pr),1).*tr,'rx')
            
            if ~strcmp(stim_name,'ImageNet')
                wf = diff(pf,[],2);
                pf(wf>0.28,:) = [];
                wr = diff(pr,[],2);
                pr(wr>0.28,:) = [];
            end

            plot(pf(:,1),ones(length(pf),1).*tf,'go') 
            plot(pr(:,1),ones(length(pr),1).*tr,'ro')

            pf = pf(:,1);
            pr = pr(:,1);

            % Manual correction (MNIST, AllenScenes)
            % Add the last falling & rising edge:

            periods1 = [pf,pr]; 
            % nstim_frames = size(stimInfo.image_order,1)*size(stimInfo.image_order,2);

            % [xf,yf] = ginput(1); 
            xf = 343.863;
            % xf = 625.5587;
            % [xr,yr] = ginput(1);
            xr = 344.210;
            % xr = 625.8698;

            pf = [pf;xf];
            pr = [pr;xr];

            % Manual correction (ImageNet)
            % pr = pr(3:end,:);
            % pf = pf(2:end,:);
            
            periods1 = [pf,pr];
            plot(pf,ones(length(pf),1).*tf,'bo') 
            plot(pr,ones(length(pr),1).*tr,'bo')
        
        case 'Movies'
            % Allen movies frame rate: 30 Hz. Image display time: ~33 ms
            % LOC movie: 24 Hz
            tf = 15000;
            tr = 2000;
            % min_ft = 0.02; % Allen movies
            min_ft = 0.03; % LOC movie
            [pf,~] = Threshold([td1,double(d1)],'<',tf,'min',min_ft);
            [pr,~] = Threshold([td1,double(d1)],'>',tr,'min',min_ft);

            plot(pf,ones(length(pf),1).*tf,'gx') 
            plot(pr,ones(length(pr),1).*tr,'rx')
            
            % For Allen movies:
            % wf = diff(pf,[],2);
            % pf(wf>0.07,:) = [];
            % wr = diff(pr,[],2);
            % pr(1:2,:) = [];

            % For LOC movie:
            pr = pr(:,1);
            pf = pf(:,1);
            pr([1,2]) = [];
            pf(1) = [];

            plot(pf(:,1),ones(length(pf),1).*tf,'go') 
            plot(pr(:,1),ones(length(pr),1).*tr,'ro')
            periods1 = [pf,pr];               
    end  
        
    % d1_midpts = mean(periods1,2); 
    % Check Photodiode is synchronized with ADC0: Count number of 
    % photodiode pulses in each sweep (of ADC0 FG signal)
    st0 = st_adc0(1:end-1);
    % PD_pulses_sweep = arrayfun(@(s) sum(d1_midpts>=st0(s) & ...
    %     d1_midpts<et_adc0(s)), 1:length(st0));

    p1_tb = table(periods1(:,1),periods1(:,2),diff(periods1,[],2),...
        'VariableNames',{'RisingEdge','FallingEdge','Difference'});
    
    periods1_vector = reshape(periods1.', [], 1);    
    periods1_vector = periods1_vector(~isnan(periods1_vector));

    % Mouse01_RSC_20250717_Shank1to4_RFmapping long PD pulse at 1696 sec
    % p1d = [periods1_vector,[diff(periods1_vector);NaN]];
    % last_pulse = find(p1d(:,2)>0.103,1);
    % periods1_vector = periods1_vector(1:last_pulse); 


    %% Function generator signal sent to Probe data stream %%
    
    if ~use_events
         td385 = (1:length(d385))/sr;
         td385 = td385';
   
        figure
        i385 = td385<=d1_dur;
        plot(td385(i385),d385(i385)), hold on
        ylim([-0.1,1.1])
        title('FG signal, Probe') 
        
        [periods385,in385] = Threshold([td385,double(d385)],'>',0.8,'min',0.2);
        d385_midpts = mean(periods385,2); 
        periods385_1col = reshape(periods385.', [], 1);

        % Testing: 
        y1 = ones(length(periods385),1)*0.9;
        y2 = ones(length(periods385),1)*0.88;
        plot(periods385,y1,'gx','MarkerSize',12)
        plot(periods385,y2,'rx','MarkerSize',12)
    
        all_periods = table(periods0_1col,periods385_1col,periods0_1col-periods385_1col,...
            'VariableNames',{'periods0','periods385','dif'});
    
        all_ts = table(tsd0,tsd385,tsd0-tsd385,'VariableNames',{'tsd0','tsd385','dif'});
    
        all_sn = table(snd0t,snd385t,snd0t-snd385t,'VariableNames',{'snd0t','snd385t','dif'});
        
        blocks1 = 0:sweep_time:d1_dur-sweep_time;
        blocks2 = sweep_time:sweep_time:d1_dur;
        
        % Function generator: 55 cycles per 20 sec sweep
        nperiods = arrayfun(@(n) sum(periods385(:,1)>blocks1(n) & periods385(:,1)<=blocks2(n)),...
            1:ns, 'UniformOutput', false);
        
        p385_tb = table(periods385(:,1),periods385(:,2),periods385(:,2)-periods385(:,1),...
            'VariableNames',{'RisingEdge','FallingEdge','Difference'});
        
        nperiods = cell2mat(nperiods);
        bf = cumsum(nperiods);
        bs = [1,bf(1:end-1)+1];
        
        [~,min_pts] = arrayfun(@(n) min(p385_tb.Difference(bs(n):bf(n))), 1:length(bs));
        min_pts = min_pts+bs-1;
        start_pts = min_pts+1; % Add 1: The next edge should be the largest pulse
        if start_pts(end)>=length(periods385)
            start_pts = start_pts(1:end-1);
        end
        
        % For each start_pt, keep selecting the next pulse if it is larger
        for i = 1:length(start_pts)
            idx = start_pts(i);
            w = p385_tb.Difference(idx);
            while p385_tb.Difference(idx+1)>w || p0_tb.Difference(idx)<0.4 
                start_pts(i) = idx+1;
                idx = idx+1;
                w = p385_tb.Difference(idx);
            end
        end
    
        % Remove duplicate points:
        [~, ~, ic] = unique(start_pts, 'stable');
        [~, first_occurrence] = unique(ic, 'first');
        all_indices = 1:numel(start_pts);
        dup = setdiff(all_indices, first_occurrence);
        start_pts(dup) = [];
        
        st_probe = p385_tb.RisingEdge(start_pts); % Leading pulse rising edge start times
        et_probe = p385_tb.RisingEdge(start_pts(2:end)); % End time: Before next leading pulse
        
        y1 = ones(length(st_probe),1)*0.9;
        y2 = ones(length(et_probe),1)*0.88;
        plot(st_probe,y1,'gx','MarkerSize',12)
        plot(et_probe,y2,'rx','MarkerSize',12)
    
        % Break here & uncomment code to fix edge detection errors:
        % Manually add last edge if missed
        % [x,y] = ginput(1);
        % x = 1601.593;
        % x = 1601.636;
        % st_probe = [st_probe;x];
        % et_probe = [et_probe;x];
        
        np = min(length(st_adc0),length(st_probe));
        st_adc0 = st_adc0(1:np);
        % et_adc0 = et_adc0(1:np);
        st_probe = st_probe(1:np);
        % et_probe = et_probe(1:np);
        
        dif_leading_edges = st_adc0-st_probe;
        dif_incr = [NaN;diff(dif_leading_edges)];
        
        tb1 = table(st_probe,st_adc0,dif_leading_edges,dif_incr,'VariableNames', ...
            {'LeadingEdgeProbe','LeadingEdgeADC0','Dif','Dif_Incr'});
        
        p_sweeptime = [NaN;diff(st_probe)];
        adc_sweeptime = [NaN;diff(st_adc0)];
        tb2 = table(st_probe,p_sweeptime,st_adc0,adc_sweeptime,'VariableNames',...
            {'LeadingEdgeProbe','P_SweepTime','LeadingEdgeADC0','ADC_SweepTime'});
        
        mean_rate = mean(dif_incr,'omitnan');
        rate1s = mean_rate/sweep_time;
    else
        st_probe = snd385t;

        ts_diff = st_adc0-st_probe;
        ts_table = table(st_adc0,st_probe,ts_diff,'VariableNames',...
            {'st_adc0','st_probe','diff'});

        dif = snd0t-snd385t;
        dif_dif = [NaN;diff(dif)];

        sn_table = table(snd0t,snd385t,dif,dif_dif,'VariableNames',...
            {'snd0t','snd385t','dif','dif_dif'});
    end
    
    max_adc0 = max(st_adc0);
    out_of_range =  periods1_vector > max_adc0;
    fprintf('%d values are out of range.\n', sum(out_of_range));

    % Align stimulus frames to Probe time:
    % Include linear extrapolation flag for values in periods1_vector >
    % st_adc0. Only use linear extrapolation if there is a linear relationship
    % between st_adc0 and st_probe
    Vq = interp1(st_adc0,st_probe,periods1_vector,'linear','extrap');
    
    % Add the total time of prior recordings
    Vq = Vq+d385_start_time; 
    
    Vq_periods = [Vq(1:2:end-1),Vq(2:2:end)];
    mean_dif = mean(diff(Vq(~isnan(Vq))));
    nPD_pulses = length(Vq_periods);

    if any(strcmp(protocol,{'Grating','Grating360'}))
        Vq = Vq_periods;
    end

end

function RFmapping_EB(ks_path,fn,trials,Vq,chan_config,channel_map,mode_sel,...
    nbins,save_figs,save_gauss,gauss_2d_data,gdi)
  
    close all
    save_dir = 'R:\Basic_Sciences\Phys\SenzaiLab\Shared\RFMapping_Analysis\RF_maps\';
    save_data = [save_dir,'Spike_Data\',fn];
    save_pdf = [save_dir,'Summary_Figures\',fn];
    save_gd = [save_dir,'GaussFit_Data\',fn];

    addpath(genpath("R:\Basic_Sciences\Phys\SenzaiLab\Yuta_Senzai\MatlabCodes\MATLAB\MyCodes"))
    addpath(genpath("D:\buzcode-master"))

    npulses = length(Vq);
    last_offset = Vq(end)+mean(diff(Vq));
    Vq_periods = [Vq,[Vq(2:end);last_offset]];
    
    nframes = height(trials);
    VisStim.periods = Vq_periods;
    VisStim.duration = VisStim.periods(:,2) - VisStim.periods(:,1);
    VisStim.PosX = [trials(1:nframes).Square_PositionX]';
    VisStim.PosY = [trials(1:nframes).Square_PositionY]';
    VisStim.Lum = [trials(1:nframes).Square_Luminance]';
    VisStim.SquareSize = trials(1).Square_Size;        
    
    sr = 30000;
    spike_times = readNPY([ks_path,'\spike_times.npy']);
    spike_times = double(spike_times)/sr;    
    spike_clusters = readNPY([ks_path,'\spike_clusters.npy']);
    unit_num = length(unique(spike_clusters));
    plot_wf = true;

    spike_templates = readNPY([ks_path,'\spike_templates.npy']);
    
    if ~isempty(find(spike_templates~=spike_clusters))
        % spike_templates will not match spike_clusters if clusters are
        % merged or split in Phy. This will result in plotting incorrect
        % spike waveforms
        display('Warning: spike_templates does not match spike_clusters')
        plot_wf = false;
    end

    % Find main channel for each unit
    templates = readNPY([ks_path,'\templates.npy']); % [nTemplates x nTimepoints x nChannels]
    
    % Manual kilosort: templates channels are ordered by channel map. e.g.
    % config 810to2250: template row 1 corresponds to channel 202, row 2 
    % is channel 204 etc.
    % Pipeline kilosort: templates channels are ordered 1 to 384. So row 1
    % is channel 0,..., row 203 is channel 202    
    
    % Pipeline kilosort: re-order template according to channel map
    if contains(ks_path,'ProbeA')
        templates = templates(:,:,channel_map+1);
    end  
    templates = permute(templates, [1 3 2]);      % Now: [nTemplates x nChannels x nTimepoints]    
           
    % Compute peak-to-peak amplitude per channel per template
    ptp = squeeze(max(templates, [], 3) - min(templates, [], 3));  % [nTemplates x nChannels]    
    % Find channel with max amplitude for each unit
    [~, max_chan_idx] = max(ptp, [], 2);  % 1-indexed

    chan_ids0 = channel_map(max_chan_idx);  
    % Add 1 for MATLAB 1-based indexing
    chan_ids1 = chan_ids0+1; 

    % Checked in Phy: channel ids are matched or close to channel ids in Phy.
    % ar = [good_units,chan_ids]; 

    % If file has been processed in Phy, analyze good units only
    labels_file = [ks_path,'\cluster_KSLabel.tsv'];
    if exist(labels_file, 'file') == 2
        unit_labels = readtable(labels_file, 'FileType', 'text', 'Delimiter', '\t');
        good_idx = (strcmp(unit_labels.KSLabel,'good'));
        good_units = unit_labels.cluster_id(good_idx); % 0-indexing, corresponding to Phy
        chan_ids0 = chan_ids0(good_idx);
        chan_ids1 = chan_ids1(good_idx);
        good_templates = templates(good_idx,:,:);
        good_max_chan_idx = max_chan_idx(good_idx);
        unit_num = length(good_units);
    end    
    ntp = size(good_templates,3);
    twf_time = (0:ntp-1)/sr;

    % Inhibitory RFs: [96,107,126]
    sel = 1:unit_num; % Select units for testing code
 
    SqDeg = VisStim.SquareSize;
    RFmap = cell(unit_num,1);    
    x_num = length(unique(VisStim.PosX));
    y_num = length(unique(VisStim.PosY));
    u_xy = length(unique([VisStim.PosX,VisStim.PosY,VisStim.Lum],'rows'));
    xy_ratio = x_num/y_num;

    % white squares and black squares each have [repnum] reps
    repnum = nframes/u_xy;  
   
    % Check if stimulus presentation is complete
    passes_dif = nframes-npulses;
    if passes_dif == 0
        display('Number of passes = number of pulses')
    elseif passes_dif < 0
        display('Number of pulses > number of passes')
    else
        new_repnum = floor(npulses/u_xy);
        sprintf('Number of pulses < number of frames, averaging over %d instead of %d repeats', ...
            new_repnum, repnum)
        repnum = new_repnum;
        nframes = u_xy*repnum;

        VisStim.periods = VisStim.periods(1:nframes,:);
        VisStim.duration = VisStim.duration(1:nframes);
        VisStim.PosX = VisStim.PosX(1:nframes);
        VisStim.PosY = VisStim.PosY(1:nframes);
        VisStim.Lum = VisStim.Lum(1:nframes);
    end
    ss = VisStim.SquareSize;

    if strcmp(mode_sel,'sum')
        sd = [save_data,'_RFmap_SumOfSpikes.mat'];
        sn = [save_pdf,'_Maps_SumOfSpikes.pdf'];
    elseif strcmp(mode_sel,'mean')
        sd = [save_data,'_RFmap_SpikeRate.mat'];
        sn = [save_pdf,'_Maps_SpikeRate.pdf'];
    end

    if exist(sd) == 2
        load(sd);
    else
        for k=1:unit_num
            RFmap{k}.ON.OnSet = zeros(y_num,x_num,nbins);
            RFmap{k}.OFF.OnSet = zeros(y_num,x_num,nbins);
        end

        for k=sel
            if exist('good_units', 'var')
                u = good_units(k);
            else
                u = k;
            end
            s = spike_times(spike_clusters == u);
            % Get mean spike sum or rate
            [sync,i] = Sync(s,Vq_periods(1,1),'durations', [-5 0]);
            [baseline,~] = SyncHist(sync,i,'durations', [-5 0],'nBins',1,'mode',mode_sel);
            if isempty(baseline)
                baseline = 0;
            end
            RFmap{k}.baseline = baseline;

            for x=1:x_num
                for y=1:y_num
                    curX = -(x_num-1)/2*SqDeg + (x-1)*SqDeg;
                    curY = -(y_num-1)/2*SqDeg + (y-1)*SqDeg;
                    IDon = VisStim.PosX==curX&VisStim.PosY==curY&VisStim.Lum==1;
                    IDoff= VisStim.PosX==curX&VisStim.PosY==curY&VisStim.Lum==0;
        
                    [sync,i] = Sync(s,VisStim.periods(IDon,1),'durations', [0 0.1]);
                    [hist,~] = SyncHist(sync,i,'durations', [0 0.1],'nBins',nbins,'mode',mode_sel);
                    if ~isempty(hist)
                        RFmap{k}.ON.OnSet(y,x,:) = hist;
                    end
        
                    [sync,i] = Sync(s,VisStim.periods(IDoff,1),'durations', [0 0.1]);
                    [hist,~] = SyncHist(sync,i,'durations', [0 0.1],'nBins',nbins,'mode',mode_sel);
                    if ~isempty(hist)
                        RFmap{k}.OFF.OnSet(y,x,:) = hist;
                    end
                end
            end 
            display(['done ' num2str(k) ' out of ' num2str(unit_num)]);
        end
        save(sd,'RFmap');
    end  
    xdeg = unique(VisStim.PosX);
    ydeg = unique(VisStim.PosY);

    % Make summary pdf
    if isempty(gauss_2d_data)

        if save_figs && exist(sn) == 2 
            delete(sn)
        end
    
        gauss_2d_data = cell(length(sel),1);
        for k=sel
            if exist('good_units', 'var')
                u = good_units(k);
                chan_idx = good_max_chan_idx(k);
            else
                u = k;
            end  
    
            % Get largest amplitude waveform and waveforms on neighbouring channels
            chans_idx = chan_idx-2:chan_idx+2;
            chans_idx = chans_idx(chans_idx > 0 & chans_idx < 385);
    
            twfs = arrayfun(@(c) squeeze(good_templates(k,c,:)), chans_idx,...
                'UniformOutput',false); % Template waveforms
            twfs = cell2mat(twfs);      
    
            fig=figure('Units','centimeters','Position',[23.336,4.3,17,17]);
            m = {RFmap{k}.ON.OnSet,RFmap{k}.OFF.OnSet};
            t = {'ON stim RF map','OFF stim RF map'};
            gauss_2d = cell(1,2);
            sb_pix = 30/ss; % Make scale bar 30°
            baseline = RFmap{k}.baseline;
            for ii = 1:2
                sp = subplot(2,2,ii);
                if strcmp(mode_sel,'sum')
                    rf = sum(m{ii},3)-baseline;
                elseif strcmp(mode_sel,'mean')
                    rf = mean(m{ii},3)-baseline;
                    mx_fr = max(rf(:));
                    min_fr = min(rf(:));                
                    if mx_fr>0 && mx_fr>abs(min_fr)
                        kc_sign = 1;
                    else
                        kc_sign = -1;
                    end
                end
                imagesc(xdeg,ydeg,rf); hold on 
                pos = sp.Position;           
                sp.Position(4) = pos(3)/xy_ratio;  
                cb = colorbar;
                cb_pos = cb.Position;
                cb.Position([1,3]) = [1.02*(pos(1)+pos(3)),cb_pos(3)*0.6]; 
                im = sp.Children; 
                axis off % axis on
                plot([sp.XLim(2),sp.XLim(2)-30],[sp.YLim(2),sp.YLim(2)],'r-',...
                    'LineWidth',1.2)
                title(t{ii});
                if ii == 1
                    text(sp.XLim(2)-15,sp.YLim(2)+5,'30°')
                elseif strcmp(mode_sel,'sum')
                    cb.Label.String = '# spikes';
                elseif strcmp(mode_sel,'mean')
                    cb.Label.String = 'imp/s';
                end
                % Compute 2-D gaussian fit
                % TO DO: 
                % 1. Adjust p0, LB, UB to better capture RFs on the edges
                % 2. Include inhibitory RFs (k: 96,107,126)
                % gauss_2d_model = @(p,x,y) ...
                %     p(1) * exp(-((( ( (x - p(2))*cos(p(7)) + (y - p(4))*sin(p(7)) ).^2 ) / (2*p(3)^2) + ...
                %       ((-(x - p(2))*sin(p(7)) + (y - p(4))*cos(p(7)) ).^2 ) / (2*p(5)^2)))) + p(6);
                gauss_2d_model = @(p,x,y) ...
                    p(1) * exp(-((( ( (x - p(2))*cos(p(6)) + (y - p(4))*sin(p(6)) ).^2 ) / (2*p(3)^2) + ...
                      ((-(x - p(2))*sin(p(6)) + (y - p(4))*cos(p(6)) ).^2 ) / (2*(p(3)*p(7))^2)))) + p(5);
    
                [x,y] = meshgrid(xdeg, ydeg);
                xy = [x(:),y(:)];
                gauss_fun = @(p) sum((gauss_2d_model(p, x(:), y(:)) - rf(:)).^2);
    
                % p(1): Amplitude
                % p(2): Center x-coordinate
                % p(3): Standard deviation along x (σₓ)
                % p(4): Center y-coordinate
                % p(5): Baseline spike sum/rate (free parameter, baseline already subtracted)
                % p(6): Gaussian orientation in radians
                % p(7): σᵧ/σₓ ratio, constrained to keep minor axis length at least 0.5 x major axis length
                % p(8): Standard deviation along y (σᵧ) determined during fit, added to pfit_cell
                kc = max(rf(:));
                cx = 0;
                sdx = 3*SqDeg;
                cy = 0;
                b = 0;
                theta = 0;
                r = 1;
                
                % Note making LB of kc negative reduces the accuracy of many
                % fits. Find a way to decide for individual cells whether to
                % use a negative kc LB.
                p0 = [kc, cx, sdx, cy, b, theta, r];
                LB = [0, min(xdeg)-2*SqDeg, SqDeg, min(ydeg)-2*SqDeg, 0, -pi/2, 0.6];
                UB = [3*kc, max(xdeg)+2*SqDeg, 0.5*max(ydeg), max(ydeg)+2*SqDeg, 0.25*kc, pi/2, 1.9];
    
                % Adjust parameters for inhibitory RF
                if kc_sign==-1
                    kc = min(rf(:));
                    p0(1) = kc;
                    LB([1,5]) = [3*kc,0.25*kc];
                    UB([1,5]) = [0,0];
                end
    
                % inequality constraint c(p) <= 0 → ensures min >= 0.5*max
                nonlcon = @(p) deal([], 0.5*max(p(3),p(5)) - min(p(3),p(5))); 
                % pfit = fmincon(gauss_fun, p0, [], [], [], [], LB, UB, nonlcon);           
                pfit = lsqcurvefit(@(p,xy) gauss_2d_model(p, xy(:,1), xy(:,2)), ...
                    p0, [x(:), y(:)], rf(:), LB, UB); 
    
                sdy = pfit(7)*pfit(3); % σᵧ
                txt = sprintf('k: %0.1f imp/s, σx: %0.1f°, σy: %0.1f°\nΘ: %0.1f°, sign: %d',...
                    pfit([1,3]),sdy,rad2deg(pfit(6)),kc_sign);
    
                gauss_2d{ii}.pfit = [pfit,sdy];
                gauss_2d{ii}.pfit_labels = {'k(imp/s)','cx(°)','σx(°)','cy(°)','baseline(imp/s)', ...
                  'theta(rad)','σy/σx','σy(°)'};
                gauss_2d{ii}.rf_fit = gauss_2d_model(pfit, x, y);           
                text(sp.XLim(1),1.5*sp.YLim(2),txt)
            end 
            colormap(gray) 
            ax = findall(fig,'Type','axes');
            cmx = max(arrayfun(@(n) ax(n).CLim(2), 1:length(ax)));
            cmin = min(arrayfun(@(n) ax(n).CLim(1), 1:length(ax)));
            arrayfun(@(n) set(ax(n),'CLim',[cmin,cmx]), 1:length(ax));
    
            ax_inorder = flipud(ax);
            t = linspace(0, 2*pi, 300);
            for jj = 1:2
                axes(ax_inorder(jj))
                pfit = gauss_2d{jj}.pfit;                    
                % Unrotated ellipse (axis-aligned)
                xe = pfit(3) * cos(t);
                ye = pfit(8) * sin(t);        
                % Rotate by theta
                theta = pfit(6);
                R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
                xy_rot = R * [xe; ye];
                % Translate to Gaussian center
                xe_rot = xy_rot(1,:) + pfit(2);
                ye_rot = xy_rot(2,:) + pfit(4);
                plot(xe_rot, ye_rot, 'b', 'LineWidth', 1);
            end      
            sgtitle(sprintf('Unit %d',u))
    
            % Plot template waveforms
            yshift = 7;
            if plot_wf
                ax1_pos = ax(1).Position;
                % wf_ax = axes('Position',ax1_pos.*[1,0.6,0.4,0.5]); 
                wf_ax = subplot(2,2,4);
                wf_ax_pos = wf_ax.Position;           
                for ww = 1:length(chans_idx)
                    if chans_idx(ww) == chan_idx
                        plot(twf_time,twfs(:,ww)+((ww-1)*yshift),'Color',rgb('DodgerBlue'))
                    else
                        plot(twf_time,twfs(:,ww)+((ww-1)*yshift),'Color',[.2,.2,.2]) 
                    end
                    hold on % cla
                end
                yl1 = wf_ax.YLim(1)-1;
                plot([0,0.001],[yl1,yl1],'k-')
                text(0.00025,yl1-2,'1 ms','FontSize',8), axis off
                wf_ax.Position = [wf_ax_pos(1:2).*[1.05,1.5], ax1_pos(3:4).*[0.4,1]];
                ty = 0:yshift:yshift*length(chans_idx);
                arrayfun(@(w) text(1.1*wf_ax.XLim(2),ty(w),num2str(channel_map(chans_idx(w)))),...
                    1:length(chans_idx))
            end
           
            chid = chan_ids0(k); % 0-indexed
            chid_pos = find(channel_map==chid); % Channel position in the channel map                
            PlotProbeConfig('RFmapping',chan_config,chid_pos)
    
            gauss_2d_data{find(sel==k)} = gauss_2d;
    
            if save_figs
                exportgraphics(fig, sn, 'Append', true);        
                close(fig);
            end
        end
        if save_gauss
            save(save_gd,'gauss_2d_data')
        end
    else
        load(gdi)
    end
    
    fig=figure('Units','centimeters','Position',[23.336,4.3,17,17]);
    colors = distinguishable_colors(length(sel)); 
    handles = gobjects(length(gauss_2d_data),1);   
    labels  = zeros(length(gauss_2d_data),1);  
    titles = {'ON stimulus','OFF stimulus'};
    for ss = 1:2
        sp = subplot(2,2,ss);
        set(sp,'XLim',xdeg([1,end]),'YLim',ydeg([1,end]))
        pos = sp.Position;           
        sp.Position(4) = pos(3)/xy_ratio;  
    
        for kk = 1:length(gauss_2d_data) 
            pfit = gauss_2d_data{kk}{ss}.pfit;
            t = linspace(0, 2*pi, 300);        
            % Unrotated ellipse (axis-aligned)
            xe = pfit(3) * cos(t);
            ye = pfit(8) * sin(t);        
            % Rotate by theta
            theta = pfit(6);
            R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
            xy_rot = R * [xe; ye];
            % Translate to Gaussian center
            xe_rot = xy_rot(1,:) + pfit(2);
            ye_rot = xy_rot(2,:) + pfit(4);
            % Invert ye: negative y values correspond to top of stim monitor  
            h = plot(xe_rot, -ye_rot, 'Color', colors(kk,:), 'LineWidth', 1); hold on
            handles(kk) = h;
            labels(kk) = kk;
        end
        set(sp,'XLim',xdeg([1,end]),'YLim',ydeg([1,end]))
        title(titles(ss))
    end
    valid = labels ~= 0;
    legend(handles(valid), string(labels(valid)),'Location','southoutside', ...
           'NumColumns', 4);
    ax = findall(gcf, 'type','axes');
    xlabel(ax(2),'deg')
    ylabel(ax(2),'deg')

    if save_figs
        exportgraphics(fig, sn, 'Append', true);        
        close(fig);
    end

end

function PlotProbeConfig(protocol,chan_config,chid_pos,varargin)

switch protocol
    case 'RFmapping'
        probe_plot = subplot(2,2,3); hold on  % cla
        probe_pos = probe_plot.Position;
        probe_plot.Position(3:4) = probe_pos(3:4).*[0.8,0.6]; 
    case 'Grating'
        probe_pos = varargin{1};
        probe_plot = axes('Position',probe_pos); hold on
end
num_shanks = 4;
shank_spacing = 200;   
shank_width = 70;      
shank_height = 10000;  
tip_height = 500;
ms1 = 1;
ms2 = 10;
       
shank_centers_x = linspace(-shank_spacing*1.5, shank_spacing*1.5, num_shanks);

for i = 1:num_shanks
    cx = shank_centers_x(i);  % x-center of the shank

    % Rectangle body (from y = tip_height to y = shank_height)
    x_rect = [cx - shank_width/2, cx + shank_width/2, cx + shank_width/2, cx - shank_width/2];
    y_rect = [tip_height, tip_height, shank_height, shank_height];
    fill(x_rect, y_rect, [0.8 0.8 0.8], 'EdgeColor', [0.8 0.8 0.8]);

    x_tip = [cx - shank_width/2, cx + shank_width/2, cx];
    y_tip = [tip_height, tip_height, 0];
    fill(x_tip, y_tip,[0.8 0.8 0.8],'EdgeColor',[0.8 0.8 0.8]);
end

xlim([-600, 600]);
ylim([-tip_height, shank_height + 500]);

% Add channels:
nc = 384;
nc1 = nc/num_shanks;
nc2 = nc1/2;
site_spacing = 20;           
column_offset = 16;          
% Y positions of sites: start from just above the tip
y_sites = tip_height + (0:nc-1) * site_spacing;  

for i = 1:num_shanks
    cx = shank_centers_x(i);                   
    x_left = cx - column_offset * ones(size(y_sites));
    x_right = cx + column_offset * ones(size(y_sites));        
    plot(x_left, y_sites, '.','Color',[0.5 0.5 0.5], 'MarkerSize', 0.8);
    plot(x_right, y_sites, '.','Color',[0.5 0.5 0.5], 'MarkerSize', 0.8);
end

% Show the channel configuration, and the channel this unit was
% recorded on
switch chan_config
    case 1 % 4 shanks, 0to1440        
    for i = 1:2
        y_sites = tip_height + (0:nc1-1) * site_spacing;
        cx = shank_centers_x(i);                   
        x_left = cx - column_offset * ones(1,nc1);     
        plot(x_left, y_sites, '.','Color',[0.1 0.1 0.1],'MarkerSize',ms1);
        cx = shank_centers_x(i+2);                   
        x_right = cx + column_offset * ones(size(y_sites));        
        plot(x_right, y_sites, '.','Color',[0.1 0.1 0.1],'MarkerSize',ms1);
    end
    if chid_pos <= nc1 % Shank 1: Channel position 0:96
        cx = shank_centers_x(1);                   
        x_left = cx - column_offset; 
        rel_pos = chid_pos;
    elseif chid_pos > nc1 && chid_pos <= nc1*2 % Shank 2: Channel position 97:192
        cx = shank_centers_x(2);                   
        x_left = cx - column_offset; 
        if chid_pos == nc1*2
            rel_pos = nc1;
        else
            rel_pos = mod(chid_pos,nc1);
        end
    elseif chid_pos > nc1*2 && chid_pos <= nc1*3 % Shank 3: Channel position 193:288
        cx = shank_centers_x(3);                   
        x_right = cx + column_offset;         
        if chid_pos == nc1*3
            rel_pos = nc1;
        else
            rel_pos = mod(chid_pos,nc1);
        end
     elseif chid_pos > nc1*3 && chid_pos <= nc1*num_shanks % Shank 4: Channel position 289:384
        cx = shank_centers_x(4);                   
        x_right = cx + column_offset;         
        if chid_pos == nc1*num_shanks
            rel_pos = nc1;
        else
            rel_pos = mod(chid_pos,nc1);
        end
    end
    plot(x_right, y_sites(rel_pos),'.','Color',rgb('dodgerblue'),'MarkerSize',ms2);

    case 2 % 4 shanks, 1440to2880
        % To do
    
    case 3 % 4 shanks, 2880to4320
        % To do
    
    case 4 % Shank4, 0to5760      
    cx = shank_centers_x(4);                   
    x_left = cx - column_offset * ones(size(y_sites));
    x_right = cx + column_offset * ones(size(y_sites)); 
    y_idx1 = logical([zeros(1,nc/2),ones(1,nc/2)]);
    y_idx2 = logical([ones(1,nc/2),zeros(1,nc/2)]);

    plot(x_left(y_idx1), y_sites(y_idx1),'.','Color',[0.1 0.1 0.1],'MarkerSize',ms1);
    plot(x_right(y_idx2), y_sites(y_idx2),'.','Color',[0.1 0.1 0.1],'MarkerSize',ms1);

    if chid_pos < (nc/2)+1
        plot(x_right(chid_pos), y_sites(chid_pos),'.','Color',rgb('dodgerblue'),'MarkerSize',ms2);
    else
        plot(x_left(chid_pos), y_sites(chid_pos),'.','Color',rgb('dodgerblue'),'MarkerSize',ms2);
    end

    case 5 % Shank1BankA
        % To do

    case 6 % 4 shanks, 150to1590   cla
    y0 = 150;
    for i = 1:2
        y_sites = tip_height + ((0:nc1-11)*site_spacing)+y0;
        y_sites96 = tip_height + ((0:nc1-1)*site_spacing)+y0;
        cx = shank_centers_x(i);                   
        x_left = cx - column_offset * ones(1,nc1-10);     
        plot(x_left, y_sites, '.','Color',[0.1 0.1 0.1],'MarkerSize',ms1);    

        cx = shank_centers_x(i+2);                   
        x_right = cx + column_offset * ones(size(y_sites));        
        plot(x_right, y_sites, '.','Color',[0.1 0.1 0.1],'MarkerSize',ms1); 

        % Top 10 channels on each shank switch columns
        y_sites2 = tip_height + ((nc1-10:nc1)*site_spacing)+y0;
        cx = shank_centers_x(i); 
        x_shiftR = cx + column_offset * ones(size(y_sites2));        
        plot(x_shiftR, y_sites2, '.','Color',[0.1 0.1 0.1],'MarkerSize',ms1); 
        
        cx = shank_centers_x(i+2);
        x_shiftL = cx - column_offset * ones(size(y_sites2));     
        plot(x_shiftL, y_sites2, '.','Color',[0.1 0.1 0.1],'MarkerSize',ms1);
    end

    if chid_pos <= nc1 % Shank 1: Channel position 0:96
        cx = shank_centers_x(1); 
        if chid_pos <= nc1-10                  
            xv = cx - column_offset; 
        else
            xv = cx + column_offset; 
        end
        rel_pos = chid_pos;

    elseif chid_pos > nc1 && chid_pos <= nc1*2 % Shank 2: Channel position 97:192
        cx = shank_centers_x(2);  
        if chid_pos <= (nc1*2)-10                 
            xv = cx - column_offset; 
        else
            xv = cx + column_offset; 
        end
        if chid_pos == nc1*2
            rel_pos = nc1;
        else
            rel_pos = mod(chid_pos,nc1);
        end

    elseif chid_pos > nc1*2 && chid_pos <= nc1*3 % Shank 3: Channel position 193:288
        cx = shank_centers_x(3); 
        if chid_pos <= (nc1*3)-10
            xv = cx + column_offset; 
        else
            xv = cx - column_offset; 
        end
        if chid_pos == nc1*3
            rel_pos = nc1;
        else
            rel_pos = mod(chid_pos,nc1);
        end
     
    elseif chid_pos > nc1*3 && chid_pos <= nc1*num_shanks % Shank 4: Channel position 289:384
        cx = shank_centers_x(4); 
        if chid_pos <= (nc1*num_shanks)-10
            xv = cx + column_offset; 
        else
            xv = cx - column_offset; 
        end
        if chid_pos == nc1*num_shanks
            rel_pos = nc1;
        else
            rel_pos = mod(chid_pos,nc1);
        end
    end
    plot(xv, y_sites96(rel_pos),'.','Color',rgb('dodgerblue'),'MarkerSize',ms2);

    case 7 % 4 shanks, 810 to 2250 (y position 54 to 150)
    y0 = 54; 
    yl = 149;
    y_sites = tip_height + ((y0:yl)*site_spacing);
      
    cx = shank_centers_x(1);                   
    xLL = cx - column_offset * ones(1,nc2-5);     
    plot(xLL, y_sites(1:nc2-5), '.','Color',[0.1 0.1 0.1],'MarkerSize',ms1);  
    xLU = cx + column_offset * ones(1,nc2+3);     
    plot(xLU, y_sites(nc2-2:end), '.','Color',[0.1 0.1 0.1],'MarkerSize',ms1);  
    
    cx = shank_centers_x(2);                   
    xRL = cx - column_offset * ones(1,nc2-5);        
    plot(xRL, y_sites(1:nc2-5), '.','Color',[0.1 0.1 0.1],'MarkerSize',ms1); 
    xRU = cx + column_offset * ones(1,nc2+3);        
    plot(xRU, y_sites(nc2-2:end), '.','Color',[0.1 0.1 0.1],'MarkerSize',ms1); 

    cx = shank_centers_x(3);                   
    xLL = cx + column_offset * ones(1,nc2-5);     
    plot(xLL, y_sites(1:nc2-5), '.','Color',[0.1 0.1 0.1],'MarkerSize',ms1);  
    xLU = cx - column_offset * ones(1,nc2+3);     
    plot(xLU, y_sites(nc2-2:end), '.','Color',[0.1 0.1 0.1],'MarkerSize',ms1);  
    
    cx = shank_centers_x(4);                   
    xRL = cx + column_offset * ones(1,nc2-5);        
    plot(xRL, y_sites(1:nc2-5), '.','Color',[0.1 0.1 0.1],'MarkerSize',ms1); 
    xRU = cx - column_offset * ones(1,nc2+3);        
    plot(xRU, y_sites(nc2-2:end), '.','Color',[0.1 0.1 0.1],'MarkerSize',ms1);   


    if chid_pos <= nc1 % Shank 1: Channel position 0:96
        cx = shank_centers_x(1); 
        if chid_pos <= nc2-5                 
            xv = cx - column_offset; 
        else
            xv = cx + column_offset; 
        end
        rel_pos = chid_pos;

    elseif chid_pos > nc1 && chid_pos <= nc1*2 % Shank 2: Channel position 97:192
        cx = shank_centers_x(2);  
        if chid_pos <= nc1+nc2-5                
            xv = cx - column_offset; 
        else
            xv = cx + column_offset; 
        end
        if chid_pos == nc1*2
            rel_pos = nc1;
        else
            rel_pos = mod(chid_pos,nc1);
        end

    elseif chid_pos > nc1*2 && chid_pos <= nc1*3 % Shank 3: Channel position 193:288
        cx = shank_centers_x(3); 
        if chid_pos <= (nc1*3)-(nc2-5)
            xv = cx + column_offset; 
        else
            xv = cx - column_offset; 
        end
        if chid_pos == nc1*3
            rel_pos = nc1;
        else
            rel_pos = mod(chid_pos,nc1);
        end
     
    elseif chid_pos > nc1*3 && chid_pos <= nc1*num_shanks % Shank 4: Channel position 289:384
        cx = shank_centers_x(4); 
        if chid_pos <= (nc1*num_shanks)-(nc2-5)
            xv = cx + column_offset; 
        else
            xv = cx - column_offset; 
        end
        if chid_pos == nc1*num_shanks
            rel_pos = nc1;
        else
            rel_pos = mod(chid_pos,nc1);
        end
    end
    plot(xv, y_sites(rel_pos),'.','Color',rgb('dodgerblue'),'MarkerSize',ms2);

end

text(shank_centers_x-50,ones(1,4).*probe_plot.YLim(2),{'s1','s2','s3','s4'})
text(gca().XLim(1),-2000,sprintf('Channel distance from tip: %d µm',...
    (y0*15)+((rel_pos-1)*15)),'FontSize',9')
axis off

end

function colors = distinguishable_colors(n)
% Generate n maximally distinct, vivid colors (no white, no gray)

    % Sample RGB cube
    M = 40;
    [r,g,b] = ndgrid(linspace(0,1,M));
    cand = [r(:), g(:), b(:)];

    % Compute brightness and saturation
    brightness = max(cand,[],2);
    sat = std(cand,0,2);   % saturation ~ channel variance (0 = gray/white)

    % Remove near-white OR near-gray/faint colors
    keep = (brightness < 0.9) & (sat > 0.05);

    cand = cand(keep,:);

    % Safety: must have candidates
    if size(cand,1) < n
        error('Not enough vivid colors; reduce n or loosen thresholds.');
    end

    % Greedy farthest-point sampling
    colors = zeros(n,3);
    colors(1,:) = cand(randi(size(cand,1)),:);  % random vivid seed

    for i = 2:n
        d = pdist2(cand, colors(1:i-1,:));
        score = min(d,[],2);
        [~,best] = max(score);
        colors(i,:) = cand(best,:);
    end
end

function filledMask = fill_short_gaps(logicalMask, maxGapSamples)
%FILL_SHORT_GAPS Fill short false gaps between true samples in a logical mask.
%
% filledMask = fill_short_gaps(logicalMask, maxGapSamples)
%
% Inputs
%   logicalMask   : logical (or numeric) vector/array; treated as logical
%   maxGapSamples : maximum gap length (in samples) to bridge (default = 600)
%
% Output
%   filledMask    : logical mask with short gaps filled

    % Linear indices where the mask is true
    filledMask = logicalMask;
    trueSampleIndices = find(filledMask);
    % Nothing to bridge if fewer than two true samples
    if numel(trueSampleIndices) < 2
        return;
    end
    % Gap length between consecutive true samples (minus 1)
    gapLengths = diff(trueSampleIndices) - 1;
    % Locations (in the diff array) where the gap is short enough to fill
    shortGapLocations = find((gapLengths > 0) & (gapLengths <= maxGapSamples));
    % Fill each short gap
    for locationIdx = 1:numel(shortGapLocations)
        gapLocation = shortGapLocations(locationIdx);
        firstFalseAfterTrue = trueSampleIndices(gapLocation) + 1;
        lastFalseBeforeNextTrue = trueSampleIndices(gapLocation + 1) - 1;
        if firstFalseAfterTrue <= lastFalseBeforeNextTrue
            filledMask(firstFalseAfterTrue:lastFalseBeforeNextTrue) = true;
        end
    end
end

function c = rgb(name)
%RGB Returns RGB triplet for a given color name string (case-insensitive).
% Usage:
%   c = rgb('SkyBlue');
%   c = rgb('LightGreen');
%   c = rgb('Lilac');
% Returns a 1x3 RGB vector in the range [0 1].

    if nargin == 0
        error('Please provide a color name.');
    end
    name = lower(string(name));

    % Define color dictionary
    cols = {
        % Blues
        {'navy', 'darkblue'},        [0 0 0.5];
        'blue',                      [0 0 1];
        'dodgerblue',                [0.12 0.56 1];
        'skyblue',                   [0.53 0.81 0.92];
        'lightblue',                 [0.68 0.85 0.9];
        'steelblue',                 [0.27 0.51 0.71];

        % Greens
        {'green','lime'},            [0 1 0];
        'forestgreen',               [0.13 0.55 0.13];
        'limegreen',                 [0.2 0.8 0.2];
        'lightgreen',                [0.56 0.93 0.56];
        'mediumseagreen',            [0.24 0.7 0.44];
        'springgreen',               [0 1 0.5];
        'chartreuse',                [0.5 1 0];

        % Reds
        'red',                       [1 0 0];
        'darkred',                   [0.55 0 0];
        'indianred',                 [0.8 0.36 0.36];
        'lightcoral',                [0.94 0.5 0.5];
        'salmon',                    [0.98 0.5 0.45];
        'tomato',                    [1 0.39 0.28];

        % Oranges
        'orange',                    [1 0.5 0];
        'darkorange',                [1 0.55 0];
        'coral',                     [1 0.5 0.31];
        'orangered',                 [1 0.27 0];

        % Yellows
        'yellow',                    [1 1 0];
        'gold',                      [1 0.84 0];
        'khaki',                     [0.94 0.9 0.55];
        'lightyellow',              [1 1 0.88];

        % Purples
        'purple',                    [0.5 0 0.5];
        'indigo',                    [0.29 0 0.51];
        'violet',                    [0.93 0.51 0.93];
        'mediumorchid',             [0.73 0.33 0.83];
        'plum',                      [0.87 0.63 0.87];

        % Pinks
        'pink',                      [1 0.75 0.8];
        'hotpink',                   [1 0.41 0.71];
        'deeppink',                  [1 0.08 0.58];
        'lightpink',                 [1 0.71 0.76];
        'palevioletred',             [0.86 0.44 0.58];

        % Browns
        'brown',                     [0.65 0.16 0.16];
        'sienna',                    [0.63 0.32 0.18];
        'saddlebrown',               [0.55 0.27 0.07];
        'chocolate',                 [0.82 0.41 0.12];
        'peru',                      [0.8 0.52 0.25];

        % Grays / Neutrals
        {'gray','grey'},             [0.5 0.5 0.5];
        'lightgray',                 [0.83 0.83 0.83];
        'darkgray',                  [0.66 0.66 0.66];
        'slategray',                 [0.44 0.5 0.56];
        'black',                     [0 0 0];
        'white',                     [1 1 1];

        % Turquoise / Teals
        {'teal','turquoise'},        [0 0.5 0.5];
        'mediumturquoise',           [0.28 0.82 0.8];
        'paleturquoise',             [0.69 0.93 0.93];

        % Magentas
        {'magenta','fuchsia'},       [1 0 1];
        'orchid',                    [0.85 0.44 0.84];
        'mediumvioletred',           [0.78 0.08 0.52];

        % Custom Aliases
        'lilac',                     [0.78 0.64 0.78];
    };

    % Search color
    c = [];
    for i = 1:size(cols, 1)
        keys = cols{i, 1};
        if iscell(keys)
            if any(strcmp(name, keys))
                c = cols{i,2};
                return;
            end
        elseif strcmp(name, keys)
            c = cols{i,2};
            return;
        end
    end

    error('Color "%s" not found.', name);
end
