%% CCG in SC -- CCG.m in buzcode

% CCG structure: wake, OpenField, REM, NREM,

cd([basepath '\' metaname]);
load([metaname '-HomeCage-OpenField-Periods']);

for ff=1:length(fbasename_pool)
    fbasename = fbasename_pool{ff};
    cd([basepath '\' fbasename]);

    load([fbasename '-UnitFeature.mat']);
    load([fbasename '.SleepState.states.mat']);

    SetCurrentSession('basename',fbasename);
    unit_pool_orig = GetUnits;
    unit_idx = unit_pool_orig(:,1)==SCshank{ff}; % SC
    unit_pool = unit_pool_orig(unit_idx,:);

    unit_num = size(unit_pool,1);

    s = GetSpikeTimes(unit_pool,'output','numbered');

    DepthVec = UnitFeature.Depth(unit_idx);
    dist_mtrx = zeros(unit_num,unit_num);

    for ii=1:unit_num
        for jj=1:unit_num
            dist_mtrx(ii,jj)=DepthVec(ii) - DepthVec(jj);
        end
    end
  
    HCper = HomeCagePeriods{ff};
    OFper = OpenFieldPeriods{ff};
    numHC = size(HCper,1);
    numOF = size(OFper,1);

    HC_periods = [];

    for hchc=1:numHC
        t_hc = (HCper(hchc,1):HCper(hchc,2))';
        [HC_wakeIdx,~,~] = InIntervals(t_hc,SleepState.ints.WAKEstate);
        [HC_periods_cur,~] = Threshold([t_hc, HC_wakeIdx],'>',0);
        HC_periods = [HC_periods; HC_periods_cur];
    end

   

    OF_periods = [];

    for ofof=1:numOF
        t_of = (OFper(ofof,1):OFper(ofof,2))';
        [OF_wakeIdx,~,~] = InIntervals(t_of,SleepState.ints.WAKEstate);
        [OF_periods_cur,~] = Threshold([t_of, OF_wakeIdx],'>',0);
        OF_periods = [OF_periods; OF_periods_cur];
    end
   
    CCG_SC = cell(4,1);
    fig=figure;

    for j=1:4
        if j==1
            curperiod = round(OF_periods);
            savename = 'WAKE_OF';
        elseif j==2
            curperiod = SleepState.ints.REMstate ;
            savename = 'REM';
        elseif j==3
            curperiod = round(HC_periods) ;
            savename = 'WAKE_HC';
        elseif j==4
            curperiod = SleepState.ints.NREMstate ;
            savename = 'NREM';
        end

        s_res = Restrict(s,[curperiod]);
        ts_res = s_res(:,1);
        gs_res = s_res(:,2);
        [ccg,t_ccg] = CCG(ts_res,gs_res,'binSize',0.02,'duration',8);
        binnum = size(ccg,1);
        normCCG_ppln=[];
        dist_ppln=[];
        preID_ppln=[];
        postID_ppln=[];

        for k=1:unit_num

            CCGmtrx = squeeze(ccg(:,k,1:(k-1)))';
            dist_cur= dist_mtrx(k,1:(k-1))';
            normCCGmtrx = CCGmtrx./repmat(median(CCGmtrx,2),[1 binnum]);
            preID = ones(k-1,1)*k;
            postID = (1:(k-1))';
            normCCG_ppln = [normCCG_ppln; normCCGmtrx];
            dist_ppln = [dist_ppln; dist_cur];
            preID_ppln = [preID_ppln; preID];
            postID_ppln = [postID_ppln; postID];
            % peakval_mtrx(:,k)= mean(normCCGmtrx(:,200:202),2);
        end

        %         [dist_rad_sorted, sortedIdx2] = sort(dist_rad_ppln);

        %         normCCG_sorted = normCCG_ppln(sortedIdx2,:);

        CCG_SC{j}.ccg = ccg;
        CCG_SC{j}.normCCG = normCCG_ppln;
        CCG_SC{j}.dist = dist_ppln;
        CCG_SC{j}.preID = preID_ppln;

        CCG_SC{j}.postID = postID_ppln;

        CCG_SC{j}.t_ccg = t_ccg;

        CCG_SC{j}.state = savename;

        CCG_SC{j}.periods = curperiod;

        CCG_SC{j}.unitID = unit_idx;

       

        subplot(1,4,j);
        PlotColorMap(normCCG_ppln(abs(dist_ppln)>1,:),'x',t_ccg);set(gca, 'YDir','reverse');
        title(savename);

       

    end

   

    save([fbasename '--CCG_SC_REM_NREM_WAKE-HC-OF'],'CCG_SC');

   % SaveFigPngPDFSvg(fig,[fbasename '--CCG-SC_REM_NREM_WAKE-HC-OF']);

   
end

 

%% add CCGppln

for ff=1:length(fbasename_pool)

    fbasename = fbasename_pool{ff};

    cd([basepath '\' fbasename]);

    load([fbasename '--CCG_SC_REM_NREM_WAKE-HC-OF']);

   

    for j=1:4

        ccg = CCG_SC{j}.ccg;

        unit_num = sum(CCG_SC{j}.unitID==1);

        CCG_ppln=[];

        for k=1:unit_num

            CCGmtrx = squeeze(ccg(:,k,1:(k-1)))';

            CCG_ppln = [CCG_ppln; CCGmtrx];

        end

        CCG_SC{j}.CCGppln = CCG_ppln;

    end

    save([fbasename '--CCG_SC_REM_NREM_WAKE-HC-OF'],'CCG_SC');

end