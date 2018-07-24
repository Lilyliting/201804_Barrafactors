%% Initialize Factors
% Get All Factors: factorsXlsFile must be updated manually periodically
%取得因子名称并保存成mat文件

% First trade day
tradeStart = '2010-01-01'; 
trycnt = 0;
%firsttradeday = factorDateMap('countryFRetMtrx');
firsttradeday = tradeStart;
while trycnt <10
    try
        posSta = tradingDaysMap(firsttradeday);
        break;
    catch
    end
    firsttradeday = datestr(datenum(firsttradeday)+1, 'yyyy-mm-dd');
    trycnt=trycnt+1;
end

trdSta = find(Timeseries==datenum(firsttradeday,'yyyy-mm-dd'));
ts = size(tradingDaysLst,1);
dur = min(ts - posSta,size(Timeseries,1)-trdSta);
posEnd = posSta+dur-1;
espSta = 506;
espEnd = posSta+dur-1;
tradeEnd = datestr(Timeseries(espEnd),'yyyy-mm-dd');

fs = dir([stkdataFile,'*.mat']);
    stknum = size(fs,1);

[~,~,C] = xlsread(factorsXlsFile);
if exist(factorsMatFile,'file')
    load(factorsMatFile);
else
    factorLst = cell(0,0);
    factorMap = containers.Map();
    factorDateMap = containers.Map();
    ts = size(tradingDaysLst,1);
    for ii = 1 : size(C,1)
        factorId = C{ii,1};
        if ~factorMap.isKey(factorId)
            factorLst{end+1,1} = factorId;
            factorMap(factorId) = ii;
            factorDateMap([factorId 'RetMtrx']) = tradeStart;
            factorDateMap([factorId 'ExpMtrx']) = tradeStart;
            factorDateMap([factorId 'PtfRetMtrx']) = tradeStart;
            
            eval([factorLst{ii} 'RetMtrx=zeros(ts,1);']);
            eval([factorLst{ii} 'ExpMtrx=zeros(ts,stknum);']);
            
        end
    end
    save(factorsMatFile, 'factorLst', 'factorMap', 'factorDateMap');
end

%% Get Factor Returns1 (country)

if exist(factorReturnMatFile,'file')
    load(factorReturnMatFile);
else
    
    %     countryFRetMtrx(size(tradingDaysLst,1),1) = nan;
    %     sizeFRetMtrx(size(tradingDaysLst,1),1) = nan;
    %     momentumFRetMtrx(size(tradingDaysLst,1),1) = nan;
    %     residualVolFRetMtrx(size(tradingDaysLst,1),1) = nan;
    %     liquidityFRetMtrx(size(tradingDaysLst,1),1) = nan;
    %     stkUniqueFRetMtrx(size(tradingDaysLst,1),size(stkCodeLst,1)) = nan;
    
    % countryF
    disp([datestr(now) ': countryFRetMtrx']);
    factorDateMap('countryFRetMtrx') = tradingDaysLst{end,1};
%     countryFRetMtrx(posSta:posEnd,1) = nansum((stkEVMtrx(espSta:espEnd,:) ./ repmat(nansum(stkEVMtrx(espSta:espEnd,:),2), 1, size(stkEVMtrx,2))) .*...
%         (stkAdjcloseMtrx(espSta:espEnd,:) ./ stkAdjcloseMtrx((espSta-1):(espEnd-1),:) - 1), 2); 
    countryFRetMtrx(2:posEnd,1) = nansum((stkEVMtrx(2:espEnd,:) ./ repmat(nansum(stkEVMtrx(2:espEnd,:),2), 1, size(stkEVMtrx,2))) .*...
        (stkAdjcloseMtrx(2:espEnd,:) ./ stkAdjcloseMtrx(1:(espEnd-1),:) - 1), 2);
end

% if ih ic
% 没有

%% Calculate Each Factor Exposure
% factorExposureMatFile
if exist(factorExposureMatFile,'file')
    load(factorExposureMatFile);
else
%     sizeFExpMtrx = zeros(0,0);
%     momentumFExpMtrx = zeros(0,0);
%     alphaFExpMtrx = zeros(0,0);
%     betaFExpMtrx = zeros(0,0);
%     residualVolFExpMtrx = zeros(0,0);
%     liquidityFExpMtrx = zeros(0,0);
%     
%     
%     nowDate = datestr(now, 'yyyy-mm-dd');
%     sizeFExpMtrx(size(tradingDaysLst,1),size(stkCodeLst,1)) = nan;
%     momentumFExpMtrx(size(tradingDaysLst,1),size(stkCodeLst,1)) = nan;
%     alphaFMtrx(size(tradingDaysLst,1),size(stkCodeLst,1)) = nan;
%     alphaFExpMtrx(size(tradingDaysLst,1),size(stkCodeLst,1)) = nan;
%     betaFMtrx(size(tradingDaysLst,1),size(stkCodeLst,1)) = nan;
%     betaFExpMtrx(size(tradingDaysLst,1),size(stkCodeLst,1)) = nan;
%     residualVolFExpMtrx(size(tradingDaysLst,1),size(stkCodeLst,1)) = nan;
%     liquidityFExpMtrx(size(tradingDaysLst,1),size(stkCodeLst,1)) = nan;
    
    disp([datestr(now) ': Factor Exposure-------------------------'])
    % sizeF
    disp([datestr(now) ': sizeFExpMtrx']);
    factorDateMap('sizeFExpMtrx') = tradeEnd;
    sizeFExpMtrx(espSta:espEnd,:) = log(stkEVMtrx(espSta:espEnd,:));
    sizeFExpMtrx(espSta:espEnd,:) = Standardize(sizeFExpMtrx(espSta:espEnd,:));
    
    % momentumF
    disp([datestr(now) ': momentumFExpMtrx']);
    factorDateMap('momentumFExpMtrx') = tradingDaysLst{end,1};
    T = 504;
    hl = 126;
    L = 21;
    wExp = ExponentialWeight(T, hl);
    for ii = espSta:espEnd
        rt = stkAdjcloseMtrx(ii-T:ii-L,:) ./ stkAdjcloseMtrx(ii-T-1:ii-L-1,:) - 1;
        rft = repmat(idxCloseMtrx(ii-T:ii-L,4) / 25200, 1, size(stkAdjcloseMtrx,2));
        tmp = sum(repmat(wExp(1:end-L), 1, size(stkAdjcloseMtrx,2)) .* (log(1+rt) - log(1+rft)),1);
        momentumFExpMtrx(ii,:) = tmp;
    end
    momentumFExpMtrx(espSta:espEnd,:) = Standardize(momentumFExpMtrx(espSta:espEnd,:));
    
    % alphaFExpMtrx & betaFExpMtrx
    disp([datestr(now) ': alphaFExpMtrx & betaFExpMtrx']);
    factorDateMap('alphaFExpMtrx') = tradingDaysLst{end,1};
    factorDateMap('betaFExpMtrx') = tradingDaysLst{end,1};
    T = 252;
    hl = 63;
    wExp = ExponentialWeight(T, hl);
    for ii = espSta:espEnd
        rt = stkAdjcloseMtrx(ii-T:ii,:) ./ stkAdjcloseMtrx(ii-T-1:ii-1,:) - 1;
        rft = repmat(idxCloseMtrx(ii-T:ii,4) / 25200, 1, size(stkAdjcloseMtrx,2));
        Rt = repmat(countryFRetMtrx(ii-T:ii,1), 1, size(stkAdjcloseMtrx,2)) - rft;
        xt = Rt;
        yt = rt - rft;
        
        [alpha, beta] = LinearFitMtrx(xt,yt,wExp);
        alphaFMtrx(ii,:) = alpha;
        alphaFExpMtrx(ii,:) = (alphaFMtrx(ii,:) - nanmean(alphaFMtrx(ii,:))) / nanstd(alphaFMtrx(ii,:));
        betaFMtrx(ii,:) = beta;
        betaFExpMtrx(ii,:) = (betaFMtrx(ii,:) - nanmean(betaFMtrx(ii,:))) / nanstd(betaFMtrx(ii,:));
    end
    
    % residualVolF
    disp([datestr(now) ': residualVolFExpMtrx']);
    factorDateMap('residualVolFExpMtrx') = tradingDaysLst{end,1};
    T = 252;
    hl = 42;
    wExpDASTD = ExponentialWeight(T, hl);
    T = 252;
    hl = 63;
    wExpHSIGMA = ExponentialWeight(T, hl);
    for ii = espSta:espEnd
        %% DASTD
        rt = stkAdjcloseMtrx(ii-T:ii,:) ./ stkAdjcloseMtrx(ii-T-1:ii-1,:) - 1;
        rft = repmat(idxCloseMtrx(ii-T:ii,4) / 25200, 1, size(stkAdjcloseMtrx,2));
        wExpDASTDExpand = repmat(wExpDASTD,1,size(stkAdjcloseMtrx,2));
        DASTD = sqrt(sum(wExpDASTDExpand .* (rt - rft).^2, 1) ./ sum(wExpDASTDExpand, 1));
        DASTD = Standardize(DASTD);
        
        %% CMRA
        rt = flipud(stkAdjcloseMtrx((ii-11*21):21:ii,:) ./ stkAdjcloseMtrx((ii-12*21):21:(ii-21),:) - 1);
        rft = flipud(repmat((idxCloseMtrx((ii-11*21):21:ii,4) + idxCloseMtrx((ii-12*21):21:(ii-21),4)) / 2 / 2100, 1, size(stkAdjcloseMtrx,2)));
        ZT = cumsum(log(1+rt) - log(1+rft),1);
        CMRA = log(abs(1+max(ZT,[],1))) - log(abs(1+min(ZT,[],1)));     % abs() for min(ZT) < -1
        CMRA = Standardize(CMRA);
        
        %% HSIGMA
        rt = stkAdjcloseMtrx(ii-T:ii,:) ./ stkAdjcloseMtrx(ii-T-1:ii-1,:) - 1;
        rft = repmat(idxCloseMtrx(ii-T:ii,4) / 25200, 1, size(stkAdjcloseMtrx,2));
        Rt = repmat(countryFRetMtrx(ii-T:ii,1), 1, size(stkAdjcloseMtrx,2)) - rft;
        et = (rt - rft) - repmat(alphaFMtrx(ii,:),size(rt,1),1) - repmat(betaFMtrx(ii,:),size(rt,1),1) .* Rt;
        wExpHSIGMAExpand = repmat(wExpHSIGMA,1,size(stkAdjcloseMtrx,2));
        HSIGMA = sqrt(sum(wExpHSIGMAExpand .* et.^2, 1) ./ sum(wExpHSIGMA, 1));
        HSIGMA = Standardize(HSIGMA);
        
        %% residualVolF
        residualVolFExpMtrx(ii,:) = 0.74 * DASTD + 0.16 * CMRA + 0.10 * HSIGMA;
    end
    
    % liquidityF
    disp([datestr(now) ': liquidityFExpMtrx']);
    factorDateMap('liquidityFExpMtrx') = tradingDaysLst{end,1};
    for ii = espSta:espEnd
        STOM_t = zeros(12,size(stkTurnOverMtrx,2));
        for jj = 1 : 12
            STOM_t(jj,:) = log(sum(stkTurnOverMtrx((ii-(jj-1)*21-20):(ii-(jj-1)*21),:),1)/100+1e-3);
        end
        
        %% STOM
        STOM = STOM_t(1,:);
                STOM = Standardize(STOM);
        
        %% STOQ
        STOQ = log(mean(exp(STOM_t(1:3,:)),1));
         STOQ = Standardize(STOQ);
        
        %% STOA
        STOA = log(mean(exp(STOM_t(1:12,:)),1));
        STOA = Standardize(STOA);
        
        %% liquidityF
        liquidityFExpMtrx(ii,:) = 0.35 * STOM + 0.35 * STOQ + 0.30 * STOA;
    end
    
    PEFExpMtrx(espSta:espEnd,:) = Standardize(stkPEMtrx(espSta:espEnd,:));  %8
    
    PETTMFExpMtrx(espSta:espEnd,:) = Standardize(stkPE_TTMMtrx(espSta:espEnd,:));
    PBFExpMtrx(espSta:espEnd,:) = Standardize(stkPBMtrx(espSta:espEnd,:));
    PSFExpMtrx(espSta:espEnd,:) = Standardize(stkPSMtrx(espSta:espEnd,:));
    PSTTMFExpMtrx(espSta:espEnd,:) = Standardize(stkPS_TTMMtrx(espSta:espEnd,:)); %12
    
    PCFNFExpMtrx(espSta:espEnd,:) = Standardize(stkPCF_NCFMtrx(espSta:espEnd,:));
    PCFNTTMFExpMtrx(espSta:espEnd,:) = Standardize(stkPCF_NCFTTMMtrx(espSta:espEnd,:));
    PCFOFExpMtrx(espSta:espEnd,:) = Standardize(stkPCF_OCFMtrx(espSta:espEnd,:)); %15
    PCFOTTMFExpMtrx(espSta:espEnd,:) = Standardize(stkPCF_OCFTTMMtrx(espSta:espEnd,:));
    
    DPFExpMtrx(espSta:espEnd,:) = Standardize(stkDivPerShareMtrx(espSta:espEnd,:));
    revgrowthFExpMtrx(espSta:espEnd,:) = Standardize(stkRevenueMtrx(espSta:espEnd,:)./stkRevenueMtrx((espSta-252):(espEnd-252),:)-1);
    OCFLYRFExpMtrx(espSta:espEnd,:) = Standardize(stkOCF_LYRMtrx(espSta:espEnd,:)); %19
%     
%     netprofitgrowthFExpMtrx(espSta:espEnd,:) = Standardize(stkPEMtrx(espSta:espEnd,:));
    EPSFExpMtrx(espSta:espEnd,:) = Standardize(stkEPSMtrx(espSta:espEnd,:));
    dilutedEPSFExpMtrx(espSta:espEnd,:) = Standardize(stkDiluted_EPSMtrx(espSta:espEnd,:)); %22
    
%     ROEgrowthFExpMtrx(espSta:espEnd,:) = Standardize(stkPEMtrx(espSta:espEnd,:));
    grossprofitFExpMtrx(espSta:espEnd,:) = Standardize(stkGrossProfitMarginMtrx(espSta:espEnd,:));
    netprofitFExpMtrx(espSta:espEnd,:) = Standardize(stkNetProfitMarginMtrx(espSta:espEnd,:));
    
    ROEFExpMtrx(espSta:espEnd,:) = Standardize(stkROEMtrx(espSta:espEnd,:)); %26
    ROAFExpMtrx(espSta:espEnd,:) = Standardize(stkROAMtrx(espSta:espEnd,:)); 
    MACD_DIFFFExpMtrx(espSta:espEnd,:) = Standardize(stkMACD_DIFFMtrx(espSta:espEnd,:));
    MACD_DEAFExpMtrx(espSta:espEnd,:) = Standardize(stkMACD_DEAMtrx(espSta:espEnd,:));
    
    debtFExpMtrx(espSta:espEnd,:) = Standardize(stkDebtMtrx(espSta:espEnd,:)); %30 
    cashRatioFExpMtrx(espSta:espEnd,:) = Standardize(stkCashRatioMtrx(espSta:espEnd,:));
    quickRatioFExpMtrx(espSta:espEnd,:) = Standardize(stkQuickRatioMtrx(espSta:espEnd,:));
    currentRatioFExpMtrx(espSta:espEnd,:) = Standardize(stkCurrentRatioMtrx(espSta:espEnd,:));
    
    
    % save factorExposure to mat file
    save(factorExposureMatFile, 'sizeFExpMtrx', 'momentumFExpMtrx', ...
        'alphaFMtrx', 'alphaFExpMtrx', 'betaFMtrx', 'betaFExpMtrx', ...
        'residualVolFExpMtrx', 'liquidityFExpMtrx','PEFExpMtrx',...
        'PETTMFExpMtrx','PBFExpMtrx','PSFExpMtrx','PSTTMFExpMtrx',...
        'PCFNFExpMtrx','PCFNTTMFExpMtrx','PCFOFExpMtrx','PCFOTTMFExpMtrx',...
        'DPFExpMtrx','revgrowthFExpMtrx','OCFLYRFMtrx',...
        'netprofitgrowthFExpMtrx','EPSFExpMtrx','dilutedEPSFExpMtrx',...
        'ROEgrowthFExpMtrx','grossprofitFExpMtrx','netprofitFExpMtrx',...
        'ROEFExpMtrx','ROAFExpMtrx','MACD_DIFFFExpMtrx','MACD_DEAFExpMtrx',...
        'debtFExpMtrx','cashRatioFExpMtrx','quickRatioFExpMtrx','currentRatioFExpMtrx');
end


%% Get Factor Returns2 (others)
disp([datestr(now) ': Factor Returns--------------------------']);
posSta = tradingDaysMap(firsttradeday);
posEnd = size(tradingDaysLst,1);
factorDateMap('sizeFRetMtrx') = tradingDaysLst{end,1};
factorDateMap('momentumFRetMtrx') = tradingDaysLst{end,1};
factorDateMap('residualVolFRetMtrx') = tradingDaysLst{end,1};
factorDateMap('liquidityFRetMtrx') = tradingDaysLst{end,1};
factorDateMap('stkUniqueFRetMtrx') = tradingDaysLst{end,1};
for ii = posSta : posEnd
    X = [sizeFExpMtrx(ii,:);momentumFExpMtrx(ii,:);residualVolFExpMtrx(ii,:);liquidityFExpMtrx(ii,:);];
    rt = stkAdjcloseMtrx(ii,:) ./ stkAdjcloseMtrx(ii-1,:) - 1;
    rft = repmat(idxCloseMtrx(ii,4) / 25200, 1, size(stkAdjcloseMtrx,2));
    fc = repmat(countryFRetMtrx(ii,1), 1, size(stkAdjcloseMtrx,2));
    Y = rt - rft - fc;
    W = ones(1, size(stkAdjcloseMtrx,2));
    
    [F, Err] = MultiFactorsRegression(X, Y, W, 'real');
    
    sizeFRetMtrx(ii,1) = F(1);
    momentumFRetMtrx(ii,1) = F(2);
    residualVolFRetMtrx(ii,1) = F(3);
    liquidityFRetMtrx(ii,1) = F(4);
    stkUniqueFRetMtrx(ii,:) = Err;
end
stkUniqueFRetMtrx = real(stkUniqueFRetMtrx);
save(factorReturnMatFile,'countryFRetMtrx','ifBasisRetFRetMtrx','ihBasisRetFRetMtrx','icBasisRetFRetMtrx','sizeFRetMtrx','momentumFRetMtrx','residualVolFRetMtrx','liquidityFRetMtrx','stkUniqueFRetMtrx');








