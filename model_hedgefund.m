%% Get Portfolio Names
tmp = dir([portfolioXlsPath '投顾业绩（*）.xlsx']); %提取投顾业绩的文件名称
portfolioXlsFiles = cell(size(tmp,1),1); %建立投顾业绩文件名下的列向量，向量行数取决于投顾业绩家数
for ii = 1 : size(tmp,1)
    portfolioXlsFiles{ii,1} = tmp(ii).name; %建立一系列投顾业绩文件
end
tmp = sort(portfolioXlsFiles); %排序投顾业绩文件
portfolioXlsFile = [portfolioXlsPath tmp{end}]; %重新定义最后一个portfolioXlsFile文件


%% Get Portfolio Factor Exposures
factorsLst = {'countryFRetMtrx';'ifBasisRetFRetMtrx';'ihBasisRetFRetMtrx';'icBasisRetFRetMtrx';'sizeFRetMtrx';'momentumFRetMtrx';'residualVolFRetMtrx';'liquidityFRetMtrx';};
factorRetMtrx = [countryFRetMtrx ifBasisRetFRetMtrx ihBasisRetFRetMtrx icBasisRetFRetMtrx sizeFRetMtrx momentumFRetMtrx residualVolFRetMtrx liquidityFRetMtrx];

tradingDaysNum = datenum(tradingDaysLst);
[~,~,C] = xlsread(portfolioXlsFile, '概览');
if exist(portfolioMatFile,'file')
    load(portfolioMatFile);
else
    portfolioLst = cell(0);
    portfolioMap = containers.Map();
    portfolioDateMap = containers.Map();
    portfolioNavMap = containers.Map();
    portfolioFactorRetMap = containers.Map();
    portfolioFactorExpMap = containers.Map();
    portfolioRetDecompMap = containers.Map();
    
    T = 500;
    hl = 63;
    wReg = ExponentialWeight(T, hl);
    for ii = 2 : size(C,1)
        ptfId = C{ii,2};
        disp([datestr(now) ': ' ptfId]);
        try
            [~,~,navInfo] = xlsread(portfolioXlsFile, ptfId);
        catch
            continue;
        end
        
        if ~portfolioMap.isKey(ptfId)
            portfolioLst{end+1,1} = ptfId;
            portfolioMap(ptfId) = size(portfolioLst,1);
        end
        
        dateAllNum = datenum(navInfo(2:end,1));
        navAll = cell2mat(navInfo(2:end,2));
        portfolioDateMap(ptfId) = datestr(dateAllNum(end,1),'yyyy-mm-dd');
        
        oneNav = zeros(size(dateAllNum,1),4) * nan;
        oneFctRet = zeros(size(dateAllNum,1),size(factorsLst,1)+1-3) * nan;
        oneFctExp = zeros(size(dateAllNum,1),size(factorsLst,1)-3) * nan;
        for jj = 1 : size(dateAllNum,1)
            oneNav(jj,1) = find(tradingDaysNum >= dateAllNum(jj,1), 1, 'first'); % trade day
            oneNav(jj,4) = tradingDaysNum(oneNav(jj,1));
            if jj ~= 1
                oneFctRet(jj,1:end-1) = prod(1 + factorRetMtrx(oneNav(jj-1,1)+1:oneNav(jj,1),[1,5:8]), 1) - 1;
            end
        end
        oneNav(:,2) = navAll;
        oneNav(2:end,3) = oneNav(2:end,2) ./ oneNav(1:end-1,2) - 1;
        
        if size(dateAllNum,1) > size(factorsLst,1)+5
            for jj = size(factorsLst,1)+5 : size(dateAllNum,1)
                posSta = find(dateAllNum >= tradingDaysNum(oneNav(jj,1) - 500), 1, 'first');
                if posSta == 1
                    posSta = 2;
                end
                posEnd = jj;
                
                tmp = oneNav(posSta:posEnd,1) - oneNav(posEnd,1) + T + 1;
                wTmp = wReg(tmp);
                X = oneFctRet(posSta:posEnd,1:end-1);
                Y = oneNav(posSta:posEnd,3);
                
                [F,Err] = MultiFactorsRegression(X, Y, wTmp, 'real');
                oneFctExp(jj,:) = F;
                oneFctRet(jj,end) = real(Err(end));
            end
        end
        
        portfolioNavMap(ptfId) = oneNav;
        portfolioFactorRetMap(ptfId) = oneFctRet;
        portfolioFactorExpMap(ptfId) = oneFctExp;
        portfolioRetDecompMap(ptfId) = [oneFctExp.*oneFctRet(:,1:end-1) oneFctRet(:,end)];
    end
    
    save(portfolioMatFile, 'portfolioLst', 'portfolioMap', 'portfolioDateMap', 'portfolioNavMap', 'portfolioFactorRetMap', 'portfolioFactorExpMap', 'portfolioRetDecompMap');
end

%% hedge fund

factorsLst = {'countryFRetMtrx';'sizeFRetMtrx';'momentumFRetMtrx';'residualVolFRetMtrx';'liquidityFRetMtrx';}; 
for ii=2:size(C,1)
    ptfId = C{ii,2};
    %ptfId = '牟合资产信合1号'; 
    oneFctRet = portfolioFactorRetMap(ptfId); 
    oneFctExp = portfolioFactorExpMap(ptfId); 
    oneRetDmp = portfolioRetDecompMap(ptfId); 
    oneNav = portfolioNavMap(ptfId);
    clf;
    figure(1);
    subplot(3,1,1);
    date = oneNav(:,4);
    accumret = cumprod(nansum(oneRetDmp,2)+1)-1; % 累计收益率
    plot(date, accumret);
    hold on; 
    plot(date(2:end), cumprod(oneNav(2:end,3)+1)-1);
    title(ptfId)
    legend({'估计累积净值', '实际累积净值'},'Location','northwest')
    dateaxis('x',1);
    
    subplot(3,1,2);
    plot(date, oneFctExp) 
    legend(factorsLst,'Location','northwest'); 
    dateaxis('x',1);
    
    retDmp = [factorsLst;{'ResidualErr'}]; 
    oneRetDmp(isnan(oneRetDmp)) = 0; 

    subplot(3,1,3);
    plot(date, cumprod(oneRetDmp+1)-1); 
    dateaxis('x',1);
    legend(retDmp,'Location','northwest');
    saveas(gcf,['.\figures\' ptfId '.jpg']);
    pause(1)
end
