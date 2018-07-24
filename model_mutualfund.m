%% Get Portfolio Factor Exposures

MutualfundMatFile = '.\mat\mutualfund.mat';

factorsLst = {'countryFRetMtrx';'ifBasisRetFRetMtrx';'ihBasisRetFRetMtrx';'icBasisRetFRetMtrx';'sizeFRetMtrx';'momentumFRetMtrx';'residualVolFRetMtrx';'liquidityFRetMtrx';};
factorRetMtrx = [countryFRetMtrx ifBasisRetFRetMtrx ihBasisRetFRetMtrx icBasisRetFRetMtrx sizeFRetMtrx momentumFRetMtrx residualVolFRetMtrx liquidityFRetMtrx];

tradingDaysNum = datenum(tradingDaysLst);
C = MutualfundCodeLst;

if exist(MutualfundMatFile,'file')
    load(MutualfundMatFile);
else
    MutualfundLst = cell(0);
    MutualfundMap = containers.Map();
    MutualfundDateMap = containers.Map();
    MutualfundNavMap = containers.Map();
    MutualfundFactorRetMap = containers.Map();
    MutualfundFactorExpMap = containers.Map();
    MutualfundRetDecompMap = containers.Map();
    
    T = 500;
    hl = 63;
    wReg = ExponentialWeight(T, hl);
    for ii = 1 : size(C,1)
        fundid = C{ii,1};
        disp([datestr(now) ': ' fundid]);
        navInfo = mfnavMap(fundid);

        if ~MutualfundMap.isKey(fundid)
            MutualfundLst{end+1,1} = fundid;
            MutualfundMap(fundid) = size(MutualfundLst,1);
        end

        dateAllNum = navInfo(:,1);
        navAll = navInfo(:,2);
        MutualfundDateMap(fundid) = datestr(dateAllNum(end,1),'yyyy-mm-dd');

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

        MutualfundNavMap(fundid) = oneNav;
        MutualfundFactorRetMap(fundid) = oneFctRet;
        MutualfundFactorExpMap(fundid) = oneFctExp;
        MutualfundRetDecompMap(fundid) = [oneFctExp.*oneFctRet(:,1:end-1) oneFctRet(:,end)];
    end

    save(MutualfundMatFile,  'MutualfundMap', 'MutualfundDateMap', 'MutualfundNavMap', 'MutualfundFactorRetMap', 'MutualfundFactorExpMap', 'MutualfundRetDecompMap');
end


%% 图

factorsLst = {'countryFRetMtrx';'sizeFRetMtrx';'momentumFRetMtrx';'residualVolFRetMtrx';'liquidityFRetMtrx';}; 
for ii=1:size(C,1)
    %ii=21;
    fundid = C{ii};
    oneFctRet = MutualfundFactorRetMap(fundid); 
    oneFctExp = MutualfundFactorExpMap(fundid); 
    oneRetDmp = MutualfundRetDecompMap(fundid); 
    oneNav = MutualfundNavMap(fundid);
    clf;
    figure(1);
    subplot(3,1,1);
    date = oneNav(:,4);
    accumret = cumprod(nansum(oneRetDmp,2)+1)-1; % 累计收益率
    plot(date, accumret);
    hold on; 
    plot(date(2:end), cumprod(oneNav(2:end,3)+1)-1);
    title(char(C{ii,2}))
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
    saveas(gcf,['.\figures\公募_' char(C{ii,2}) '.jpg'])
    pause(1)
    
end

