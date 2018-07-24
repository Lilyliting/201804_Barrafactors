%% Get Factor Portfolio Return Matrix
ngroup = input('how many groups you want: ');

daysLst = [1, 5, 10, 20];
fctExptrLst = {'sizeFExpMtrx', 'momentumFExpMtrx', 'alphaFExpMtrx', 'betaFExpMtrx', 'residualVolFExpMtrx', 'liquidityFExpMtrx'};
fctPtfRetLst = cell(0,0);
for ii = 1 : size(fctExptrLst,2)
    fctKey = strrep(fctExptrLst{ii},'ExpMtrx','PtfRetMtrx');
    fctPtfRetLst{1,end+1} = fctKey;
end

% %factorPortfolioMatFile
% if exist(factorPortfolioMatFile,'file')
%     load(factorPortfolioMatFile);
% else
%     factorPortfolioRetMap = containers.Map();
%     save(factorPortfolioMatFile,'factorPortfolioRetMap');
% end

for ii = 1 : size(fctPtfRetLst,2)
    fctKey = fctPtfRetLst{ii};
    factorPortfolioRetMap(fctKey) = zeros(0,0,0);
    
    onePtfRetMtrx = factorPortfolioRetMap(fctKey);
    onePtfRetMtrx(size(daysLst,2),size(tradingDaysLst,1),ngroup) = 0;  % (4 * T * 10)
    
    disp([datestr(now) ': ' fctKey]);
    posSta = tradingDaysMap(firsttradeday);
    posEnd = size(tradingDaysLst,1);
    factorDateMap(fctKey) = tradingDaysLst{end,1};
    for jj = posSta : posEnd %时间
        for kk = 1 : size(daysLst,2) %四个日期
            dayN = daysLst(kk);
            rt = stkCloseFMtrx(jj,:) ./ stkCloseFMtrx(jj-1,:) - 1;
            
            onePtfRetMtrx(kk,jj,:) = 0;
            for mm = 1 : dayN
                oneExp = real(eval([fctExptrLst{ii} '(' num2str(jj - mm - 1) ',:)']));
                oneExp(1,stkVolumeMtrx(jj - mm,:) == 0) = nan;  % 去除交易量为零的股票（停牌、涨跌停）
                prtl = prctile(oneExp,0:(100/ngroup):100); % group stocks by exposure
                wt = stkEVMtrx(jj-mm-1,:);
                for ll = 1 : ngroup
                    if ll == 1
                        stkGroupPos = oneExp <= prtl(ll+1);
                    else
                        stkGroupPos = oneExp > prtl(ll) & oneExp <= prtl(ll+1);
                    end
                    rtNow = rt(stkGroupPos);
                    wtNow = wt(stkGroupPos) / nansum(wt(stkGroupPos));
                    onePtfRetMtrx(kk,jj,ll) = onePtfRetMtrx(kk,jj,ll) + nansum(rtNow .* wtNow) / dayN;
                end
            end
        end
    end
    
    factorPortfolioRetMap(fctKey) = onePtfRetMtrx;
end

lgd = cell(1,ngroup);
for ii = 1:ngroup
    lgd{ii}=num2str(ii);
end
%% 单因子检验
for kk = 1:4
    disp([num2str(daysLst(kk)) ' day(s) accumulated']);
    for ii=1:size(fctPtfRetLst,2)
        onePtfRetMtrx = factorPortfolioRetMap(fctPtfRetLst{ii});
        onePtfRetMtrx = onePtfRetMtrx(:,1300:end,:);
        tmp = reshape(onePtfRetMtrx(4,:,:),size(onePtfRetMtrx,2),size(onePtfRetMtrx,3));
        ast = cumprod(tmp+1,1)-1;
        col = jet(ngroup);
        clf;
        figure(1);
        
        for jj=1:ngroup
            plot(ast(:,jj),'color',col(jj,:));hold on;
        end
        titlename = [fctPtfRetLst{ii} '(' num2str(daysLst(kk)) 'days)'];
        title(titlename);
        legend(lgd, 'Location','northwest');
        saveas(gcf,['.\figures\factor_' titlename '.jpg'])
        % legend({'1','2','3','4','5','6','7','8','9','10'},'Location','northwest');
        pause(1)
    end
end


%% factorPortfolioOptmWgtMatFile
% fctExptrLst = {'sizeFExpMtrx', 'momentumFExpMtrx', 'alphaFExpMtrx', 'betaFExpMtrx', 'residualVolFExpMtrx', 'liquidityFExpMtrx'};
fctPtfRetMtrx = [];
for ii = 1 : size(fctPtfRetLst,2)
    fctKey = fctPtfRetLst{ii};
    onePtfRetMtrx = factorPortfolioRetMap(fctKey);
    fctPtfRetMtrx = [fctPtfRetMtrx reshape(onePtfRetMtrx(3,:,:),size(onePtfRetMtrx,2),size(onePtfRetMtrx,3))];
end

% mean-variance efficient frontier
T = 150;
hl = 60;
wExpPtfRet = ExponentialWeight(T, hl);
tmpRetMtrx = fctPtfRetMtrx(end-T:end,:) .* repmat(wExpPtfRet, 1, size(fctPtfRetMtrx,2));
ExpReturn = sum(tmpRetMtrx,1);
ExpCovariance = cov(tmpRetMtrx);
NumPorts = 100;
[PortRisk, PortReturn, PortWts] = frontcon(ExpReturn,ExpCovariance, NumPorts);
tmp = fctPtfRetMtrx(end-T:end,:) * PortWts';
tmp = [sum(tmp,1);sqrt(T)*std(tmp,1)];
figure(2);
plot(tmp(2,:),tmp(1,:));
xlabel('risk');ylabel('return');

%
% if exist(factorPortfolioOptmWgtMatFile,'file')
%     load(factorPortfolioOptmWgtMatFile);
% else
%
% end
