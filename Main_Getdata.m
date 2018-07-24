%% Define wind API
% w = windmatlab;

%% Connect to mysql
conn = database('wind','hulitingali','Hu.lt@2018','com.mysql.jdbc.Driver','jdbc:mysql://rm-2zey1z6px42nits51.mysql.rds.aliyuncs.com/wind');

%% Define work path
% userpath('C:\Users\Administrator\Desktop\FactorAnalysis');
% savepath

%% Parameters
tradingDaysMatFile = '.\mat\tradingDays.mat'; %放所有交易日以MATLAB格式
stkLstXlsFile = '.\xls\stockList.xlsx'; %放所有股票代码
stkLstMatFile = '.\mat\stockList.mat'; %放所有股票代码以MATLAB格式
stkQuotMatFile = '.\mat\stkQuotations.mat';
idxLstMatFile = '.\mat\indexList.mat';
idxQuotMatFile = '.\mat\idxQuotations.mat';
factorsXlsFile = '.\xls\factors.xlsx';
factorsMatFile = '.\mat\factors.mat';
factorExposureMatFile = '.\mat\factorExposure.mat';
factorReturnMatFile = '.\mat\factorReturn.mat';
factorPortfolioMatFile = '.\mat\factorPortfolio.mat';
factorPortfolioOptmWgtMatFile = '.\mat\factorPortfolioOptmWgt.mat';
portfolioMatFile = '.\mat\portfolio.mat';
portfolioXlsPath = '.\xls\';
bondXlsPath ='.\bond\';

staDate = '2007-01-01';
stadate = datestr(staDate,'yyyymmdd');

%ed = '2017-11-01';
ed = now;
enddate = datestr(ed,'yyyymmdd');
endDate = datestr(ed, 'yyyy-mm-dd');
endDate0 = endDate;

%% Get Trading days
%确定交易日开始的第一天
if exist(tradingDaysMatFile,'file')
    load(tradingDaysMatFile);
else
    %清空tradingDaysLst的mat文件，并重新确定开始日期
    tradingDaysLst = cell(0,0);
    tradingDaysMap = containers.Map();
    futExpDayMtrx = zeros(0,0);
    stadate = staDate;
    
    %提取万得交易日数据并定义交易日序列
    query = ['select TRADE_DAYS from asharecalendar where S_INFO_EXCHMARKET = ''SSE'' and TRADE_DAYS > ''' , ...
        stadate, ''' and TRADE_DAYS < ''', enddate, ''' order by TRADE_DAYS;'];
    curs = exec(conn, query);
    fetchcurs = fetch(curs);
    Tradedays = fetchcurs.data;
    w_tdays_times = datenum(Tradedays,'yyyymmdd');


    if iscell(Tradedays)
        posEnd = size(tradingDaysLst,1);
        for ii = 1 : size(Tradedays,1)
            tradingDaysLst{end+1,1} = datestr(w_tdays_times(ii),'yyyy-mm-dd');
            tradingDaysMap(tradingDaysLst{ii,1}) = posEnd + ii;
        end

        %计算期货合约到期天数并存成tradingdays的mat格式
        futExpDayMtrx(size(tradingDaysLst,1),1) = nan;
        this1st = [datestr(staDate,'yyyy-mm-') '01'];
        this3rdFriday = datenum(this1st) + 20 - rem(weekday(this1st),7);
        for ii = 1 : size(w_tdays_times,1)
            expDays = this3rdFriday - w_tdays_times(ii);
            if expDays < 0
                this1st = [datestr(w_tdays_times(ii)+25,'yyyy-mm-') '01'];
                this3rdFriday = datenum(this1st) + 20 - rem(weekday(this1st),7);
                expDays = this3rdFriday - w_tdays_times(ii);
            end
            futExpDayMtrx(posEnd + ii,1) = expDays;
        end
        save(tradingDaysMatFile, 'tradingDaysLst', 'tradingDaysMap', 'futExpDayMtrx');
    end
end





%% Get All Future Code
if exist(futLstMatFile,'file')
    load(futLstMatFile);
else
    futCodeLst = {'IF00.CFE';'IF01.CFE';'IF02.CFE';'IF03.CFE';'IH00.CFE';'IH01.CFE';'IH02.CFE';'IH03.CFE';'IC00.CFE';'IC01.CFE';'IC02.CFE';'IC03.CFE';};
    futCodeMap = containers.Map();
    futCodeDateMap = containers.Map();
    for ii = 1 : size(futCodeLst,1)
        futCodeMap(futCodeLst{ii}) = ii;
        futCodeDateMap(futCodeLst{ii}) = staDate;
    end
    save(futLstMatFile, 'futCodeLst', 'futCodeMap', 'futCodeDateMap');
end

% Get Future Quotations: prices, volumes
if exist(futQuotMatFile,'file')
    load(futQuotMatFile);
else
    futSettleMtrx = zeros(0,0);
    futOiMtrx = zeros(0,0);
    futVolumeMtrx = zeros(0,0);
end
% nowDate = datestr(now, 'yyyy-mm-dd');
futSettleMtrx(size(tradingDaysLst,1),size(futCodeLst,1)) = nan;
futOiMtrx(size(tradingDaysLst,1),size(futCodeLst,1)) = nan;
futVolumeMtrx(size(tradingDaysLst,1),size(futCodeLst,1)) = nan;

%% 提取期货数据 未完成
% 
% for ii = 1 : size(futCodeLst,1)
%     disp([datestr(now) ': ' futCodeLst{ii,1} '  ' futCodeDateMap(futCodeLst{ii,1}) '  ' endDate]);
%     % CAUTION: this data has now data, but the quotation is at yesterday
%     
%     query3 = ['select TRADE_DT,S_DQ_TURN, S_VAL_MV,TOT_SHR_TODAY,FREE_SHARES_TODAY from ashareeodderivativeindicator where S_INFO_WINDCODE = ''' ...
%     stkCode ''' and TRADE_DT > ''' staDate ''' and TRADE_DT < ''' endDate ''' and S_DQ_TURN is not null order by TRADE_DT;'];
%     curs = exec(conn, query3);
%     fetchcurs = fetch(curs);
%     Futuredata = fetchcurs.data;
%     
%     [w_wsd_data,w_wsd_codes,w_wsd_fields,w_wsd_times,w_wsd_errorid,w_wsd_reqid]=w.wsd(futCodeLst{ii,1},'His_settle,His_oi,His_volume',futCodeDateMap(futCodeLst{ii,1}), endDate,'unit=1');
%     if w_wsd_errorid == 0
%         timeStrLst = datestr(w_wsd_times, 'yyyy-mm-dd');
%         dateSta = timeStrLst(1,:);
%         dateEnd = timeStrLst(end,:);
%         posSta = tradingDaysMap(dateSta);
%         posEnd = tradingDaysMap(dateEnd);
%         posStk = futCodeMap(futCodeLst{ii,1});
%         
%         futSettleMtrx(posSta:posEnd,posStk) = w_wsd_data(:,1);
%         futOiMtrx(posSta:posEnd,posStk) = w_wsd_data(:,2);
%         futVolumeMtrx(posSta:posEnd,posStk) = w_wsd_data(:,3);
%         
%         futCodeDateMap(futCodeLst{ii,1}) = dateEnd;
%     end
% end

save(futLstMatFile, 'futCodeLst', 'futCodeMap', 'futCodeDateMap');
save(futQuotMatFile, 'futSettleMtrx', 'futOiMtrx', 'futVolumeMtrx');

%% Get All Indices Code

if exist(idxLstMatFile,'file')
    load(idxLstMatFile);
else
    idxCodeLst = {'000300.SH';'000016.SH';'000905.SH';'SHIBOR1Y.IR';};
    idxCodeMap = containers.Map();
    idxCodeDateMap = containers.Map();
    for ii = 1 : size(idxCodeLst,1)
        idxCodeMap(idxCodeLst{ii}) = ii;
        idxCodeDateMap(idxCodeLst{ii}) = staDate;
    end
    save(idxLstMatFile, 'idxCodeLst', 'idxCodeMap', 'idxCodeDateMap');
end

% Get Indices Quotations: prices, volumes
if exist(idxQuotMatFile,'file')
    load(idxQuotMatFile);
else
    idxCloseMtrx = zeros(0,0);
    idxAmtMtrx = zeros(0,0);
    % nowDate = datestr(now, 'yyyy-mm-dd');
    idxCloseMtrx(size(tradingDaysLst,1),size(idxCodeLst,1)) = nan;
    idxAmtMtrx(size(tradingDaysLst,1),size(idxCodeLst,1)) = nan;
    
    for ii = 1 : (size(idxCodeLst,1)-1)
        disp([datestr(now) ': ' idxCodeLst{ii,1} '  ' idxCodeDateMap(idxCodeLst{ii,1}) '  ' endDate]);
        % CAUTION: this data has now data, but the quotation is at yesterday

        query4 = ['select TRADE_DT, S_DQ_CLOSE, S_DQ_AMOUNT from aindexeodprices where S_INFO_WINDCODE = ''' ...
        idxCodeLst{ii,1} ''' and TRADE_DT > ''' stadate ''' and TRADE_DT < ''' enddate ''' order by TRADE_DT;'];
        curs = exec(conn, query4);
        fetchcurs = fetch(curs);
        indexdata = fetchcurs.data;
        
        w_wsd_times = datenum(indexdata(:,1), 'yyyymmdd');
        timeStrLst = datestr(w_wsd_times, 'yyyy-mm-dd');
        
        posarray = zeros(size(timeStrLst,1),1);
        for jj=1:size(timeStrLst,1)
            posarray(jj) = tradingDaysMap(timeStrLst(jj,:));
        end
        
%         dateSta = timeStrLst(1,:);
%         dateEnd = timeStrLst(end,:);
%         tradingDaysMap(timeStrLst)
%         posSta = tradingDaysMap(dateSta);
%         posEnd = tradingDaysMap(dateEnd);
        posStk = idxCodeMap(idxCodeLst{ii,1});

        temp = cell2mat(indexdata(:,2));
        idxCloseMtrx(posarray,posStk) = temp;
        temp = cell2mat(indexdata(:,3));
        idxAmtMtrx(posarray,posStk) = temp;
        idxCodeDateMap(idxCodeLst{ii,1}) = dateEnd;
        
    %     [w_wsd_data,w_wsd_codes,w_wsd_fields,w_wsd_times,w_wsd_errorid,w_wsd_reqid]=w.wsd(idxCodeLst{ii,1},'close,amt',idxCodeDateMap(idxCodeLst{ii,1}), endDate,'unit=1');
    %     if w_wsd_errorid == 0
    %         timeStrLst = datestr(w_wsd_times, 'yyyy-mm-dd');
    %         dateSta = timeStrLst(1,:);
    %         dateEnd = timeStrLst(end,:);
    %         posSta = tradingDaysMap(dateSta);
    %         posEnd = tradingDaysMap(dateEnd);
    %         posStk = idxCodeMap(idxCodeLst{ii,1});
    %         
    %         idxCloseMtrx(posSta:posEnd,posStk) = w_wsd_data(:,1);
    %         idxAmtMtrx(posSta:posEnd,posStk) = w_wsd_data(:,2);
    %         
    %         idxCodeDateMap(idxCodeLst{ii,1}) = dateEnd;
    %     end
    end
    
    disp([datestr(now) ': ' idxCodeLst{4,1} '  ' idxCodeDateMap(idxCodeLst{4,1}) '  ' endDate]);
        
    tmp = dir([bondXlsPath '*年中债国债收益率曲线标准期限信息.xlsx']); %提取国债收益曲线的文件名称
    tot = size(tmp,1);
    riskfree = [];
    for ii = 1 : tot
        tmpname = tmp(ii).name;
        tmppath = [bondXlsPath tmpname];
        [~,~,C] = xlsread(tmppath);
        idx = find(cell2mat(C(2:end,3))==1);
        if ii==12
            dateseq = datenum(C(idx+1,1),'yyyy-mm-dd');
        else
            dateseq = datenum(C(idx+1,1),'yyyy/mm/dd');
        end
        
        bondquote = cell2mat(C(idx+1,4));
        tmpbond = [dateseq bondquote];
        riskfree = [riskfree; tmpbond];
    end
    
    kk = 1;
    for ii = 1:size(riskfree,1)
        try
            t = find(w_wsd_times==riskfree(ii,1));
            idxCloseMtrx(kk,4) = riskfree(t,2);
            kk = kk+1;
        catch
        end
    end

    idxAmtMtrx(:,4)=NaN;
    save(idxQuotMatFile, 'idxCloseMtrx', 'idxAmtMtrx');
end