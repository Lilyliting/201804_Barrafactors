%%
% Connect to mysql.wind
conn = database('wind','hulitingali','Hu.lt@2018','com.mysql.jdbc.Driver','jdbc:mysql://rm-2zey1z6px42nits51.mysql.rds.aliyuncs.com/wind');


% parameters
stkLstXlsFile = '.\xls\stockList.xlsx'; %放所有股票代码
stkLstMatFile = '.\mat\stockList.mat'; %放所有股票代码以MATLAB格式
stkCodeDateMatFile = '.\mat\stockDay.mat';
stkdataFile = '.\mat\data\';

stadate = '20070101';
setend = now;
enddate = datestr(setend,'yyyymmdd');

% ashareeodprices
info{1} = {'close','adjclose','volume'};
qry{1} = ['select TRADE_DT,S_INFO_WINDCODE,S_DQ_CLOSE,S_DQ_ADJCLOSE,'...
    'S_DQ_VOLUME from ashareeodprices '];


% ashareeodderivativeindicator
info{2} = {'turn','EV','totalshare','freefloatshare','PE','PE_TTM','PB',...
    'PCF_OCF','PCF_OCFTTM','PCF_NCF','PCF_NCFTTM','PS','PS_TTM',...
    'divpershare','OCF_TTM','OCF_LYR'};
qry{2} = ['select TRADE_DT,S_INFO_WINDCODE,S_DQ_TURN, S_VAL_MV,TOT_SHR_TODAY,' ...
    'FREE_SHARES_TODAY,S_VAL_PE,S_VAL_PE_TTM,S_VAL_PB_NEW,S_VAL_PCF_OCF,'...
    'S_VAL_PCF_OCFTTM,S_VAL_PCF_NCF,S_VAL_PCF_NCFTTM,S_VAL_PS,S_VAL_PS_TTM,'...
    'S_PRICE_DIV_DPS,NET_CASH_FLOWS_OPER_ACT_TTM,NET_CASH_FLOWS_OPER_ACT_LYR',...
    ' from ashareeodderivativeindicator'];

% asharefinancialindicator
info{3} = {'grossprofitmargin','netprofitmargin','ROE','ROA','debt',...
    'cashratio','quickratio','currentratio'};
qry{3} = ['select ANN_DT,S_INFO_WINDCODE,S_FA_GROSSPROFITMARGIN,'...
    'S_FA_NETPROFITMARGIN,S_FA_ROE,S_FA_ROA2,S_FA_DEBTTOASSETS,'...
    'S_FA_CASHRATIO,S_FA_QUICK,S_FA_CURRENT from asharefinancialindicator'];

% ashareincome
info{4} = {'revenue','netprofit','EBITDA','EPS','diluted_EPS'};
qry{4} = ['select ANN_DT,S_INFO_WINDCODE,OPER_REV,NET_PROFIT_EXCL_MIN_INT_INC,'...
    'EBITDA,S_FA_EPS_BASIC,S_FA_EPS_DILUTED from ashareincome'];

% asharetechindicators
info{5} = {'MACD_DIFF','MACD_DEA'};
qry{5} = ['select TRADE_DT,S_INFO_WINDCODE,MACD_DIFF,MACD_DEA '...
    'from ashareintensitytrendadj'];

%% Get All Stocks Code: stkLstXlsFile must be updated manually periodically
% 检验没存放股票代码的文件是否存在加载文件，或者创造新的mat文件并清空
[~,~,C] = xlsread(stkLstXlsFile);

stkCodeLst = cell(0,0);
stkCodeMap = containers.Map();
stkCodeDateMap = containers.Map();
validstk = containers.Map();
%读取每个股票代码并存成mat格式文件
jj = 1;
for ii = 1 : size(C,1)
    stkId = C{ii,1};
    if isnan(stkId)
        break;
    elseif ~isempty(strfind(C{ii,2},'*'))
        continue;
    else %~stkCodeMap.isKey(stkId)
        stkCodeLst{end+1,1} = stkId;
        stkCodeMap(stkId) = jj;
        jj = jj+1;
        %         stkCodeDateMap(stkId) = stadate;
    end
end
save(stkLstMatFile, 'stkCodeLst', 'stkCodeMap');
stocknum = size(stkCodeMap,1);

%% 获取股票信息并展示进度
% fetchstkdatabase(stkCodeLst,enddate,stkinfo1,qry1,10,1)
% fetchstkdatabase(stkCodeLst,enddate,stkinfo2,qry2,10,2);
% fetchstkdatabase(stkCodeLst,enddate,stkinfo3,qry3,10,3);
% fetchstkdatabase(stkCodeLst,enddate,stkinfo4,qry4,10,4);
% fetchstkdatabase(stkCodeLst,enddate,stkinfo5,qry5,10,5);

n = 20;
opt = 1;
stkinfo = info{opt};
query = qry{opt};




errCnt = 0;
totNum = size(stkCodeLst,1);
if exist(stkCodeDateMatFile, 'file')
    load(stkCodeDateMatFile, 'stkCodeDateMap');
end

yesorno = input('download stocks?y/n:');
if yesorno == 'y'
    for ii = 1: floor(totNum/n)
        tic;
        %     try
        %         begindate = stkCodeDateMap([stkCodeLst{(ii-1)*n+1},num2str(opt)]);
        %         if strcmp(enddate,begindate)
        %             continue;
        %         end
        %     catch
        %         begindate = stadate;
        %     end
        begindate = stadate;
        % CAUTION: this data has now data, but the quotation is at yesterday
        stkcodes = strjoin(stkCodeLst(((ii-1)*n+1):ii*n),''',''');
        
        
        if strfind(query,'ANN_DT')
            dstr = 'ANN_DT';
        else
            dstr = 'TRADE_DT';
        end
        wholequery = [query ' where S_INFO_WINDCODE in (''' stkcodes ...
            ''') and ' dstr ' >= ''' begindate ''' and ' dstr ' < ''' enddate ...
            '''order by ' dstr ';'];
        curs = exec(conn, wholequery);
        fetchcurs = fetch(curs);
        mysqldata = fetchcurs.data;
        
        if strcmp(mysqldata{1},'No Data')
            continue;
        end
        
        for jj = 1:n
            try
                idx = (ii-1)*n+jj;
                matFileName = [stkdataFile stkCodeLst{idx,1}(1:6) '.mat'];
                
                stkCode = stkCodeLst{idx,1};
                disp([datestr(now) ': ' stkCode '  ' num2str(idx) '/' num2str(totNum) ' ' ...
                    begindate  ' - ' enddate '   errCnt = ' num2str(errCnt)])
                
                stkidx = sum(strcmp(mysqldata,stkCode),2);
                stkdata = mysqldata(stkidx==1,:);
                
                tsA = datenum(stkdata(:,1), 'yyyymmdd');
                
                stkAdd = zeros(size(tsA,1),size(stkinfo,2)+1);
                stkAdd(:,1) = tsA;
                for kk = 1:size(stkinfo,2)
                    stkAdd(:,kk+1)=cell2mat(stkdata(:,kk+2));
                end
                
                
                if exist(matFileName, 'file')
                    %                 if opt == 1
                    %                     stockdata = cell();
                    %                 else
                    %                     load(matFileName, 'stockdata');
                    %                 end
                    load(matFileName, 'stockdata');
                    stockdata{opt} = stkAdd;
                else
                    stockdata = cell(0);
                    stockdata{opt} = stkAdd;
                end
                
                
                stkCodeDateMap([stkCodeLst{idx,1},num2str(opt)]) = enddate;
                save(stkCodeDateMatFile, 'stkCodeDateMap');
                save(matFileName,'stockdata')
                validstk(ii) = stkCode;
            catch
                errCnt = errCnt + 1;
            end
        end
        toc;
    end
    
    for ii= (totNum-mod(totNum,n)+1) : totNum
        tic;
        
        % CAUTION: this data has now data, but the quotation is at yesterday
        stkcode = stkCodeLst{ii};
        %     try
        %         begindate = stkCodeDateMap([stkcode,num2str(opt)]);
        %         if strcmp(enddate,begindate)
        %             continue;
        %         end
        %     catch
        %         begindate = stadate;
        %     end
        begindate = stadate;
        if strfind(query,'ANN_DT')
            dstr = 'ANN_DT';
        else
            dstr = 'TRADE_DT';
        end
        wholequery = [query ' where S_INFO_WINDCODE in (''' stkcode ...
            ''') and ' dstr ' >= ''' begindate ''' and ' dstr ' < ''' enddate ...
            '''order by ' dstr ';'];
        curs = exec(conn, wholequery);
        fetchcurs = fetch(curs);
        stkdata = fetchcurs.data;
        
        if strcmp(stkdata{1},'No Data')
            continue;
        end
        
        try
            matFileName = [stkdataFile stkCodeLst{ii,1}(1:6) '.mat'];
            disp([datestr(now) ': ' stkcode '  ' num2str(ii) '/' num2str(totNum) ' ' ...
                begindate  ' - ' enddate '   errCnt = ' num2str(errCnt)])
            
            tsA = datenum(stkdata(:,1), 'yyyymmdd');
            
            stkAdd = zeros(size(tsA,1),size(stkinfo,2)+1);
            stkAdd(:,1) = tsA;
            for kk = 1:size(stkinfo,2)
                stkAdd(:,kk+1)=cell2mat(stkdata(:,kk+2));
            end
            
            if exist(matFileName, 'file')
                %                 if opt == 1
                %                     stockdata = cell();
                %                 else
                %                     load(matFileName, 'stockdata');
                %                 end
                load(matFileName, 'stockdata');
                stockdata{opt} = stkAdd;
            else
                stockdata{opt} = stkAdd;
            end
            
            stkCodeDateMap([stkCodeLst{idx,1},num2str(opt)]) = enddate;
            save(stkCodeDateMatFile, 'stkCodeDateMap');
            save(matFileName,'stockdata')
            validstk(ii) = stkcode;
        catch
            errCnt = errCnt + 1;
        end
    end
    toc;
    disp(['---------------------Database ' num2str(opt) ' completed'])
end


%% 提取股票信息
% Get Stocks Quotations: prices, volumes
if exist(stkQuotMatFile,'file')
    load(stkQuotMatFile);
else
    fs = dir([stkdataFile,'*.mat']);
    totNum2 = size(fs,1);
    load([stkdataFile fs(3).name],'stockdata')
    Timeseries = stockdata{1}(:,1);
    tt = size(Timeseries,1);
    
    %建立2*2的0矩阵
    stkCloseMtrx = zeros(tt,totNum2); %1
    stkAdjcloseMtrx = zeros(tt,totNum2);
    stkVolumeMtrx = zeros(tt,totNum2);
    stkTurnOverMtrx = zeros(tt,totNum2);
    stkEVMtrx = zeros(tt,totNum2); %5
    stkTotalShareMtrx = zeros(tt,totNum2);
    stkFreeFloatMtrx = zeros(tt,totNum2);
    stkPEMtrx = zeros(tt,totNum2);
    stkPE_TTMMtrx = zeros(tt,totNum2);
    stkPBMtrx = zeros(tt,totNum2); %10
    stkPCF_OCFMtrx = zeros(tt,totNum2);
    stkPCF_OCFTTMMtrx = zeros(tt,totNum2);
    stkPCF_NCFMtrx = zeros(tt,totNum2);
    stkPCF_NCFTTMMtrx = zeros(tt,totNum2);
    stkPSMtrx = zeros(tt,totNum2); %15
    stkPS_TTMMtrx = zeros(tt,totNum2);
    stkDivPerShareMtrx = zeros(tt,totNum2);
    stkOCF_TTMMtrx = zeros(tt,totNum2);
    stkOCF_LYRMtrx = zeros(tt,totNum2);
    stkGrossProfitMarginMtrx = zeros(tt,totNum2); %20
    stkNetProfitMarginMtrx = zeros(tt,totNum2);
    stkROEMtrx = zeros(tt,totNum2);
    stkROAMtrx = zeros(tt,totNum2);
    stkDebtMtrx = zeros(tt,totNum2);
    stkCashRatioMtrx = zeros(tt,totNum2); %25
    stkQuickRatioMtrx = zeros(tt,totNum2);
    stkCurrentRatioMtrx = zeros(tt,totNum2);
    stkRevenueMtrx = zeros(tt,totNum2);
    stkNetProfitMtrx = zeros(tt,totNum2);
    stkEBITDAMtrx = zeros(tt,totNum2); %30
    stkEPSMtrx = zeros(tt,totNum2);
    stkDiluted_EPSMtrx = zeros(tt,totNum2);
    stkMACD_DIFFMtrx = zeros(tt,totNum2);
    stkMACD_DEAMtrx = zeros(tt,totNum2);
    
    for ii = 1 : totNum2
        disp([datestr(now) ': ' slist(ii).name '  ' num2str(ii) '/' num2str(totNum2)]);
        % CAUTION: this data has now data, but the quotation is at yesterday
        matFileName = ['C:\Users\Administrator\Desktop\FactorAnalysis\mat\data\' slist(ii).name];
        load(matFileName,'stockdata');
        
        stkCode = [slist(ii).name(:,1:6) '.SZ'];
        
        % table 1
        [yn,locb] = ismember(Timeseries,stockdata{1}(:,1));
        valididx = find(yn==1);
        locb(locb==0) = [];
        stktmp = stockdata{1}(locb,2); %1
        stkCloseMtrx(valididx,ii) = stktmp;
        stktmp = stockdata{1}(locb,3);
        stkAdjcloseMtrx(valididx,ii) = stktmp;
        stktmp = stockdata{1}(locb,4); %3
        stkVolumeMtrx(valididx,ii) = stktmp;
        
        % table 2
        [yn,locb] = ismember(Timeseries,stockdata{2}(:,1));
        valididx = find(yn==1);
        locb(locb==0) = [];
        stktmp = stockdata{2}(locb,2); %4
        stkTurnOverMtrx(valididx,ii) = stktmp;
        stktmp = stockdata{2}(locb,3);
        stkEVMtrx(valididx,ii) = stktmp;
        stktmp = stockdata{2}(locb,4); %6
        stkTotalShareMtrx(valididx,ii) = stktmp;
        stktmp = stockdata{2}(locb,5);
        stkFreeFloatMtrx(valididx,ii) = stktmp;
        stktmp = stockdata{2}(locb,6); %8
        stkPEMtrx(valididx,ii) = stktmp;
        stktmp = stockdata{2}(locb,7);
        stkPE_TTMMtrx(valididx,ii) = stktmp;
        stktmp = stockdata{2}(locb,8); %10
        stkPBMtrx(valididx,ii) = stktmp;
        stktmp = stockdata{2}(locb,9);
        stkPCF_OCFMtrx(valididx,ii) = stktmp;
        stktmp = stockdata{2}(locb,10); %12
        stkPCF_OCFTTMMtrx(valididx,ii) = stktmp;
        stktmp = stockdata{2}(locb,11);
        stkPCF_NCFMtrx(valididx,ii) = stktmp;
        stktmp = stockdata{2}(locb,12); %14
        stkPCF_NCFTTMMtrx(valididx,ii) = stktmp;
        stktmp = stockdata{2}(locb,13);
        stkPSMtrx(valididx,ii) = stktmp;
        stktmp = stockdata{2}(locb,14); %16
        stkPS_TTMMtrx(valididx,ii) = stktmp;
        stktmp = stockdata{2}(locb,15);
        stkDivPerShareMtrx(valididx,ii) = stktmp;
        stktmp = stockdata{2}(locb,16); %18
        stkOCF_TTMMtrx(valididx,ii) = stktmp;
        stktmp = stockdata{2}(locb,17);
        stkOCF_LYRMtrx(valididx,ii) = stktmp;
        
        
        % table 3
        [yn,locb] = ismember(Timeseries,stockdata{3}(:,1));
        valididx = find(yn==1);
        locb(locb==0) = [];
        stktmp = stockdata{3}(locb,2); %20
        stkGrossProfitMarginMtrx(valididx,ii) = stktmp;
        stktmp = stockdata{3}(locb,3);
        stkNetProfitMarginMtrx(valididx,ii) = stktmp;
        stktmp = stockdata{3}(locb,4); %22
        stkROEMtrx(valididx,ii) = stktmp;
        stktmp = stockdata{3}(locb,5);
        stkROAMtrx(valididx,ii) = stktmp;
        stktmp = stockdata{3}(locb,6); %24
        stkDebtMtrx(valididx,ii) = stktmp;
        stktmp = stockdata{3}(locb,7);
        stkCashRatioMtrx(valididx,ii) = stktmp;
        stktmp = stockdata{3}(locb,8); %26
        stkQuickRatioMtrx(valididx,ii) = stktmp;
        stktmp = stockdata{3}(locb,9);
        stkCurrentRatioMtrx(valididx,ii) = stktmp;
        
        % table 4
        [yn,locb] = ismember(Timeseries,stockdata{4}(:,1));
        valididx = find(yn==1);
        locb(locb==0) = [];
        stktmp = stockdata{4}(locb,2); %28
        stkRevenueMtrx(valididx,ii) = stktmp;
        stktmp = stockdata{4}(locb,3);
        stkNetProfitMtrx(valididx,ii) = stktmp;
        stktmp = stockdata{4}(locb,4); %30
        stkEBITDAMtrx(valididx,ii) = stktmp;
        stktmp = stockdata{4}(locb,5);
        stkEPSMtrx(valididx,ii) = stktmp;
        stktmp = stockdata{4}(locb,6); %32
        stkDiluted_EPSMtrx(valididx,ii) = stktmp;
        
        % table 5
        [yn,locb] = ismember(Timeseries,stockdata{5}(:,1));
        valididx = find(yn==1);
        locb(locb==0) = [];
        stktmp = stockdata{5}(locb,2); %33
        stkMACD_DIFFMtrx(valididx,ii) = stktmp;
        stktmp = stockdata{5}(locb,3);
        stkMACD_DEAMtrx(valididx,ii) = stktmp;
        
    end
    
    %stkEVMtrx(isnan(stkEVMtrx))=0;
    
    stkGrossProfitMarginMtrx = nandatafill(stkGrossProfitMarginMtrx); %20
    stkNetProfitMarginMtrx = nandatafill(stkNetProfitMarginMtrx);
    stkROEMtrx = nandatafill(stkROEMtrx);
    stkROAMtrx = nandatafill(stkROAMtrx);
    stkDebtMtrx = nandatafill(stkDebtMtrx);
    stkCashRatioMtrx = nandatafill(stkCashRatioMtrx); %25
    stkQuickRatioMtrx = nandatafill(stkQuickRatioMtrx);
    stkCurrentRatioMtrx = nandatafill(stkCurrentRatioMtrx);
    stkRevenueMtrx = nandatafill(stkRevenueMtrx);
    stkNetProfitMtrx = nandatafill(stkNetProfitMtrx);
    stkEBITDAMtrx = nandatafill(stkEBITDAMtrx); %30
    stkEPSMtrx = nandatafill(stkEPSMtrx);
    stkDiluted_EPSMtrx = nandatafill(stkDiluted_EPSMtrx);
    stkMACD_DIFFMtrx = nandatafill(stkMACD_DIFFMtrx);
    stkMACD_DEAMtrx = nandatafill(stkMACD_DEAMtrx);
    
    
    posNan = stkAdjcloseMtrx==0;
    stkAdjcloseMtrx(posNan) = NaN;
    stkVolumeMtrx(posNan) = NaN;
    
    stkTurnOverMtrx(posNan) = NaN;
    stkEVMtrx(posNan) = NaN; %5
    stkTotalShareMtrx(posNan) = NaN;
    stkFreeFloatMtrx(posNan) = NaN;
    
    stkPEMtrx(posNan) = NaN;
    stkPE_TTMMtrx(posNan) = NaN;
    stkPBMtrx(posNan) = NaN; %10
    stkPCF_OCFMtrx(posNan) = NaN;
    
    stkPCF_OCFTTMMtrx(posNan) = NaN;
    stkPCF_NCFMtrx(posNan) = NaN;
    stkPCF_NCFTTMMtrx(posNan) = NaN;
    
    stkPSMtrx(posNan) = NaN; %15
    stkPS_TTMMtrx(posNan) = NaN;
    stkDivPerShareMtrx(posNan) = NaN;
    stkOCF_TTMMtrx(posNan) = NaN;
    
    stkOCF_LYRMtrx(posNan) = NaN;
    stkGrossProfitMarginMtrx(posNan) = NaN; %20
    stkNetProfitMarginMtrx(posNan) = NaN;
    stkROEMtrx(posNan) = NaN;
    stkROAMtrx(posNan) = NaN;
    stkDebtMtrx(posNan) = NaN;
    stkCashRatioMtrx(posNan) = NaN; %25
    stkQuickRatioMtrx(posNan) = NaN;
    stkCurrentRatioMtrx(posNan) = NaN;
    stkRevenueMtrx(posNan) = NaN;
    stkNetProfitMtrx(posNan) = NaN;
    stkEBITDAMtrx(posNan) = NaN; %30
    stkEPSMtrx(posNan) = NaN;
    stkDiluted_EPSMtrx(posNan) = NaN;
    stkMACD_DIFFMtrx(posNan) = NaN;
    stkMACD_DEAMtrx(posNan) = NaN;
    
    save(stkQuotMatFile, 'stkCloseMtrx', 'stkAdjcloseMtrx', 'stkVolumeMtrx',...
        'stkTurnOverMtrx', 'stkEVMtrx', 'stkTotalShareMtrx', 'stkFreeFloatMtrx',...
        'stkPEMtrx','stkPE_TTMMtrx','stkPBMtrx','stkPCF_OCFMtrx',...
        'stkPCF_OCFTTMMtrx','stkPCF_NCFMtrx','stkPCF_NCFTTMMtrx',...
        'stkPSMtrx','stkPS_TTMMtrx','stkDivPerShareMtrx','stkOCF_TTMMtrx',...
        'stkOCF_LYRMtrx','stkGrossProfitMarginMtrx','stkNetProfitMarginMtrx','stkROEMtrx',...
        'stkROAMtrx','stkDebtMtrx','stkCashRatioMtrx','stkQuickRatioMtrx',...
        'stkCurrentRatioMtrx','stkRevenueMtrx','stkNetProfitMtrx','stkEBITDAMtrx',...
        'stkEPSMtrx','stkDiluted_EPSMtrx','stkMACD_DIFFMtrx','stkMACD_DEAMtrx','Timeseries');
end