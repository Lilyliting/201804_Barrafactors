%% 获取公募基金数据
MutualfundLstMatFile = '.\mat\mutualfundList.mat'; %放所有基金代码 matlab格式
mfnavMatFile = '.\mat\mfnav.mat';
regstadate = '2010-01-01';
regenddate = '2017-11-01';

if exist(MutualfundLstMatFile,'file')
    load(MutualfundLstMatFile);
else
    query = 'select F_INFO_WINDCODE,S_INFO_SECTOR from chinamutualfundsector where S_INFO_SECTOR in (''2001010101000000'',''2001010201000000'',''2001010204000000'') and F_INFO_WINDCODE not like ''F%'' and CUR_SIGN=1 limit 0,100;';
    curs = exec(conn, query);
    fetchcurs = fetch(curs);
    Mutualfundname = fetchcurs.data;
    MutualfundCodeLst = Mutualfundname(:,:);
    MutualfundCodeMap = containers.Map();
    MutualfundCodeDateMap = containers.Map();
    jj=1;
    for ii = 1:size(Mutualfundname,1)
        fundid = Mutualfundname{ii,1};
        if isnan(fundid)
            break;
        elseif ~MutualfundCodeMap.isKey(fundid)
            MutualfundCodeMap(fundid) = jj;
            MutualfundCodeDateMap(fundid) = regstadate;   
            jj = jj+1;
        end
    end
    save(MutualfundLstMatFile, 'MutualfundCodeLst', 'MutualfundCodeMap', 'MutualfundCodeDateMap');
end

%% Get mutual fund nav
if exist(mfnavMatFile, 'file')
    load(mfnavMatFile)
else
    mfnavMap = containers.Map();
    for ii=1:size(MutualfundCodeLst,1)
        fundid = char(MutualfundCodeLst(ii,1));
        query = ['select ANN_DATE, F_NAV_ADJUSTED from chinamutualfundnav where F_INFO_WINDCODE = ''', fundid ,''' and ANN_DATE > ''20100101'' and ANN_DATE<''20171101'' order by ANN_DATE;'];
        curs = exec(conn, query);
        fetchcurs = fetch(curs);
        Mutualfunddata = fetchcurs.data;

        dateAllNum = datenum(Mutualfunddata(:,1), 'yyyymmdd');
        navAll = cell2mat(Mutualfunddata(:,2));

        mfnavMap(fundid) = [dateAllNum, navAll];
    end

    save(mfnavMatFile, 'mfnavMap')
end



