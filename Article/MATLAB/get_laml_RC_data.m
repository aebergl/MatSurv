cgdsURL = 'http://www.cbioportal.org/public-portal/';
try
    c = getclinicaldata(cgdsURL,'laml_tcga_pub_all');
catch
    disp('Get MSKCC CGDS Cancer Genomics Toolbox from http://www.mathworks.com/matlabcentral/fileexchange/31297-mskcc-cgds-cancer-genomics-toolbox')
    return
end

laml_RC_TimeVar = c.data(:,strcmp('OS_MONTHS',c.clinVariable));
laml_RC_EventVar = c.data(:,strcmp('OS_STATUS',c.clinVariable));
laml_RC_GroupVar =c.data(:,strcmp('RISK_CYTO',c.clinVariable));
EmptyIndx = cellfun('isempty',laml_RC_TimeVar);
if any(EmptyIndx)
    laml_RC_TimeVar(EmptyIndx)={'NaN'};
end
laml_RC_TimeVar = sscanf(sprintf('%s,',laml_RC_TimeVar{:}),'%f,');
clear c EmptyIndx cgdsURL
[p,fh,stats]=MatSurv(laml_RC_TimeVar,laml_RC_EventVar,laml_RC_GroupVar,'GroupsToUse',{'Good','Intermediate','Poor'},'Xstep',24);
