cgdsURL = 'http://www.cbioportal.org/public-portal/';
try
    c = getclinicaldata(cgdsURL,'laml_tcga_pub_all');
catch
    disp('Get MSKCC CGDS Cancer Genomics Toolbox from http://www.mathworks.com/matlabcentral/fileexchange/31297-mskcc-cgds-cancer-genomics-toolbox')
    return
end

laml_HGF_gene_TimeVar = c.data(:,strcmp('OS_MONTHS',c.clinVariable));
laml_HGF_gene_EventVar = c.data(:,strcmp('OS_STATUS',c.clinVariable));
EmptyIndx = cellfun('isempty',laml_HGF_gene_TimeVar);
if any(EmptyIndx)
    laml_HGF_gene_TimeVar(EmptyIndx)={'NaN'};
end
laml_HGF_gene_TimeVar = sscanf(sprintf('%s,',laml_HGF_gene_TimeVar{:}),'%f,');

Clinical_SampleId = c.caseId;
HGF_RNAseq = getprofiledata(cgdsURL, 'laml_tcga_pub_all','laml_tcga_pub_rna_seq_v2_mrna','HGF', true);

[~,indx_surv,indx_rnaseq]=intersect(Clinical_SampleId,HGF_RNAseq.caseId,'stable');

HGF_gene = HGF_RNAseq.data(indx_rnaseq);
laml_HGF_gene_TimeVar = laml_HGF_gene_TimeVar(indx_surv);
laml_HGF_gene_EventVar = laml_HGF_gene_EventVar(indx_surv);

HGF_gene = log2(HGF_gene + 1);
nan_indx = isnan(HGF_gene);
HGF_gene(nan_indx) = [];
laml_HGF_gene_TimeVar(nan_indx) = [];
laml_HGF_gene_EventVar(nan_indx) = [];

clear c EmptyIndx cgdsURL  Clinical_SampleId HGF_RNAseq nan_indx indx_surv indx_rnaseq
[p,fh,stats]=MatSurv(laml_HGF_gene_TimeVar,laml_HGF_gene_EventVar,HGF_gene,'Xstep',12,'InvHR',1);
[p,fh,stats]=MatSurv(laml_HGF_gene_TimeVar,laml_HGF_gene_EventVar,HGF_gene,'Xstep',12,'InvHR',1,'CutPoint','quartile');
[p,fh,stats]=MatSurv(laml_HGF_gene_TimeVar,laml_HGF_gene_EventVar,HGF_gene,'Xstep',12,'InvHR',1,'CutPoint',[6 12]);
