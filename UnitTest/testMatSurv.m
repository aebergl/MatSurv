%  MatSurv test

tol = 1e-8;

% preconditions
assert(exist('MatSurv','file') == 2,'Could not find MatSurv function')

%% Test 1: Freireich data
try 
    load Freireich_Results
catch
    error('Could not load Freireich_Results data')
end

try
    [p,fh,stats] = MatSurv([],[],[],'NoPlot',false,'Print',false);
catch
    error('Could not run MatSurv')
end

assert(strcmp(get(fh, 'type'), 'figure') & isa(fh, 'matlab.ui.Figure') );

assert(abs(p - p_true) < tol)
assert(abs(stats.p_MC - stats_true.p_MC) < tol)
assert(abs(stats.Chi2_MC - stats_true.Chi2_MC) < tol)
assert(abs(stats.HR_logrank - stats_true.HR_logrank) < tol)
assert(abs(stats.HR_logrank_Inv - stats_true.HR_logrank_Inv) < tol)
assert(abs(stats.HR_MH - stats_true.HR_MH) < tol)
assert(abs(stats.HR_MH_Inv - stats_true.HR_MH_Inv) < tol)

close(fh);

 
%% Test 2: laml_RC_data default setting
try 
    load laml_RC_data
catch
    error('Could not load laml_RC_data')
end

try 
    load laml_RC_Results_01
catch
    error('Could not load laml_RC_Results_01 data')
end

try
    [p,fh,stats]=MatSurv(laml_RC_TimeVar,laml_RC_EventVar, laml_RC_GroupVar,...
        'GroupsToUse', {'Good','Intermediate','Poor'},'Xstep',24,'Print',false);
catch
    error('Could not run MatSurv')
end

assert(strcmp(get(fh, 'type'), 'figure') & isa(fh, 'matlab.ui.Figure') );

assert(abs(p - p_true) < tol)
assert(abs(stats.p_MC - stats_true.p_MC) < tol)
assert(abs(stats.Chi2_MC - stats_true.Chi2_MC) < tol)

close(fh);

%% Test 3: laml_RC_data merged groups and TimeMax=60
try 
    load laml_RC_data
catch
    error('Could not load laml_RC_data')
end

try 
    load laml_RC_Results_02
catch
    error('Could not load laml_RC_Results_02 data')
end

try
     [p,fh,stats]=MatSurv(laml_RC_TimeVar,laml_RC_EventVar, laml_RC_GroupVar,...
         'GroupsToUse', {'Good','Intermediate',{'ND+Poor','Poor','N.D.'}},...
         'Xstep',24,'Print',false,'Xstep',12,'TimeMax',60);
catch
    error('Could not run MatSurv')
end

assert(strcmp(get(fh, 'type'), 'figure') & isa(fh, 'matlab.ui.Figure') );

assert(abs(p - p_true) < tol)
assert(abs(stats.p_MC - stats_true.p_MC) < tol)
assert(abs(stats.Chi2_MC - stats_true.Chi2_MC) < tol)

close(fh);

%% Test 4: laml_RC_data pairwise comparisons
try 
    load laml_RC_data
catch
    error('Could not load laml_RC_data')
end

try 
    load laml_RC_Results_03
catch
    error('Could not load laml_RC_Results_03 data')
end

try
     [p,fh,stats]=MatSurv(laml_RC_TimeVar,laml_RC_EventVar, laml_RC_GroupVar,...
         'PairwiseP',true,'Xstep',24,'Print',false,'Xstep',12,...
         'GroupsToUse', {'Good','Intermediate','Poor'});
catch
    error('Could not run MatSurv')
end

assert(strcmp(get(fh, 'type'), 'figure') & isa(fh, 'matlab.ui.Figure') );

assert(abs(p - p_true) < tol)
assert(abs(stats.p_MC - stats_true.p_MC) < tol)
assert(abs(stats.Chi2_MC - stats_true.Chi2_MC) < tol)

for i=1:3
    assert(abs(p - p_true) < tol)
    assert(abs(stats.ParwiseStats(i).p_MC - stats_true.ParwiseStats(i).p_MC) < tol)
    assert(abs(stats.ParwiseStats(i).Chi2_MC - stats_true.ParwiseStats(i).Chi2_MC) < tol)
    assert(abs(stats.ParwiseStats(i).HR_logrank - stats_true.ParwiseStats(i).HR_logrank) < tol)
    assert(abs(stats.ParwiseStats(i).HR_logrank_Inv - stats_true.ParwiseStats(i).HR_logrank_Inv) < tol)
    assert(abs(stats.ParwiseStats(i).HR_MH - stats_true.ParwiseStats(i).HR_MH) < tol)
    assert(abs(stats.ParwiseStats(i).HR_MH_Inv - stats_true.ParwiseStats(i).HR_MH_Inv) < tol) 
end


close(fh);

%% Test 5: laml_HGF_gene_data median cut
try 
    load laml_HGF_gene_data.mat
catch
    error('Could not load laml_HGF_gene_data.mat')
end

try 
    load laml_HGF_gene_Results_01
catch
    error('Could not load laml_HGF_gene_Results_01 data')
end

try
    [p,fh,stats]=MatSurv(laml_HGF_gene_TimeVar,laml_HGF_gene_EventVar,HGF_gene,...
        'Xstep',12,'InvHR',1,'Print',false); 
catch
    error('Could not run MatSurv')
end

assert(strcmp(get(fh, 'type'), 'figure') & isa(fh, 'matlab.ui.Figure') );

assert(abs(p - p_true) < tol)
assert(abs(stats.p_MC - stats_true.p_MC) < tol)
assert(abs(stats.Chi2_MC - stats_true.Chi2_MC) < tol)
assert(abs(stats.HR_logrank - stats_true.HR_logrank) < tol)
assert(abs(stats.HR_logrank_Inv - stats_true.HR_logrank_Inv) < tol)
assert(abs(stats.HR_MH - stats_true.HR_MH) < tol)
assert(abs(stats.HR_MH_Inv - stats_true.HR_MH_Inv) < tol)

close(fh);

%% Test 6: laml_HGF_gene_data tertile
try 
    load laml_HGF_gene_data.mat
catch
    error('Could not load laml_HGF_gene_data.mat')
end

try 
    load laml_HGF_gene_Results_02
catch
    error('Could not load laml_HGF_gene_Results_02 data')
end

try
    [p,fh,stats]=MatSurv(laml_HGF_gene_TimeVar,laml_HGF_gene_EventVar,HGF_gene,...
    'CutPoint','tertile','Xstep',12,'InvHR',1,'GroupOrder',[2 1 3],'Print',false); 
 
catch
    error('Could not run MatSurv')
end

assert(strcmp(get(fh, 'type'), 'figure') & isa(fh, 'matlab.ui.Figure') );
close(fh);

assert(abs(p - p_true) < tol)
assert(abs(stats.p_MC - stats_true.p_MC) < tol)
assert(abs(stats.Chi2_MC - stats_true.Chi2_MC) < tol)


%% Test 7: laml_HGF_gene_data tertile
try 
    load laml_HGF_gene_data.mat
catch
    error('Could not load laml_HGF_gene_data.mat')
end

try 
    load laml_HGF_gene_Results_03
catch
    error('Could not load laml_HGF_gene_Results_03 data')
end

try
    [p,fh,stats]=MatSurv(laml_HGF_gene_TimeVar,laml_HGF_gene_EventVar,HGF_gene,...
    'CutPoint',[6 12],'Xstep',12,'InvHR',1,'RT_KMplot',true,'Print',false); 
 
catch
    error('Could not run MatSurv')
end

assert(strcmp(get(fh, 'type'), 'figure') & isa(fh, 'matlab.ui.Figure') );

assert(abs(p - p_true) < tol)
assert(abs(stats.p_MC - stats_true.p_MC) < tol)
assert(abs(stats.Chi2_MC - stats_true.Chi2_MC) < tol)

close(fh);