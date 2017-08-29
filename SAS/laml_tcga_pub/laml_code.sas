* data downloaded from Cbioportal at http://www.cbioportal.org/study?id=laml_tcga_pub#summary;
* data was imported using SAS Import Wizard;

data lamlv2;
set laml;
where Risk__Cyto_ ^= "N.D.";
if Overall_Survival_Status = "LIVING" then Surv=0;
else if Overall_Survival_Status = "DECEASED" then Surv=1;
else Surv=.;
run;

proc lifetest data=lamlv2(where=(Risk__Cyto_ ^= 'N.D.')) plots=survival(atrisk test);
time Overall_Survival__Months_*Surv(0);
strata Risk__Cyto_/test=logrank;
run;
