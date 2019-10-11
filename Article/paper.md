Summary
=======

Survival analysis is a statistical tool for evaluating time-to-event
data that is widely applied across research disciplines. Commonly
reported elements of survival analysis include log-rank tests, hazard
ratios (HR) and Kaplan-Meier (KM) curves. KM-curves are used to compare
suvival durations between two groups and give users a particular
estimate of survival probability at a given time; log-rank tests are
used to conduct statistical inference on suvival durations between
groups; and HRs provide a ratio of the hazard rates between groups.
MATLAB (MATLAB 2018A) currently lacks functions to easily create
KM-plots with accompanying risk tables. Furthermore, MATLAB does not
have a built-in log-rank test, nor is one available in any of the
existing toolboxes, including the Statistics and Machine Learning
Toolbox. Our goal for MatSurv is to provide an easy-to-use tool that
creates publication quality KM-plots with corresponding risk tables. The
statistical procedures built into MatSurv can be used to compare two or
multiple groups. In addition, MatSurv allows the user to easily modify
the appearance of the created figure. The graphics were inspried by the
`survminer` R-package (Kassambara 2018).

Use
===

MatSurv uses the Mantel-Cox, sometimes called the Mantel-Haenszel,
log-rank test. Users have two options for calculating HRs: the log-rank
or Mantel-Haneszel approach. In the log-rank approach, HR =
(O<sub>a</sub>/E<sub>a</sub>)/(O<sub>b</sub>/E<sub>b</sub>), where
O<sub>a</sub> & O<sub>b</sub> are the observed events in each group and
E<sub>a</sub> & E<sub>b</sub> are the number of expected events. In the
Mantel-Haenszel approach, HR = exp((O<sub>1</sub>-E<sub>1</sub>)/V),
where O<sub>1</sub> is the number of observed events in a group,
E<sub>1</sub> is the expected number of events in the same group and V
is the total variance. Results from the log-rank approach will give
slightly different results when compared to the Mantel-Haneszel or Cox
regression approach, which is commonly used in R.

In order to use MatSurv, simply put MatSurv.m in any directory of your
choice and make sure it is added to your path. At a minimum, the user
should provide `TimeVar`, a vector with numeric time to event, either
observed or censored, `EventVar`, a vector or cell array defining events
or censoring, and `GroupVar`, a vector or cell array defining the
comparison groups (see example code below).

    [p,fh,stats]=MatSurv([], [], [], 'Xstep', 4, 'Title', 'MatSurv KM-Plot');

The function returns three pieces `p`, the log-rank p-value, `fh`, the
KM-plot figure handle, and `stats`, which are additional statistics from
the log-rank test. The user can further customize the style of their
KM-plot (line colors, labels, ticks, etc.) by making changes to the
figure handle.

When MatSurv is creating the groups based in the median value, the
default option uses values less than the median compared to all other
values, however this is a parameter that can be changed by the user.

The MatSurv software has no dependencies on toolboxes and runs
completely on base MATLAB functions.

Comparison
==========

The MatSurv output is comparable to that from `proc lifetest` in SAS and
`ggsurvplot` in R. Code for reproducing similar output in R and SAS are
shown below as well as the output from all 3 statistical programs (R,
SAS and MatSurv).

### R

    fit <- survfit(survobj ~ RISK_CYTO, data=dat)

    ggsurvplot(fit, risk.table=TRUE, pval=TRUE, risk.table.y.text.col=TRUE, risk.table.y.text=FALSE, break.time.by=24)

### SAS

    proc lifetest data=lamlv2(where=(RISK_CYTO ^= 'N.D.')) intervals=(0 to 120 by 24) timelist = (0 to 120 by 24)  plots=survival(atrisk=0 to 120 by 24 test);
    time OS_MONTHS*Surv(0);
    strata RISK_CYTO/test=logrank;
    run;

### MatSurv

    load laml_RC_data.mat

    [p,fh,stats]=MatSurv(laml_RC_TimeVar, laml_RC_EventVar,  laml_RC_GroupVar,... 'GroupsToUse', {'Good', 'Intermediate', 'Poor'}, 'Xstep', 24);

![](figure_20181022.png)

The results from MatSurv have been compared against both SAS and R and
found to return similar estimates. The Chi-Sq values and p-calues for a
long-rank test in MatSurv, SAS, and R are provided below (Table 1).

<table>
<thead>
<tr class="header">
<th>Data</th>
<th>Groups</th>
<th style="text-align: center;">MatSurv</th>
<th style="text-align: center;">MatSurv</th>
<th style="text-align: center;">SAS</th>
<th style="text-align: center;">SAS</th>
<th style="text-align: center;">Survminer</th>
<th style="text-align: center;">Survminer</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td></td>
<td></td>
<td style="text-align: center;"><strong>chi-sq</strong></td>
<td style="text-align: center;"><strong>p</strong></td>
<td style="text-align: center;"><strong>chi-sq</strong></td>
<td style="text-align: center;"><strong>p</strong></td>
<td style="text-align: center;"><strong>chi-sq</strong></td>
<td style="text-align: center;"><strong>p</strong></td>
</tr>
<tr class="even">
<td>Freireich</td>
<td>Groups</td>
<td style="text-align: center;">16.79</td>
<td style="text-align: center;">4.17E-5</td>
<td style="text-align: center;">16.79</td>
<td style="text-align: center;">4.17E-5</td>
<td style="text-align: center;">16.8</td>
<td style="text-align: center;">4.17E-5</td>
</tr>
<tr class="odd">
<td>LAML</td>
<td>RISK_CYTO</td>
<td style="text-align: center;">24.85</td>
<td style="text-align: center;">4.02E-6</td>
<td style="text-align: center;">24.85</td>
<td style="text-align: center;">&lt; 0.001</td>
<td style="text-align: center;">24.8</td>
<td style="text-align: center;">4.02E-6</td>
</tr>
<tr class="even">
<td>LAML</td>
<td>HGF Median</td>
<td style="text-align: center;">6.63</td>
<td style="text-align: center;">0.01</td>
<td style="text-align: center;">6.63</td>
<td style="text-align: center;">0.01</td>
<td style="text-align: center;">6.6</td>
<td style="text-align: center;">0.01</td>
</tr>
<tr class="odd">
<td>LAML</td>
<td>HGF Quartiles</td>
<td style="text-align: center;">13.01</td>
<td style="text-align: center;">3.09E-4</td>
<td style="text-align: center;">13.01</td>
<td style="text-align: center;">3.09E-4</td>
<td style="text-align: center;">13.0</td>
<td style="text-align: center;">30.9E-4</td>
</tr>
<tr class="even">
<td>LAML</td>
<td>HGF [6,12]</td>
<td style="text-align: center;">16.78</td>
<td style="text-align: center;">2.27 E-4</td>
<td style="text-align: center;">16.78</td>
<td style="text-align: center;">2.27E-4</td>
<td style="text-align: center;">16.8</td>
<td style="text-align: center;">2.27E-4</td>
</tr>
</tbody>
</table>

Acknowledgements
================

This project has recieved funding from

References
==========

Kassambara, A. 2018. “Survminer.” *GitHub Repository*.
<https://github.com/kassambara/survminer>; GitHub.
