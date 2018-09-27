# MatSurv
MatSurv is a simple survival analysis function for MATLAB that creates a KM-plot with risk table. Survival statistics, such as log rank p-value and Hazard Ration (HR) are also calculated. The log rank test has been evaluated to give the same results as SAS and R. The style of the KM-plot is easily changed with input parameters. No additional toolboxes are needed or depended upon. MatSurv was inspired by the [survminer R-package](https://github.com/kassambara/survminer).
The general usage is:
```matlab
[p, fh, stats] = MatSurv(TimeVar, EventVar, GroupVar,'param', value, …)
```

## Table of contents ##
- [Citing MatSurv](#Citing-MatSurv)
- [Simple Example](#Simple-Example)
- [Using MatSurv](#Using-Matsurv)
- [More Examples](#features)
- [List of all input options](#List-of-all-input-options)


## Citing MatSurv ##

Coming soon!

## Simple Example ##

The following code loads the data from "Freireich, EJ et al. 1963, Blood, 21, 699-716)" and creates a KM-plot with risk table. The time unit is weeks and the X-axis step length is changed to 4. The risk table shows how many are at risk (alive) for each time point. Censored points are marked with a vertical line. 
 
```matlab

[p,fh,stats]=MatSurv([], [], [],'Xstep',4,'Title','MatSurv KM-Plot');

```
<img src="/figures/Example_01.png" alt="MatSurv example" width="600">
 
## Using MatSurv ##

### Installation ###
Simply put MatSurv.m in any directory of your choice and make sure it is added to your path.

### Usage ###
  MatSurv(TimeVar, EventVar, GroupVar,'param', value, ...) creates a Kaplan-Meier plot with
  a risk table and calculates a log rank p-value.

  [p] = MatSurv( ... ) returns the log rank p-value
  
  [p, fh] = MatSurv( ... ) returns both p-value and figure handle
  
  [p, fh, stats] = MatSurv( ... ) returns additions stats from the log rank test
  
  [p, fh, stats] = MatSurv([], [], [], ... ) loads a test dataset from "Freireich, EJ et al. 1963, Blood, 21, 699-716"

INPUTS:
* 'TimeVar' is a vector with numeric time to event, either observed or
  censored. Values equal or less than zero will be removed by default

* 'EventVar' is a vector or cell array defining events or censored
  observation. Events are defined with a 1 and censored point with a 0. By
  default 'Dead', 'Deceased', 'Relapsed', 'Yes' are considered as events.
  'Alive', 'Living', 'Not Relapsed', 'DiseaseFree', 'No' are considered as censored.
  'EventDefinition' can be used to define other types of events. 

* 'GroupVar' is a vector or cell array defining the different groups.
  If 'GroupVar' is a continuous variable, median cut will be used as a default.

OUTPUTS:
* p       : log rank p-value
* fh      : figure handle for KM-plot figure
* stats   : Additional statistics from the log rank test


## More Examples ##

### Additional options ###
Below are some examples for how to create different styles of KM plots and also how one can make changes using the figure handle.

In the example below, we show how we can change some of the properties of the KM plot via various name-value pair arguments. 

```matlab
 
[p,fh,stats]=MatSurv([],[],[],'Xstep',4,'Title','MatSurv_KM Plot',...
'TitleOptions',{'Color','r','Interpreter','none'},'InvHR',1,...
'Xlim',32,'XMinorTick',3,'LineColor',[0 0 1;1 0 1],'LineStyle',{'-',':'},...
'LineWidth',3,'CensorLineColor','k','RT_Title','Risk Table');

```
<img src="/figures/Example_02.png" alt="MatSurv example" width="600">

### Example with multiple groups ###

This example is taken from the TCGA laml data set. Obtaining the data from cBioPortal can be found in the MatSurv/Article/MATLAB/get_laml_RC_data.m script. For this example we will load the data directly.

```matlab
load laml_RC_data.mat

[p,fh,stats]=MatSurv(laml_RC_TimeVar, laml_RC_EventVar,  laml_RC_GroupVar,...
'GroupsToUse', {'Good','Intermediate','Poor'},'Xstep',24);
```
<img src="/figures/laml_Risk_Cyto.png" alt="Multiple groups MatSurv example" width="600">

### Example with gene expression data ###

This example is also taken from the TCGA laml dataset but we also get the RNAseq gene expression data for the HGF gene. Obtaining the data from cBioPortal can be found in the MatSurv/Article/MATLAB/get_laml_HGF_gene_data.m script. For this example we will load the data directly.

```matlab
load laml_HGF_gene_data.mat

% Using median cut
[p,fh,stats]=MatSurv(laml_HGF_gene_TimeVar,laml_HGF_gene_EventVar,HGF_gene,'Xstep',12,'InvHR',1);

% Using qurtile
[p,fh,stats]=MatSurv(laml_HGF_gene_TimeVar,laml_HGF_gene_EventVar,HGF_gene,'Xstep',12,'InvHR',1,'CutPoint','quartile');

% Using Two Cut points
[p,fh,stats]=MatSurv(laml_HGF_gene_TimeVar,laml_HGF_gene_EventVar,HGF_gene,'Xstep',12,'InvHR',1,'CutPoint',[6 12]);

```
#### Median cut ####

<img src="/figures/laml_HGF_gene_Median.png" alt="Median MatSurv example" width="600">


#### Quartile ####

<img src="/figures/laml_HGF_gene_Quartile.png" alt="Quartile MatSurv example" width="600">


#### Two cut points ####

<img src="/figures/laml_HGF_gene_TwoCutPoints.png" alt="Two Cut points MatSurv example" width="600">






## List of all input options ##

* 'NoPlot': A true/false value which, if true, no figure is created
  (default: false)

* 'NoRiskTable': A true/false value which, if true, no risk table is
  included in the KM plot. (default: false)

* 'CutPoint': Either a string or scalar/vector with cut points to be used
  for defining groups based on a continuous 'GroupVar' input variable
  Allowed names are: 'Median' or 'Quartile'
  If a scalar or vector is given, the groups will be defined based on the
  cut points. (default: 'median')

* 'GroupsToUse': Cell array defining what groups to use from the GroupVar
  variable. Only works if GrouVar is a cell array. (default: all groups are used)

* 'GroupOrder': A cell array defining the group order to be used in the
  legend. (default: Groups are sorted alphabetically)

* 'EventDefinition': Two element cell array where the first cell defines
  the event and the second defines censored values. Example {'Dead,'Alive'}

* 'TimeMin': Scalar defining minimum valid time point. Subjects with time
  values below this will be removed. (default: 0)

* 'TimeMAx': Scalar value defining righ censoring time. Subjects with
  TimeVar > TimeMax will be set to TimeMax and considered as censored.
  (default: [])

* 'PairWiseP': A true/false value for calculating pairwise log rank test
  between group pairs; useful if there are more than two groups. (default: false)

* 'NoWarnings': A true/false value which, if true, no warnings are printed
  if subjects are removed. (default: false)

KM plot options
* 'LineColor': Either a matrix of size numLevels-by-3 representing the
   colormap to be used, or a string for a MATLAB colormap (lines, parula,
   cool, prism) or 'JCO' 'nejm' 'Lancet' 'Science' 'Nature','lines' for a
   set of Journal dependent palettes or my default 'aeb01' (default:'aeb01')

* 'FlipGroupOrder': Flips the order of the groups in the legend.
  (default: false)

* 'FlipColorOrder': Flips the color order of the groups.
  (default: false)

* 'KM_position': Vector defining the KM axes for the KM plot.
  (default: [0.3 0.4 0.68 0.45])

* 'RT_position': Vector defining the Risk Table axes for the KM plot.
  (default: [0.3 0.05 0.68 0.20])

* 'TimeUnit': String defining the time unit displayed on the x-axis.
  (default: 'Months')

* 'BaseFontSize': Base font size for all text in the plot.
  (default: 16)

* 'DispP': A true/false value which, if true, the log rank test p-value
  is displayed on the KM plot. (default: true)

* 'DispHR': A true/false value which, if true, the Hazard Ratio (HR)
  is displayed on the KM plot. (default: true)

* 'InvHR': A true/false value which, if true, the inverted HR value
  is displayed on the KM plot. (default: false)

* 'XLim': Vector defining the XLim. Does not affect the log rank test.
  (default: automatic)

* 'LineWidth': Scalar defining the line width used in the KM plot.
  (Default: 2)

* 'LineStyle': Cell array defining the line style for the KM plot.
  If an array is given each group will have different linestyle, for example
  'LineStyle',{'-','--',':','-.'}
  (Default: {'-'})

* 'CensorLineWidth': Scalar defining the line width of the censored ticks.
  (default: 2)

* 'CensorLineLength': Scalar defining the length of the censored ticks.
  (Default: 0.02)

* 'CensorLineColor': Text string defining color of censor ticks. 'same'
  gives the same colors as the lines while 'k' makes them black
  (Default: 'same')

* 'Xstep': Scalar defining the X tick step length.
  (defaut: automatic)

* 'XTicks': Vector defining the position of the X tick marks. 
  (Default: automatic)

* 'XMinorTick': Scalar defining the number of minor ticks between major X
  ticks. (Default: 1)

* 'Xlabel': Text string for X label (Default: 'Time(Months)')

* 'XlabelOptions': MATLAB Name-value pair arguments for X label. (Default: '')

* 'XLabelFontSize': Scalar describing X label font size change compared
  to base font size. (Default: 0)

* 'XTickFontSize': Scalar describing X tick font size change compared
  to base font size. (Default: -2)

* 'YTicks': Vector defining the position of the X tick marks.
  (Default: [0:0.2:1])

* 'YMinorTick': Scalar defining the number of minor ticks between major Y
  ticks. (Default: 1)

* 'Ylabel': Text string for Y label. (Default: 'Survival Probability' )

* 'YlabelOptions': MATLAB Name-value pair arguments for Y label. (Default: '')

* 'YLabelFontSize': Scalar describing Y label font size change compared
  to base font size. (Default: 0)

* 'YTickFontSize': Scalar describing Y tick font size change compared
  to base font size. (Default: -2)

* 'Title': Text string for Title. (Default: ' ' )

* 'TitleOptions': MATLAB Name-value pair arguments for Title. (Default: '')

* 'LegendFontSize': Scalar describing Legend font size change compared
  to base font size. (Default: -2)

* 'PvalFontSize': Scalar describing p-value font size change compared
  to base font size. (Default: 0)

Risk table plot options
* 'RT_FontSize': Scalar describing Risk Table font size change compared
  to base font size. (Default: 0)

* 'RT_Color': Text string defining color of Risk Table text. 'same'
  gives the same colors as the groups in the KM plot while 'k' would make
  them black. (Default: 'same')

* 'RT_Title': Text string for Risk Table title. (Default: '' )

* 'RT_TitleOptions': MATLAB Name-value pair arguments for Risk Table title. (Default: '')

* 'RT_YLabel': A true/false value for displaying the group names on the Risk Table
  Y axis. (Default: True )


