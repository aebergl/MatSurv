# MatSurv
MatSurv is a simple survival analysis function for MATLAB that creates a KM-plot with risk table. Survival statistics, such as log rank p-value Hazard Ration (HR) are also calculated. The log rank test has been tested to give the same results as SAS and R. The style of the KM-plot is easily changed with input parameters. No additional toolboxes are needed. MatSurv was inspired by the [survminer R-package](https://github.com/kassambara/survminer).
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

The following code loads the data from "Freireich, EJ et al. 1963, Blood, 21, 699-716)" and creates a KM-plot with risk table. The time unit is weeks and the X-axis step length is changed to 4. The risk table shows how many is at risk (alive) for each time point. Censored points are marked with a vertical line. 
 
```matlab

[p,fh,stats]=MatSurv([], [], [],'Xstep',4,'Title','MatSurv KM-Plot');

```

<img src="/html/Example_01.png" alt="MatSurv example" width="800">

## Using MatSurv ##

### Installation ###
Simply put MatSurv.m in any directory of your choice and make sure its added to your path

### Usage ###
  MatSurv(TimeVar, EventVar, GroupVar,'param', value, ...) creates a Kaplan-Meier plot,
  a risk table and calculates a log rank p-value

  [p] = MatSurv( ... ) returns the log rank p-value
  
  [p, fh] = MatSurv( ... ) returns both p-value and figure handle
  
  [p, fh, stats] = MatSurv( ... ) returns additions stats from log rank test
  
  [p, fh, stats] = MatSurv([], [], [], ... ) loads test dataset

INPUTS:
* 'TimeVar' is a vector with numeric time to event, either observed or
  censored. Values equal or less than zero will be removed by default

* 'EventVar' is a vector or cell array defining events or censored
  observation. Events are defined with a 1 and censored point with a 0. By
  default 'Dead', 'Deceased', 'Relapsed', 'Yes' are considered as events.
  'Alive', 'Living', 'Not Relapsed', 'DiseaseFree', 'No' are considers as censored
  'EventDefinition' can be used to define other types of events

* 'GroupVar' is a vector or cell array defining the different groups.
  If 'GroupVar' is a continues variable, median cut will be used as a default.

OUTPUTS:
* p       : log rank p-value
* fh      : figure handle to KM-plot figure
* stats   : Additional statistics from the log rank test


## More Features ##

Below is some examples given how to create different style of the KM-Plot and also how one can make change using the figure handle

## List of all input options ##

* 'NoPlot': A true/false value which, if true, no figure is created
  (default: false)

* 'NoRiskTable': A true/false value which, if true, no risk table is
  included in the KM-plot. (default: false)

* 'CutPoint': Either a string or scalar/vector with cut points to be used
  for defining groups based on a continuous 'GroupVar' input variable
  Allowed names are: 'Median' or 'Quartile'
  If a scalar or vector is given the groups will be defined based on the
  cut points. (default: 'median')

* 'GroupsToUse': Cell array defining what groups to use from the GroupVar
  variable. Works only if GrouVar is a cell array. (default: all groups are used)

* 'GroupOrder': A cell array defining the group order to be used in the
  legend. (default: Groups are sorted alphabetically)

* 'EventDefinition': Two element cell array where the first cell defines
  the event and the second censored values. Example {'Dead,'Alive'}

* 'TimeMin': Scalar defining minimum valid time point. Subjects with time
  values below this will be removed. (default: 0)

* 'TimeMAx': Scalar value defining righ censoring time. Subjects with
  TimeVar > TimeMax will be set to TimeMax and considered as censored.
  (default: [])

* 'PairWiseP': A true/false for caulculating pairwise log rank test
  between group pairs, useful if there is more than two groups. (default: false)

* 'NoWarnings': A true/false value which, if true, no warnings are printed
  id subjects are removed. (default: false)

KM plot options
* 'FlipGroupOrder': Flips the order of the groups in the legend.
  (default: false)

* 'FlipColorOrder': Flips the color order of the groups.
  (default: false)

* 'KM_position': Vector defining the KM axes for the KM plot
  (default: [0.3 0.4 0.68 0.45])

* 'RT_position': Vector defining the Risk Table axes for the KM plot
  (default: [0.3 0.05 0.68 0.20])

* 'TimeUnit': String defning time unit displayd on the x-axis.
  (default: 'Months')

* 'BaseFontSize': Base font size for all text in the plot
  (default: 16)

* 'DispP': A true/false value which, if true, log rank test p-value
  is displayed on the KM-plot. (default: true)

* 'DispHR': A true/false value which, if true, Hazard ration (HR)
  is displayed on the KM-plot. (default: true)

* 'InvHR': A true/false value which, if true, the iverted HR value
  is displayed on the KM-plot. (default: false)

* 'XLim': Vector defining the XLim. Do not affect the log rank test
  (default: automatic)

* 'LineColor': Either a matrix of size numLevels-by-3 representing the
  colormap to be used or a string for a MATLAB colormap
  (default: 'lines')

* 'LineWidth': Scalar defining the line width used in the KM-plot
  (Default: 2)

* 'LineStyle': Cell array defining the linestyle for the KM-plot.
  If an array is given each group will have different linestyle, for example
  'LineStyle',{'-','--',':','-.'}
  (Default: {'-'})

* 'CensorLineWidth': Scalar defining the linewith of the censored ticks
  (default: 2)

* 'CensorLineLength': Scalar defining the length of the censored ticks
  (Default: 0.02)

* 'CensorLineColor': Text string defining color of censor ticks. 'same
  gives the same colors as the lines while 'k' would make them black
  (Default: 'same')

* 'Xstep': Scalar defining the X tick step length.
  (defaut: automatic)

* 'XTicks': Vector defining the position of the X-tick marks
  (Default: automatic)

* 'XMinorTick': Scalar defining the number of minor ticks between major X
  ticks (Default: 1)

* 'Xlabel': Text string for X-label (Default: 'Time(Months)' )

* 'XlabelOptions': MATLAB Name-value pair arguments for xlabel (Default: '')

* 'XLabelFontSize': Scalar describing Xlabel font size change compared
  to base font size (Default: 0)

* 'XTickFontSize': Scalar describing Xtick font size change compared
  to base font size (Default: -2)

* 'YTicks': Vector defining the position of the X-tick marks
  (Default: [0:0.2:1])

* 'YMinorTick': Scalar defining the number of minor ticks between major Y
  ticks (Default: 1)

* 'Ylabel': Text string for Y-label (Default: 'Survival Probability' )

* 'YlabelOptions': MATLAB Name-value pair arguments for ylabel (Default: '')

* 'YLabelFontSize': Scalar describing Ylabel font size change compared
  to base font size (Default: 0)

* 'YTickFontSize': Scalar describing Ytick font size change compared
  to base font size (Default: -2)

* 'Title': Text string for Title (Default: '' )

* 'TitleOptions': MATLAB Name-value pair arguments for Title (Default: '')

* 'LegendFontSize': Scalar describing Legend font size change compared
  to base font size (Default: -2)

* 'PvalFontSize': Scalar describing p-value font size change compared
  to base font size (Default: 0)

Risk table plot options
* 'RT_FontSize': Scalar describing Risk Table font size change compared
  to base font size (Default: 0)

* 'RT_Color': Text string defining color of Risk table text. 'same
  gives the same colors as the groups in the KM plot while 'k' would make
  them black (Default: 'same')

* 'RT_Title': Text string for Risk Table Title (Default: '' )

* 'RT_TitleOptions': MATLAB Name-value pair arguments for Risk Table Titel (Default: '')

* 'RT_YLabel': True/False for displaying the group names on the Risk table
  Y-axis (Default: True )


