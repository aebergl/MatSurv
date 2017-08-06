# MatSurv
Survival analysis in MATLAB

USAGE:

  MatSurv(TimeVar, EventVar, GroupVar,'param', value, ...) creates a Kaplan-Meier plot
  
  [p] = MatSurv( ... ) returns the log rank p-value
  
  [p, fh] = MatSurv( ... ) returns both p-value and figure handle
  
  [p, fh, stats] = MatSurv( ... ) returns additions stats from log rank test
  
  [p, fh, stats] = MatSurv([], [], [], ... ) loads toy dataset

INPUTS:
* 'TimeVar' is a vector with numeric time to event either observed or
  censored. Values equal or less than zero will be removed by default

* 'EventVar' is a vector or cell array defining events or censored
  observation. Events are defined with a 1 and censored point with a 0. By
  default 'Dead', 'Deceased', 'Relapsed', 'Yes' are considered as events. 
  'Alive', 'Living', 'Not Relapsed', 'DiseaseFree', 'No' are considers as censored
  'EventDefinition' can be used to define other types of events

* 'GroupVar' is a vector or cell array that defines the different groups.
  If it is a continues variable median cut will be used as a default. 

OUTPUTS:
* p       : log rank p-value
* fh      : figure handle to KM-plot figure
* stats   : Additional statistics from the log rank test

OTHER PARAMETERS (passed as parameter-value pairs)
* 'NoPlot': A true/false value which, if true, no figure is created
  (default: false)

* 'CalcP': A true/false value which, if true, a log rank test is
  performed and displayed on the KM-plot. (default: true)

* 'CutPoint': Either a string or scalar/vector with cut points to be used
  for defining groups based on a continuous 'GroupVar' input variable
  Allowed names are: 'Median' or 'Quartile'
  If a scalar or vector is used the groups will be defined based on the
  cut points. (default: 'median')

* 'GroupsToUse': Cell array defining what groups to use from the GroupVar
  variable. Works only if GrouVar is a cell array. (default: all groups are used)

* 'GroupOrder': A cell array defining the group order to be used in the
  legend. (default: Groups are sorted alphabetically)

* 'EventDefinition': Two element cell array where the first cell defines
  the event and the second censored values. Example {'Dead,'Alive'}
  (default = true)

* 'FlipGroupOrder': Flips the order of the groups in the legend. 
  default: false)

* 'FlipColorOrder': Flips the color order of the groups. 
  (default: false)

* 'BaseFontSize': Base font size for all text in the plot
  (default: 16)

* 'KM_position': Vector defining the KM axes for the KM plot
  (default: [0.2 0.4 0.75 0.45])

* 'RT_position': Vector defining the Risk Table axes for the KM plot
  (default: [0.2 0.05 0.75 0.20])

KM plot options
* 'RT_position': Vector defining the Risk Table axes for the KM plot
  (default: [0.2 0.05 0.75 0.20]

* 'XLim': Vector defining the XLim. Do not affects the log rank test
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

* 'CensorLineLength':Scalar defining the length of the censored ticks
  (Default: 0.02)

* 'CensorLineColor': Text string defining color of censor ticks. 'same
  gives the same colors as the lines while 'k' would make them black
  (Default: 'same')

* 'Xstep': Scalar defining the X tick step length. 
  (defaut: automatic)

* 'XTicks': Vector defining the position of the X-tick marks
  (Default: 16)

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

  EXAMPLES:

*** Anders Berglund ***

