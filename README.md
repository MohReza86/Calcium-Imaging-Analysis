# Calcium Imaging Analysis 

To identify calcium events, we initially de-trended raw calcium traces by subtracting a 
local  minimum  value  from  each  time  point.  The  local  minimum  was  defined  as  the 
minimum  value  within  a  13  s  period  of  the  trace.  We  then  extracted  fluorescence 
changes over time (ΔF/F) as follows:

ΔF/F = (F-F0)/F0

where F is the fluorescent value at a given time and  F0 is the mean fluorescent value in 
the non-event portion of the trace, defined as the the minimal averaged signal across 8 
s  periods throughout  the  measurement.  The  de-trended ΔF/F  calcium  traces  were 
further  smoothed  using  a  moving  average  of  0.8  s,  scaled  by  a  value  of  1.2.  Due  to 
motion  artifact,  the  first  10  s  of  the  fluorescence  traces  was  removed.  All  points 
separated by at least 1.5 s with fluorescence signal exceeding five standard deviations 
above the baseline, calculated as the mean of the non-event portion of the trace, were 
identified as calcium events. 

