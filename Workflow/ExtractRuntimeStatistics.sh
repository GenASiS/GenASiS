#!/bin/bash

#-- This script extracts value of interest from ShowStatistics output 
#   and write them to files compatible for plotting with Jenkins

OUTPUT=$1
ID=$2

tac $OUTPUT | grep -m 1 "Elapsed wall time" \
  | awk -F'=' '{print $2}' | awk '{print "YVALUE=",$1}' \
  | sed 's/ //g' > walltime.$ID

tac $OUTPUT | grep -m 1 "Computational" \
  | awk -F'=' '{print $2}' | awk '{print "YVALUE=",$1}' \
  | sed 's/ //g' > compute.$ID

tac $OUTPUT | grep -m 1 "Input/output" \
  | awk -F'=' '{print $2}' | awk '{print "YVALUE=",$1}' \
  | sed 's/ //g' > io.$ID

tac $OUTPUT | grep -m 1 "Across processes max HWM" \
  | awk -F'=' '{print $2}' | awk '{print "YVALUE=",$1}' \
  | sed 's/ //g' > memory.$ID
