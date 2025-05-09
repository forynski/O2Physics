OPTION="-b --configuration json://myConfig.json --pipeline track-propagation:2"

o2-analysis-tracks-extra-v002-converter ${OPTION} | 
o2-analysistutorial-mm-my-example-task-pid ${OPTION} | 
o2-analysis-track-propagation ${OPTION} | 
o2-analysis-pid-tpc ${OPTION} |
o2-analysis-pid-tpc-base ${OPTION} |
o2-analysis-pid-tof ${OPTION} |
o2-analysis-pid-tof-base ${OPTION} |
o2-analysis-pid-tof-beta ${OPTION} |
o2-analysis-timestamp ${OPTION} |
o2-analysis-event-selection ${OPTION} |
o2-analysis-multiplicity-table ${OPTION} |
o2-analysis-mccollision-converter ${OPTION}