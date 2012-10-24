# Set the task name
TASK = fid_drift_mon

# Uncomment the correct choice indicating either SKA or TST flight environment
FLIGHT_ENV = SKA

# Set the names of all files that get installed
SHARE = fid_drift_mon.pl plot_drift.py plot_drift_model.py
DATA = task_schedule.cfg

include /proj/sot/ska/include/Makefile.FLIGHT

install:
#  Uncomment the lines which apply for this task
	mkdir -p $(INSTALL_SHARE)
	rsync --times --cvs-exclude $(SHARE) $(INSTALL_SHARE)/
	rsync --times --cvs-exclude $(DATA) $(INSTALL_DATA)/
