#!/usr/bin/python
import os
import string
import sys

os.system("/bin/rm -rf DICE_2.0.tar.gz dice/src/*.o dice/src/*.so")
os.system("R CMD build dice")
os.system("R CMD INSTALL dice")

## For a local installation - comment the above line and uncomment the lines starting with 'mkdir', 'export' and 'os.system'
## Create a directory in the home directory, called for example R_libs:
## mkdir /home/your_username/R_libs
## set a variable to point R at that directory: 
## export R_LIBS="/home/your_username/R_libs"
## and then install the package ih this directory:
#os.system("R CMD INSTALL $HOME/R_LIBS DICE_2.0.tar.gz")




# If the '-test' flag is used, run the test script
if len(sys.argv) > 1:
   if any(x=="-test" for x in sys.argv):
      print("\nDICE build complete. Begin testing.\n")
      os.system("Rscript test/TestChanges.R")
      
