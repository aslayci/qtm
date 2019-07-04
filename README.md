# qtm
Query the model

# Compilation  

**make -f Makefile_preprocessing
**Make -f Makefile_queryProcessing

# Configuration of input files and parameters 

**config.txt:** sample configuration file that needs to be edited and passed as a command line argument.  
It contains all the necessary information on arguments required as input.  

# Running from command line

**materialization:** ./main_preprocess -c config.txt 

**query processing:** ./main_queryProcess -c config.txt 



