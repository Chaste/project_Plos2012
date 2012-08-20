README.txt

 Created on: 8 May 2012
     Author: Gary Mirams
     
This project contains all the code required to recreate the results in the PloS 2012 paper.

It must be checked out into

<chaste directory>/projects/Plos2012 

in order for file paths to be picked up correctly, and it should be used with Chaste v3.1

To run the examples you should do (for example):
cd <chaste directory>
scons build=GccOpt test_suite=projects/Plos2012/test/TestCryptsAndVillusLiteratePaper.hpp

There is a Chaste wiki page associated with this project that gives links to more information on Chaste installation,
and the example simulations that are in the 'test' folder:
https://chaste.cs.ox.ac.uk/cgi-bin/trac.cgi/wiki/PaperTutorials/Plos2012
