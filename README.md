# Applied Python Tutorial - Ando Lab
This is a simple example of how we can use python to perform some basic and some semi-fancy analysis.
Note: This tutorial will attempt to touch on all basic concepts covered in Darren Xu's lecture (Available upon request)

**WARNING:**

If you choose to modify this code you can always go back and access the original code here on Github. Currently, only Rob Miller can modify the 'master branch'.

_For example:_
**Potential User Issues**   
(1) Attempting to run code with Python2 - you will get lots of print and other errors.  
(2) Non-Anaconda Python may result in missing packages - use 'pip' to install necessary modules  
_for example:_ 
```shell script
pip3 install matplotlib  
pip3 install numpy  
pip3 install scipy  
```
(3) 
     
## How The Code Works
A theoritical scattering curve was determined from the crystal structure, PDBID=6MT9, via the software FoxS. I've provided the scattering profile in this repository.

From this theoritical profile, we will learn how to:  
(1) Import data  
(2) Import modules/classes  
(3) Define functions/objects
(4) Review conditionals and loops
(5) Perform least squares minimization and Guiner analysis
(6) Calculate a pair-distance distribution function
(7) Exporting data & so much more! 

**Details of the backend of the code:**

1. For the calculation for the pair-wise distance distribution function, see: !REF


_I turned off RuntimeWarnings caused by division by 0. Beware!!_
```python
import warnings

warnings.filterwarnings("ignore", category=RuntimeWarning) 
```

Good luck. 

Rob


## To access Spyder via anaconda
1. Open terminal & issue the following command

 ```shell script 
 cd
 spyder
 ```

*** Make sure you have python3 as your default interpreter in Spyder
