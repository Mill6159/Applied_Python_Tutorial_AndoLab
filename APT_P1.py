#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  3 17:32:11 2020

@author: Rob
"""

'''
We will build a script that imports FoxS generated SAXS profile from
a crystal structure and does some analysis. PDBID: XXXX. Familiar?

Learning Objectives:
    (1) How to load in modules (& install them if necessary):
        (a) pip or conda
    (2) How to import .txt files
    (3) Assign data to objects
    (4) Transform data and check data types
    (5) Build a useful plot class()
    (6) Basic function building/analysis - Least squares minimization (LSM)
    (7) Pair distance distribution function
    (8) Data export
'''


'''
Start with importing packages:
    (1) numpy: basic data structure manipulation, mathematics, etc
    (2) matplotlib: Plot module (PlotClass)
    (3) scipy: Great data science module
    (4) Our own class!
'''

import numpy as np # this is a math package
from scipy.optimize import curve_fit as cf
from PlotClass import *
from matplotlib import pyplot as plt
import csv
import warnings

warnings.filterwarnings("ignore", category=RuntimeWarning) # deals with dividing through by zero but be careful!!


'''
Import FoxS Data: I generated this for us from the crystal structure PDBID: 6MT9 
'''

data=np.loadtxt('6mt9.pdb.dat',
                dtype={'names': ('q', 'I(q)','I(q)_Error'), 'formats': (np.float,np.float,np.float)},
                comments='#')


'''
How can we access the data?
#1 - the 'head' of the data(i.e. first 5 lines)?
'''
N=5
# Method 1:
filename='6mt9.pdb.dat'
with open(filename, "r") as f: # important to use with command when dealing with files
    counter = 0
    print('File: %s' % filename)
    for line in f:
        print (line)
        counter += 1
        if counter == N: break

# Method 2:
with open(filename, "r") as f:
    print('File: %s' % filename)
    for i in range(N):
        line = next(f)
        print (line)
        
# Method 3: Can you think of another one?

#N = 10
#with open("file.txt", "a") as file:  # the a opens it in append mode
#    for i in range(N):
#        line = next(file).strip()
#        print(line)

## Mention stack exchange
'''
But... we used a convenient python module 'numpy' to open the file
So life is even easier but we may have lost information..
Observe
'''

print('First 10 rows of the dataframe (No column labels): ',data[0:10])

'''
We have lost the header information, beware! If you are doing more advanced data processing and require
header file information in order to process the data, you will need to be more clever 
than the simple import we are using.
But.. There are some benefits. We can now access each column by its 'name'
'''

print('The first 5 I(q) points: ',data['I(q)'][0:5])


'''
Now before we do any analysis, we need to define a few functions:
    (1) A basic linear model function. We will pass this function into a module function for least squares minimization
    (2) We will define our own function for least squares regression
    
First we will start with a basic linear model with a slope and an intercept
'''

print('-------------------------------------------------------------------------------')
print('')

def lineModel(x,m,b):
    return m*x+b

'''
Now our function that will both fit the data and extract the coefficients/error (LSM)
'''

def lsq_w_sigma(X, Y, SIG):

    """ This is a special version of linear least squares that uses 
    known sigma values of Y to calculate standard error in slope and 
    intercept. The actual fit is also SIGMA-weighted"""

    XSIG2 = np.sum(X / (SIG ** 2))
    X2SIG2 = np.sum((X ** 2) / (SIG ** 2))
    INVSIG2 = np.sum(1.0 / (SIG ** 2))

    XYSIG2 = np.sum(X * Y / (SIG ** 2))
    YSIG2 = np.sum(Y / (SIG ** 2))

    delta = INVSIG2 * X2SIG2 - XSIG2 ** 2

    sby = X2SIG2 / delta
    smy = INVSIG2 / delta

    by = (X2SIG2 * YSIG2 - XSIG2 * XYSIG2) / delta
    my = (INVSIG2 * XYSIG2 - XSIG2 * YSIG2) / delta
    # outputs slope(my), intercept(by), sigma slope(smy), sigma intercept (sby)
    return my, by, smy, sby



'''
Recall, the ultimate goal is to perform a Guiner analysis of the theoritical scattering curve
Guiner approximation: ln(I(q)) = m*x+b
where the slope(m)=sqrt(-3b), x = q^2, & the intercept(b)=I(0)

We will use scipy and the lineModel to determine Rg & compare the results to our 'homebaked' function
For this tutorial, we will not discuss automated methods for determining the optimal q-range
Instead I will provide it (nmin,nmax): (0,120)
Recall, ideally, ~(qmin,qmax): (0.3,1.3) but this is empirically determined (i.e. a suggestion more than a rule).
'''

nmin=8
nmax=34

'''
Note, LSM algorithms can take input "guesses"(g) as a starting point for minimization
The package curve_fit from scipy will output coefficients (c) and covariance (cov) matrices separately

Nuisance points:
    (1) Several data transformations will be performed:
        (a) Transformation of the signal errors in to log-space
        (b) Transformation of the slope to calculate Rg and the error, transforming back out of log space
'''
g=[-40,10] # slope, intercept

c,cov = cf(lineModel,data['q'][nmin:nmax]**2,np.log(data['I(q)'][nmin:nmax]),g,
           sigma=data['I(q)_Error'][nmin:nmax]/data['I(q)'][nmin:nmax])

scipyI0=np.exp(c[1])
scipyRg = np.sqrt(-3*c[0])
scipyRg_Err = np.absolute(scipyRg * np.sqrt(cov[0,0])/(2*c[0]))
scipy_qminRg = scipyRg*data['q'][nmin]
scipy_qmaxRg = scipyRg*data['q'][nmax]

'''
Perform the analysis with our "homebaked" function
'''
hb_slope,hb_inter,hb_sigslope,hb_siginter = lsq_w_sigma(data['q'][nmin:nmax]**2,
                                                        np.log(data['I(q)'][nmin:nmax]), 
                                                        data['I(q)_Error'][nmin:nmax]/data['I(q)'][nmin:nmax])
hbI0=np.exp(hb_inter)
hbRg = np.sqrt(-3*hb_slope)
hbRg_Err = np.absolute(hbRg* np.sqrt(hb_sigslope)/(2*hb_slope))
hb_qminRg = hbRg*data['q'][nmin]
hb_qmaxRg = hbRg*data['q'][nmax]

'''
We will create a dictionary of all of the useful data values
& review how convenient dictionaries are

Extract the following information:
    (1) Rg/Relative Error in Rg
    (2) qminRg,qmaxRg (acceptable Guiner region verification)
    (3) Comparison of the two results


Note: We need to convert our results in order to determine Rg and the relative error in Rg
For references:
    (1) https://github.com/Mill6159/SAXSProf_Desktop_GUI/blob/master/Analytical_Derivation_Final.pdf
    (2) Numerical Recipes: Press, William H., et al. Numerical recipes 3rd edition: The art of scientific computing. Cambridge university press, 2007.
'''
scipy_dict = {'scipy_qminRg':scipy_qminRg,'scipy_qmaxRg':scipy_qmaxRg,
              'scipyRg':scipyRg,'scipyRg_Err':scipyRg_Err}

#print(scipy_dict.keys()) # very useful property of dictionaries
#print(scipy_dict.values())

'''
Example of how to clean up the output a bit
'''

print('from scipy fit, Rg and Rg Error: %.2f %.8f'%(scipy_dict['scipyRg'],
                                                             scipy_dict['scipyRg_Err']))

hb_dict = {'hb_qminRg':hb_qminRg,'hb_qmaxRg':hb_qmaxRg,'hbRg':hbRg,'hbRg_Err':hbRg_Err}
#print(hb_dict.keys()) # very useful property of dictionaries
#print(hb_dict.values())

'''
Example of how to clean up the output a bit
'''

print('from hb fit, Rg and Rg Error: %.2f %.8f'%(hb_dict['hbRg'],
                                                             hb_dict['hbRg_Err']))

print('qminRg,qmaxRg: %.4f %.4f'%(hb_dict['hb_qminRg'],hb_dict['hb_qmaxRg']))

'''
Generate an automated output to "warn" us if the relative error is 
greater than 1%
'''

if (100*(hb_dict['hbRg_Err']/hbRg)) > 1.0:
    print('Percent relative error in Rg is greater than 1%:' + ' ' + '%.2f'%(100*(hb_dict['hbRg_Err']/hbRg))+'%')
    print('Review quality of fit')
else:
    print('Percent relative error in Rg is less than 1%:' + ' ' + '%.2f'%(100*(hb_dict['hbRg_Err']/hbRg))+'%')
    print(' ')

'''
Stop!
We need to build a plot class, that you can save and use later, that will generate some beautiful plots for us
Name the file: PlotClass.py

Now we can use our beautiful plot class to generate a few plots to visualize our data
'''

##### End lecture 1 #####

visuals = PlotClass() # assign class to an object. The object will now have all characteristics of the class (def __init__(self):)


visuals.semilogyPlot(data['q'],data['I(q)'],
                  plotlabel='Scattering Profile',savelabel='Example_ScatteringCurve',
                  xlabel='q ($\AA^{-1}$)',
                  ylabel='I(q)')

'''
Compare theoritical Guiner region to Guiner model
'''

visuals.twoPlot(X=data['q'][nmin:nmax]**2,Y1=lineModel(data['q'][nmin:nmax]**2,hb_slope,hb_inter),
                  Y2=np.log(data['I(q)'][nmin:nmax]),
                  plotlabel1='Model',plotlabel2='Theoretical Data',
                  savelabel='Example_Guiner2',
                  xlabel='q$^{2}$ ($\AA^{-2}$)',
                  ylabel='ln(I(q))')


'''
Look at the model in reciprocal space (rather than log reciprocal space)
How to convert back?
'''


def gaussianGuiner(X,I0,Rg):
    return I0*np.exp(-(1/3)*(X**2)*(Rg**2))
    
gaussModel=gaussianGuiner(data['q'][nmin:nmax], hbI0, hbRg)


visuals.twoPlot(X=data['q'][nmin:nmax],Y1=gaussModel,
                  Y2=data['I(q)'][nmin:nmax],
                  plotlabel1='Model',plotlabel2='Theoretical Data',
                  savelabel='Example_Guiner3',
                  xlabel='q ($\AA^{-1}$)',
                  ylabel='I(q)')

##############################################################################

'''
Pair distance distribution calculation:
    (1) We will define a new function for this "PDDF":
        (a) I left the potential capability to add other shape profiles
        where you will have to deal with properly interpolating q onto I(q)
'''

def PDDF(shape,Dmax,I,q):
    """plot pair distance distribution"""
    # reference:"Svergen&Koch,Rep.Phys.Prog 66(2003) 1735-82
    r_range=np.arange(0,Dmax*1.4,Dmax/50)
    if shape == "FoxS":
        if q[0] == 0:
            q = q[1:]
            I = I[1:]
            # Taking the first point in exp_q out if it's 0, avoiding dividing by 0 problem
        else:
            q = q

    '''
    Initialize an array to stick the data in
    Note: this MUST be outside of the loop, otherwise we will clear the array for each iteration of the loop.
    '''
    P_r = np.array([], dtype=float)
    
    '''
    Can you think of a better, or different, way to write this loop?
    '''
    
    for r in r_range:
        p_r = np.sum(q ** 2 * I * np.sin(q * r) / (q * r) * 0.02) * (r ** 2) / (2.0 * np.pi ** 2)
        P_r = np.append(P_r, p_r)

    return P_r,r_range


Pr,r=PDDF(shape='FoxS',Dmax=80,I=data['I(q)'],q=data['q'])

for entry in Pr:
    if np.isnan(entry) == True:
        print("Uh oh, we've generated an nan value")
        nanIndex=np.where(entry)
        print("The index of the nan value is: %s"%str(nanIndex[0]))


'''
Oh no, there is an nan value generated.. Why so? We can (1) fix the code

Or...

(2) Hack our way out of it. Both are good to know how to do.
'''

nanCount = np.isnan(Pr).sum()
print('# of nan values: %s'%nanCount)
Pr = Pr[np.logical_not(np.isnan(Pr))]
r=r[nanCount:] # cutoff all points that were nan in Pr dataframe (vectors must be of equal length)

'''
Normalizing the Signal
Basic function for normalizing a signal between 0-1 - Note we could bake this into PDDF
Extra credit: Write a function that normalizes integral from rmin-rmax to 1
    Hint: try numpy.trapz()
'''

def quickNormalize(sig):
    n=len(sig)
    sig=np.absolute(sig)
    maxY=np.nanmax(sig)
    minY=min(sig)
    if np.isnan(minY):
        minY=0
    Ynew=np.empty(n)
    for i in range(len(sig)):
        Ynew[i]=(sig[i]-minY)/(maxY-minY)
    return Ynew 


Pr_norm = quickNormalize(Pr) # creating a new data object (for testing).


'''
Note: Here we are adding a baseline(Y2) to improve the visualization
How can we write that?
Hint: How can you create a list of a single value repeated n times? where n = len(Pr_norm).
baseline = Y2
'''
visuals.twoPlot(X=r,Y1=Pr_norm,Y2=[0]*len(Pr_norm),savelabel='Example_Pr',plotlabel1='Pair Distance Distribution',plotlabel2='Baseline',
                  xlabel='r($\AA$)',ylabel='P(r)',linewidth=4)


'''
How can we approximate Dmax?
    (1) Iterateover the PDDF function
    (2) Save to memory the value of Dmax that gives P(rmin) & P(rmax) = 0
'''


'''Set parameters for iterations of PDDF'''

interval=10 # sets number of iterations
dmin,dmax=45,85
Dmax=np.linspace(dmin,dmax,interval)
print('Final Dmax Value in the Dmax list: %s'%Dmax[9])

'''
We need to create an empty array of n columns that we can append to throughout the iterations
'''

PDDF_list=np.empty((interval,0)).tolist()
r_range=np.empty((interval,0)).tolist()


'''
Now iterate over both the lists and append the output of PDDF() for each Dmax value
'''

for i in PDDF_list and r_range:
    for i in range(0,(interval)):
        PDDF_list[i],r_range[i] = (PDDF(shape='FoxS',Dmax=Dmax[i],I=data['I(q)'],q=data['q']))

'''
Normalize the PDDF signal between 0 and 1
'''

for i in range(len(PDDF_list)):
    PDDF_list[i]=quickNormalize(PDDF_list[i])

'''
Quick example of how to 'toss' together a plot when it is necessary to do so quickly
or it isn't worth the effort (for some reason) to make the plots publication quality

*** Matplotlib will allow us to stack as many as we want!

(1) Overlaid with the proper r_range
(2) Plotted as a function of the number of points (fixed, 99 for every profile) for better visualization
'''

# plt.plot(np.linspace(0,Dmax[0]*1.4,n),PDDF_list[0],
#          label='Dmax: %.1f'%Dmax[0] + ' ' + '$\AA$')
plt.plot(r_range[0],PDDF_list[0],
         label='Dmax: %.1f'%Dmax[0] + ' ' + '$\AA$')
plt.plot(r_range[1],PDDF_list[1],
         label='Dmax: %.1f'%Dmax[1] + ' ' + '$\AA$')
plt.plot(r_range[2],PDDF_list[2],
         label='Dmax: %.1f'%Dmax[2] + ' ' + '$\AA$')
plt.plot(r_range[3],PDDF_list[3],
         label='Dmax: %.1f'%Dmax[3] + ' ' + '$\AA$')
plt.plot(r_range[4],PDDF_list[4],
         label='Dmax: %.1f'%Dmax[4] + ' ' + '$\AA$')
plt.plot(r_range[5],PDDF_list[5],
         label='Dmax: %.1f'%Dmax[5] + ' ' + '$\AA$')
plt.plot(r_range[6],PDDF_list[6],
         label='Dmax: %.1f'%Dmax[6] + ' ' + '$\AA$')
plt.plot(r_range[7],PDDF_list[7],
         label='Dmax: %.1f'%Dmax[7] + ' ' + '$\AA$')
plt.plot(r_range[8],PDDF_list[8],
         label='Dmax: %.1f'%Dmax[8] + ' ' + '$\AA$')
plt.plot(r_range[9],PDDF_list[9],
         label='Dmax: %.1f'%Dmax[9] + ' ' + '$\AA$')
plt.ylabel('P(r)',size=14)
plt.xlabel('r ($\\AA$)',size=14)
plt.legend(loc='best')
plt.ylim(bottom=0,top=1.09)
#plt.xlim(-5,90)
#plt.ylim(bottom=0,top=4000000)
plt.savefig('VariableDmax_PDDF.png',bbox_inches='tight',dpi=300,
            format='png')
plt.show()



plt.plot(PDDF_list[0],
         label='Dmax: %s'%Dmax[0] + ' ' + '$\AA$')
plt.plot(PDDF_list[1],
         label='Dmax: %.1f'%Dmax[1] + ' ' + '$\AA$')
plt.plot(PDDF_list[2],
         label='Dmax: %.1f'%Dmax[2] + ' ' + '$\AA$')
plt.plot(PDDF_list[3],
         label='Dmax: %.1f'%Dmax[3] + ' ' + '$\AA$')
plt.plot(PDDF_list[4],
         label='Dmax: %.1f'%Dmax[4] + ' ' + '$\AA$')
plt.plot(PDDF_list[5],
         label='Dmax: %.1f'%Dmax[5] + ' ' + '$\AA$')
plt.plot(PDDF_list[6],
         label='Dmax: %.1f'%Dmax[6] + ' ' + '$\AA$')
plt.plot(PDDF_list[7],
         label='Dmax: %.1f'%Dmax[7] + ' ' + '$\AA$')
plt.plot(PDDF_list[8],
         label='Dmax: %.1f'%Dmax[8] + ' ' + '$\AA$')
plt.plot(PDDF_list[9],
         label='Dmax: %.1f'%Dmax[9] + ' ' + '$\AA$')
plt.legend(loc='best')
plt.ylim(bottom=0,top=1.09)
plt.xlim(-5,115)
plt.ylabel('P(r)',size=14)
plt.xlabel('No. of Points',size=14)
#plt.ylim(bottom=0,top=4000000)
plt.savefig('VariableDmax_PDDF2.png',bbox_inches='tight',dpi=300,
            format='png')
plt.show()


'''
Best approximation of Dmax?
'''

def dMax_Searcher(pddf):
    tosort=[]
    for i in pddf:
        n=len(i)
        tosort.append(i[n-1])
    minVal = np.nanmin(tosort)
    minIndex= tosort.index(minVal)

    return minVal, minIndex

dMaxVal, dMaxIndex = dMax_Searcher(PDDF_list)

print('The most likely Dmax value is: %.2f'%Dmax[dMaxIndex] + ' ' + 'Angstroms')

'''
How to export data? What if a collaborator wants only a data frame from your analysis, not the entire script and repository?
'''
r0,Pr0=r_range[0],PDDF_list[0]
entryCount=0

with open('example_export.csv','w',newline='') as csvfile:
    fieldnames=['Index','r0(A)','P(r)']
    thewriter=csv.DictWriter(csvfile,fieldnames=fieldnames)
    thewriter.writeheader()
    n=len(r0)
    for i in range(n):
        entryCount += 1
        thewriter.writerow({'Index':entryCount,'r0(A)':r0[i],'P(r)':Pr0[i]})

'''
Can you write a loop that exports all of the r/P(r) data?
'''