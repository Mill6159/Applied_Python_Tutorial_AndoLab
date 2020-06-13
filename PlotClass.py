#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 14:47:22 2020

@author: Rob
"""

from matplotlib import pyplot as plt

class PlotClass:
    
    def __init__(self):
        '''
        Describe the class here.
        The initialization of the class sets a few global parameters
        (1) Axes width
        (2) Line width
        (3) Set font to Helvetica. Note, if when running locally you receive a traceback error
            along the lines of "falling back to Arial" you need to enable the base Helvetica font
            provided by MacOS. If this is important to you, let me know, it's a pretty quick fix!
            Rob: rcm347@cornell.edu
        '''
        print('--------------------------------------------------------------')
        print('Plot class was called')
        print('--------------------------------------------------------------')
        
        self.axes=plt.rc('axes',linewidth=2)
        self.lines=plt.rc('lines',markeredgewidth=2)
        self.font=plt.rc('font',**{'sans-serif': ['Helvetica']})
        
        
    def basicPlot(self,X,Y,plotlabel='',savelabel='',xlabel='',ylabel=''):
        '''
        As simple as it gets
        '''
        fig=plt.figure(figsize=(10,8)) # set figure dimensions
        ax1=fig.add_subplot(1,1,1) # allows us to build more complex plots
        for tick in ax1.xaxis.get_major_ticks():
            tick.label1.set_fontsize(20) # scale for publication needs
            tick.label1.set_fontname('Helvetica')
        for tick in ax1.yaxis.get_major_ticks():
            tick.label1.set_fontsize(20) # scale for publication needs
            tick.label1.set_fontname('Helvetica')
        plt.ylabel(ylabel,size=22)
        plt.xlabel(xlabel,size=22)
        plt.plot(X,Y,'-',label=plotlabel,
                 color='#E55334') # best to use HTML color codes: https://htmlcolorcodes.com
        plt.legend(numpoints=1,fontsize=18,loc='best')
        fig.tight_layout()
        plt.savefig(savelabel+'.png',format='png',
                    bbox_inches='tight',dpi=300)
        plt.show()
        
    def semilogyPlot(self,X,Y,plotlabel='',savelabel='',xlabel='',ylabel='',linewidth=4):
        '''
        As simple as it gets
        '''
        fig=plt.figure(figsize=(10,8)) # set figure dimensions
        ax1=fig.add_subplot(1,1,1) # allows us to build more complex plots
        for tick in ax1.xaxis.get_major_ticks():
            tick.label1.set_fontsize(20) # scale for publication needs
            tick.label1.set_fontname('Helvetica')
        for tick in ax1.yaxis.get_major_ticks():
            tick.label1.set_fontsize(20) # scale for publication needs
            tick.label1.set_fontname('Helvetica')
        plt.ylabel(ylabel,size=22)
        plt.xlabel(xlabel,size=22)
        plt.semilogy(X,Y,'-',label=plotlabel,
                     linewidth=linewidth,
                     color='#E55334')
        plt.legend(numpoints=1,fontsize=18,loc='best')
        fig.tight_layout()
        plt.savefig(savelabel+'.png',format='png',
                    bbox_inches='tight',dpi=300)
        plt.show()
        
        
    def twoPlot(self,X,Y1,Y2,plotlabel1='',plotlabel2='',savelabel='',xlabel='',ylabel='',linewidth=2):
        '''
        As simple as it gets P2
        '''
        fig=plt.figure(figsize=(10,8)) # set figure dimensions
        ax1=fig.add_subplot(1,1,1) # allows us to build more complex plots
        for tick in ax1.xaxis.get_major_ticks():
            tick.label1.set_fontsize(20) # scale for publication needs
            tick.label1.set_fontname('Helvetica')
        for tick in ax1.yaxis.get_major_ticks():
            tick.label1.set_fontsize(20) # scale for publication needs
            tick.label1.set_fontname('Helvetica')
        plt.ylabel(ylabel,size=22)
        plt.xlabel(xlabel,size=22)
        plt.plot(X,Y1,'-',label=plotlabel1,
                 linewidth=linewidth,
                 color='#E55334')
        plt.plot(X,Y2,'o',label=plotlabel2,
                 color='#1283BC')
        plt.legend(numpoints=1,fontsize=18,loc='best')
        fig.tight_layout()
        plt.savefig(savelabel+'.png',format='png',
                    bbox_inches='tight',dpi=300)
        plt.show()