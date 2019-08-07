# ----------- Data Visualisation ----------- #
# Written by Benjamin J. Griffiths
# Created on Friday 26th October, 2018
# ------------------------------------------ #

# %% -- import modules ------------------------------------------------------ #
import matplotlib as mpl
import numpy as np
import pandas
import ptitprince as pt
import seaborn as sns
from scipy import stats
from matplotlib import pyplot
import matplotlib.font_manager as font_manager

# %% -- define functions ---------------------------------------------------- #
# plot raincloud
def eeg_raincomparison(data,colour,axes,labels,ylim,condition=None):
       
    # set condition
    if type(condition) != np.ndarray:
        condition = np.append(np.zeros(len(data)))
    
    # plot violin
    axes=pt.half_violinplot(data = data,bw = "scott",inner = None, scale = "count",
                          width = 0.5, linewidth = 1, cut = 1, palette = colour,
                          orient='flip',ax = ax, edgecolor = [0,0,0],flip=True)
    
    # plot mean and confidence intervals
    axes=sns.boxplot(data = data, palette = colour, width = 0.1, ax = axes, linewidth = 1, fliersize = 1)
    
    # plot subjects
    val = data.values
    for i in range(len(data)):
        if condition[i]==0:
            m='o'
            fc=colour
        else:
            m='^'
            fc=[[1,1,1],[1,1,1]]
            
        pyplot.scatter([0.18,0.18],[val[i,0],val[i,0]],marker=m,c=fc[0],edgecolors=colour[0],linewidths=0.5)
        pyplot.scatter([0.82,0.82],[val[i,1],val[i,1]],marker=m,c=fc[1],edgecolors=colour[1],linewidths=0.5)
        pyplot.plot([0.18,0.82],val[i,:],linewidth=.5,c=[0.7,0.7,0.7],zorder=-2)
    
    # add title
    if ('title' in labels):
        pyplot.title(labels['title'],fontdict = {'fontsize':6,'fontname':'Arial','fontweight':'bold'})
          
    # aesthetics
    axes.set_ylabel(labels['ylabel'],fontname='Arial',fontsize=6,labelpad=3,fontweight='light')   # add Y axis label
    axes.set_ylim(ylim)                  # set Y axis range to 0 -> 1
    axes.set_xlim(-0.65,-0.5+len(data.columns))                  # set Y axis range to 0 -> 1
    axes.tick_params(axis='x',          # change X tick parameters
                   which='both',          # affect both major and minor ticks
                   bottom=False,          # turn off bottom ticks
                   labelbottom=True,  # keep bottom labels
                   pad=2.5,
                   width=1)   
    axes.tick_params(axis='y',          # change X tick parameters
                       pad=2,
                       width=1,
                       length=2.5)
    axes.set_yticks(labels['yticks'])
    axes.set_xticklabels(labels['xticklabel'],fontname='Arial',fontweight='light',fontsize=6)
    axes.set_yticklabels(labels['yticklabel'],fontname='Arial',fontweight='light',fontsize=6)
    
    # change axes
    axes.spines['top'].set_visible(False)
    axes.spines['right'].set_visible(False)
    axes.spines['bottom'].set_visible(False)
    axes.spines['left'].set_linewidth(1)
    
    
def eeg_anovaboxplot(data,colour,condition,ax,labels,ylim,p):
        
    # fix colour for box plot
    boxcol = [colour[0],colour[1],colour[0],colour[2],colour[3]]
    
    # create boxplot data
    y = np.reshape(data.values,-1)
    x = np.reshape(np.tile([0, 1, 3, 4],(1,len(data))),-1)
    boxdata = pandas.DataFrame(dict(x=x,y=y))
    
    # set condition
    if type(condition) != np.ndarray:
        condition = np.append(np.zeros(len(data)))
    
    # plot subjects
    val = data.values
    for i in range(len(data)):
        if condition[i]==0:
            m='o'
            fc=colour
        else:
            m='^'
            fc=[[1,1,1],[1,1,1]]
            
        # plot encoding scatter
        pyplot.scatter([-1.2,-1.2],[val[i,0],val[i,0]],marker=m,c=fc[0],edgecolors=colour[0],linewidths=0.5)
        pyplot.scatter([-0.8,-0.8],[val[i,1],val[i,1]],marker=m,c=fc[1],edgecolors=colour[1],linewidths=0.5)
        
        # plot retrieval scatter
        pyplot.scatter([4.8,4.8],[val[i,2],val[i,2]],marker=m,c=fc[0],edgecolors=colour[0],linewidths=0.5)
        pyplot.scatter([5.2,5.2],[val[i,3],val[i,3]],marker=m,c=fc[1],edgecolors=colour[1],linewidths=0.5)
        
    # plot interaction lines
    med = np.median(val,axis=0)
    pyplot.plot([0.2,2.8],[med[0],med[2]],color=colour[0],zorder=-1)
    pyplot.plot([1.2,3.8],[med[1],med[3]],color=colour[1],zorder=-1)
    
    # plot box
    axes=sns.boxplot(x='x',y='y',data = boxdata, palette = boxcol, width = 0.4, 
                     ax = ax, linewidth = 1, fliersize = 1, order = np.arange(5))
    
    # plot significance
    if (p <= 0.05) & (p > 0.01):
        px = np.linspace(0.5,3.5,19)
        px = px[9]
        pyplot.plot([0.5,3.5],[ylim[1]*0.9,ylim[1]*0.9],c=[0.2,0.2,0.2])
        pyplot.scatter([px,px],[ylim[1]*0.95,ylim[1]*0.95],s=1,c=[0.2,0.2,0.2],marker='*',edgecolors=None)
    elif (p <= 0.05) & (p > 0.001):
        px = np.linspace(0.5,3.5,20)
        px = px[9:11]
        pyplot.plot([0.5,3.5],[ylim[1]*0.9,ylim[1]*0.9],c=[0.2,0.2,0.2])
        pyplot.scatter(px,[ylim[1]*0.95,ylim[1]*0.95],s=1,c=[0.2,0.2,0.2],marker='*',edgecolors=None)
    elif (p <= 0.001):
        px = np.linspace(0.5,3.5,19)
        px = px[8:11]
        pyplot.plot([0.5,3.5],[ylim[1]*0.9,ylim[1]*0.9],c=[0.2,0.2,0.2])
        pyplot.scatter(px,[ylim[1]*0.95,ylim[1]*0.95,y[1]*0.9],s=1,c=[0,0,0],marker='*',edgecolors=None)   
    
    # add title
    if ('title' in labels):
        pyplot.title(labels['title'],fontdict = {'fontsize':6,'fontname':'Arial','fontweight':'bold'})
          
    # aesthetics
    axes.set_ylabel(labels['ylabel'],fontname='Arial',fontsize=6,labelpad=0,fontweight='light')
    axes.set_xlim([-2,6])  
    axes.set_ylim(ylim)  
    axes.tick_params(axis='x',          # change X tick parameters
                   which='both',          # affect both major and minor ticks
                   bottom=False,          # turn off bottom ticks
                   labelbottom=True,  # keep bottom labels
                   pad=1,
                   width=1)   
    axes.tick_params(axis='y',          # change X tick parameters
                       pad=2,
                       width=1,
                       length=2.5)
    axes.set_yticks(labels['yticks'])
    axes.set_xlabel('')
    axes.set_xticks([0.5,3.5])
    axes.set_xticklabels(labels['xticklabel'],fontname='Arial',fontweight='light',fontsize=6)
    axes.set_yticklabels(labels['yticklabel'],fontname='Arial',fontweight='light',fontsize=6)
    
    # set legend
    font = font_manager.FontProperties(family='Arial',weight='light',size=6)
    ax.legend(labels['legend'],frameon=False,fontsize=6, 
              loc='lower right', prop=font)
    
    # change axes
    axes.spines['top'].set_visible(False)
    axes.spines['right'].set_visible(False)
    axes.spines['bottom'].set_visible(False)
    axes.spines['left'].set_linewidth(1)
               
    
def eeg_timeseriesplot(data,variables,limits,labels,colour,ax,sig=None):
 
   # calculate y scale limits
   if ('ylim' in limits) & (type(limits['ylim'])!=list):  
      
      # calculate mean/sem
      if ('hue' in variables):
         for i in np.unique(data['condition']):
            dat = data['signal'][data['condition']==i]
            ns  = np.shape(np.unique(data['subj']))
            dat = np.reshape(dat.values,(ns[0],int(np.shape(dat)[0]/ns[0])))
            ymean = np.mean(dat,axis=0)
            ysem = np.std(dat,axis=0) / np.sqrt(np.size(dat,axis=0))
            ypos = max(ymean + ysem)
            yneg = min(ymean - ysem)
            yabs = max(np.abs([ypos,yneg]))
            
      else:
         dat = data['signal']
         ns  = np.shape(np.unique(data['subj']))
         dat = np.reshape(dat.values,(ns[0],int(np.shape(dat)[0]/ns[0])))
         ymean = np.mean(dat,axis=0)
         ysem = np.std(dat,axis=0) / np.sqrt(np.size(dat,axis=0))
         ypos = max(ymean + ysem)
         yneg = min(ymean - ysem)
         yabs = max(np.abs([ypos,yneg]))

      # if minmax request
      if limits['ylim'] == 'minmax':
         limits['ylim'] = [yneg-(yneg*.1),ypos+(ypos*.1)]
         
      # else if maxabs requested   
      elif limits['ylim'] == 'maxabs':
         limits['ylim'] = [yabs*-1.1,yabs*1.1]
         
      # else if zero to max requested   
      elif limits['ylim'] == 'zeromax':
         limits['ylim'] = [0,ypos*1.1]
         
   # plot significance
   if type(sig)==list:
       
       # check type of first list element
       if (type(sig[0])==int) | (type(sig[0])==float):
           sigout = [sig]
       else:
           sigout = sig
       
       # cycle through each inputted significance value
       for j in range(len(sigout)):
        
           # get key variables
           x = sigout[j][0:2]
           xdiff = x[1] - x[0]
           xbar  = [x[0]+(xdiff*0.1),x[1]-(xdiff*0.1)]
           p = sigout[j][2]
           y = limits['ylim']
           
           # plot
           pyplot.fill_between(x,[y[0],y[0]],[y[1],y[1]],
                               edgecolor=[0.9,0.9,0.9],facecolor=[0.9,0.9,0.9])
           if (p <= 0.05) & (p > 0.01):
               px = np.linspace(xbar[0],xbar[1],3)
               px = px[1]
               pyplot.plot(xbar,[y[1]*0.8,y[1]*0.8],c=[0.2,0.2,0.2])
               pyplot.scatter([px,px],[y[1]*0.9,y[1]*0.9],s=1,c=[0.2,0.2,0.2],marker='*',edgecolors=None)
           elif (p <= 0.05) & (p > 0.001):
               px = np.linspace(xbar[0],xbar[1],4)
               px = px[1:3]
               pyplot.plot(xbar,[y[1]*0.8,y[1]*0.8],c=[0.2,0.2,0.2])
               pyplot.scatter(px,[y[1]*0.9,y[1]*0.9],s=1,c=[0.2,0.2,0.2],marker='*',edgecolors=None)
           elif (p <= 0.001):
               px = np.linspace(xbar[0],xbar[1],5)
               px = px[1:4]
               pyplot.plot(xbar,[y[1]*0.8,y[1]*0.8],c=[0.2,0.2,0.2])
               pyplot.scatter(px,[y[1]*0.9,y[1]*0.9,y[1]*0.9],s=1,c=[0,0,0],marker='*',edgecolors=None)
   
   # plot multicondition timeseries
   if ('hue' in variables):
        sns.lineplot(x=variables['x'],
                     y=variables['y'],
                     data=data,
                     hue=variables['hue'],
                     hue_order=[0,1],
                     palette=colour,
                     ax = ax,
                     ci='sem',
                     linewidth=1)
      
      # add leened
        if ('legend' in labels):
            font = font_manager.FontProperties(family='Arial',weight='light',size=7)
            ax.legend(prop=font)
            ax.legend(labels['legend'],frameon=False,fontsize=6,bbox_to_anchor=(0., 1.02, 1., .102), 
                      loc=3, ncol=2, mode="expand", borderaxespad=0.,prop=font)
        else:
            ax.get_legend().remove()   
      
   # plot single condition timeseries   
   else:
       # add hue condition
       data = data.assign(null=np.zeros(len(data)))
       
       # plot
       sns.lineplot(x=variables['x'],
                    y=variables['y'],
                    data=data,
                    color=colour,
                    ax = ax,
                    ci='sem',
                    linewidth=1)
               
      
   # add title
   if ('title' in labels):
       pyplot.title(labels['title'],fontdict = {'fontsize':6,'fontname':'Arial','fontweight':'bold'})
      
   # add horizontal line
   if limits['xline']:
      ax.axvline(x=0, linewidth = 1, color = [0,0,0], linestyle='--')
   if limits['yline']:
      ax.axhline(y=0, linewidth = 1, color = [0,0,0], linestyle='-')
   
   # set x scale limits
   if ('xlim' in limits):   
      
      # set limits
      ax.set_xlim(limits['xlim'])
          
   # set scale
   ax.set_ylim(limits['ylim'])
   
   # set scaling factor
   if ('xscale' in limits): 
      ax.set_yscale(value=limits['xscale'])
   if ('yscale' in limits): 
      ax.set_yscale(value=limits['yscale'])
   
   # set x ticks
   if ('xticks' in limits):
      
      # if single value specified, get linspaced ticks
      if (type(limits['xticks']) == int) & ('xlim' in limits):
         limits['xticks'] = np.linspace(limits['xlim'][0],limits['xlim'][1],limits['xticks'])
      if type(limits['xticks']) == int:
         limits['xticks'] = np.linspace(np.min(data[variables['x']]),np.max(data[variables['x']]),limits['xticks'])
      
      # set ticks
      ax.set_xticks(limits['xticks'])
      
      # set ticks to precision
      if ('xprecision' in limits):
         limits['xticks'] = np.round(limits['xticks'],limits['xprecision'])      
      
      # set tick labels
      ax.set_xticklabels(limits['xticks'],fontname='Arial',fontweight='light',fontsize=6)
         
   # set x ticks
   if ('yticks' in limits):
      
      # if single value specified, get linspaced ticks
      if (type(limits['yticks']) == int) & ('ylim' in limits):
         limits['yticks'] = np.linspace(limits['ylim'][0],limits['ylim'][1],limits['yticks'])
      elif (type(limits['yticks']) == int):      
         limits['yticks'] = np.linspace(np.min(data[variables['y']]),np.max(data[variables['y']]),limits['yticks'])
        
      # set ticks
      ax.set_yticks(limits['yticks'])
         
      # set ticks to precision
      if ('yprecision' in limits):
         limits['yticks'] = np.round(limits['yticks'],limits['yprecision'])   
      
      # set tick labels
      ax.set_yticklabels(limits['yticks'],fontname='Arial',fontweight='light',fontsize=6)
         
   
   # aesthetics
   ax.set_ylabel(labels['ylabel'],fontname='Arial',fontsize=6,labelpad=2,fontweight='light')   # add Y axis label
   ax.set_xlabel(labels['xlabel'],fontname='Arial',fontsize=6,labelpad=2,fontweight='light')   # add Y axis label
   ax.tick_params(axis='both',          # change X tick parameters
                      pad=3,
                      length=2.5)
   
   # change axes
   ax.spines['top'].set_visible(False)
   ax.spines['right'].set_visible(False)


# %% -- define key parameters ----------------------------------------------- #
# get current directory
wdir = 'E:/bjg335/projects/hippcortex-interactions/'

# set context
sns.set_context("paper")

# set plot defaults
mpl.rcParams['xtick.major.width'] = 1
mpl.rcParams['ytick.major.width'] = 1
mpl.rcParams['xtick.color'] = [0,0,0]
mpl.rcParams['ytick.color'] = [0,0,0]
mpl.rcParams['lines.linewidth'] = 1
mpl.rcParams['axes.linewidth'] = 1
mpl.rcParams['lines.markersize'] = 2

# get colour dict
RdBu = sns.color_palette("RdBu_r", 10)
Purp = sns.color_palette("Purples", 10)
colpal = {'Red':RdBu[8],
          'Blue':RdBu[1],
          'Gray':[0.7,0.7,0.7],
          'Purple':Purp[7]}
del(RdBu,Purp)


# %% -- Figure 1a -- BRAINS!!!!!!! ------------------------------------------ #
cod = pandas.read_csv(wdir + "data/fig1_coords.csv", delimiter=",",header=None)
cod = cod.values
ddim = np.shape(cod)

# create figure
f,ax = pyplot.subplots(1,1)
f.set_figheight(2/3) # convert size to inches
f.set_figwidth(2.4/3)
f.set_dpi(1000)

# load data data
mri = pandas.read_csv(wdir + "data/fig1_mrixy.csv", delimiter=",",header=None)
mri = mri.values
mri[mri==0] = 130
ax.imshow(mri,cmap='gray')
for i in range(ddim[0]):
    if cod[i,3] == 1:
        ax.scatter([cod[i,0],cod[i,0]],[cod[i,1],cod[i,1]],s=1,marker='o',edgecolors=colpal['Red'],c=colpal['Red'],linewidths=.05)
    else:        
        ax.scatter([cod[i,0],cod[i,0]],[cod[i,1],cod[i,1]],s=1,marker='o',edgecolors=colpal['Blue'],c=colpal['Blue'],linewidths=.05)
ax.set_yticks([])
ax.set_xticks([])
ax.spines['top'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.spines['right'].set_visible(False)
pyplot.savefig(wdir + "figures/fig1a.tif",bbox_inches='tight',transparent=True,dpi='figure')

# create figure
f,ax = pyplot.subplots(1,1)
f.set_figheight(2/3) # convert size to inches
f.set_figwidth(2.4/3)
f.set_dpi(1000)

# load data data
mri = pandas.read_csv(wdir + "data/fig1_mrixz.csv", delimiter=",",header=None)
mri = mri.values
mri[mri==0] = 130
ax.imshow(mri,cmap='gray')
for i in range(ddim[0]):
    if cod[i,3] == 1:
        ax.scatter([cod[i,0],cod[i,0]],[cod[i,2],cod[i,2]],s=1,marker='o',edgecolors=colpal['Red'],c=colpal['Red'],linewidths=.05)
    else:
        ax.scatter([cod[i,0],cod[i,0]],[cod[i,2],cod[i,2]],s=1,marker='o',edgecolors=colpal['Blue'],c=colpal['Blue'],linewidths=.05)
ax.set_yticks([])
ax.set_xticks([])
ax.spines['top'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.spines['right'].set_visible(False)
pyplot.savefig(wdir + "figures/fig1b.tif",bbox_inches='tight',transparent=True,dpi='figure')

# create figure
f,ax = pyplot.subplots(1,1)
f.set_figheight(2/3) # convert size to inches
f.set_figwidth(2.4/3)
f.set_dpi(1000)

# load data data
mri = pandas.read_csv(wdir + "data/fig1_mriyz.csv", delimiter=",",header=None)
mri = mri.values
mri[mri==0] = 130
ax.imshow(mri,cmap='gray')
for i in range(ddim[0]):
    if cod[i,3] == 1:
        ax.scatter([cod[i,1],cod[i,1]],[cod[i,2],cod[i,2]],s=1,marker='o',edgecolors=colpal['Red'],c=colpal['Red'],linewidths=.05)
    else:
        ax.scatter([cod[i,1],cod[i,1]],[cod[i,2],cod[i,2]],s=1,marker='o',edgecolors=colpal['Blue'],c=colpal['Blue'],linewidths=.05)
ax.set_yticks([])
ax.set_xticks([])
ax.spines['top'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.spines['right'].set_visible(False)
pyplot.savefig(wdir + "figures/fig1c.tif",bbox_inches='tight',transparent=True,dpi='figure')

# %% -- Figure 1b -- BRAINS!!!!!!! ------------------------------------------ #
data = pandas.read_csv(wdir + "data/fig1_count.csv", delimiter=",",header=None,names=['count'])
data = data.assign(subj=np.tile(np.arange(0,12),2),roi=np.concatenate((np.zeros(12),np.ones(12))))

# create figure
f,ax = pyplot.subplots(1,1)
f.set_figheight(1.6/2.54) # convert size to inches
f.set_figwidth(4.4/2.54)
f.set_dpi(1000)

# plot bar
sns.barplot(x='subj',y='count',hue='roi',data=data,ax=ax,
            palette=[colpal['Blue'],colpal['Red']])

ax.set_ylabel('Elec. Count',fontname='Arial',fontsize=6,labelpad=2,fontweight='light')   # add Y axis label
ax.set_xlabel('Participant',fontname='Arial',fontsize=6,labelpad=2,fontweight='light')   # add Y axis label
ax.tick_params(axis='both',pad=3,length=2.5)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
        
# aesthetics
ax.set_xlim([-1,12])  
ax.set_ylim([0,8])  
ax.tick_params(axis='x',          # change X tick parameters
               which='both',          # affect both major and minor ticks
               labelbottom=True,  # keep bottom labels
               pad=1,
               width=1)   
ax.tick_params(axis='y',          # change X tick parameters
                   pad=2,
                   width=1,
                   length=2.5)
ax.set_yticks(np.arange(0,9,2))
ax.set_xticks(range(12))
ax.set_xticklabels(np.arange(12)+1,fontname='Arial',fontweight='light',fontsize=6)
ax.set_yticklabels(['0','','4','','8'],fontname='Arial',fontweight='light',fontsize=6)

# set legend
font = font_manager.FontProperties(family='Arial',weight='light',size=6)
ax.legend(['ATL','Hipp.'],frameon=False,fontsize=6,bbox_to_anchor=(0., 1.02, 1., .102), 
                      loc=3, ncol=2, mode="expand", borderaxespad=0.,prop=font)
pyplot.savefig(wdir + "figures/fig1d.tif",bbox_inches='tight',transparent=True,dpi='figure')

# %% -- Figure 2a -- Resonating Freqs. -------------------------------------- #
# -- Prep Data -- #
# load data data
dataATL = pandas.read_csv(wdir + "data/fig2_resATL.csv", delimiter=",",header=None)
dataHipp = pandas.read_csv(wdir + "data/fig2_resHipp.csv", delimiter=",",header=None)
ddim = np.shape(dataATL)

# normalise signal
sigATL = stats.zscore(dataATL.values,axis=1)
sigHipp = stats.zscore(dataHipp.values,axis=1)

# get raw data
signal = np.hstack((np.reshape(sigATL,np.prod(ddim)),np.reshape(sigHipp,np.prod(ddim))))
subj = np.tile(np.repeat(np.arange(0,ddim[0]),[ddim[1]]),2)
freq = np.tile(np.linspace(1.5,100,ddim[1]),ddim[0]*2)
condition = np.hstack(np.array([np.zeros(np.prod(ddim)),np.ones(np.prod(ddim))]))

# reformat data as frame
data = pandas.DataFrame(data=np.rot90([signal,subj,freq,condition]),columns=['signal','subj','freq','condition']) # reshape signal
del(signal,subj,freq,condition,ddim,sigATL,sigHipp,dataATL,dataHipp)

# -- Plot -- #
# create figure
f,ax = pyplot.subplots(1,1)
f.set_figheight(3.6/3) # convert size to inches
f.set_figwidth(6.1/3)
f.set_dpi(1000)

# define labels and variables
labels = {'legend':['ATL','Hippocampus'],'ylabel':'1/f Corrected Power (z)','xlabel':'Frequency (Hz.)'}
variables = {'x':'freq','y':'signal','hue':'condition'}
limits = {'xline':False,'yline':False,'xlim':[1.5,64],'ylim':[-1,4],
          'xticks':[4,8,16,32,64],'yticks':[-2,0,2,4],'xprecision':1,'yprecision':1}

# plot
eeg_timeseriesplot(data,variables,limits,labels,[colpal['Blue'],colpal['Red']],ax)
pyplot.savefig(wdir + "figures/fig2a.tif",bbox_inches='tight',transparent=True,dpi='figure')
del(data,variables,limits,labels)


# %% -- Figure 2b -- Enc > Ret Spectrum ------------------------------------- #
# -- Prep Data -- #
# load raincloud data
data = pandas.read_csv(wdir + "data/fig2_encRetGammaSpecHits.csv", delimiter=",",header=None)
ddim = np.shape(data)

# normalise signal
signal = np.reshape(data.values,np.prod(ddim))

# get raw data
subj = np.repeat(np.arange(0,ddim[0]),[ddim[1]])
freq = np.tile(np.linspace(30,100,ddim[1]),ddim[0])

# reformat data as frame
data = pandas.DataFrame(data=np.rot90([signal,subj,freq]),columns=['signal','subj','freq']) # reshape signal
del(ddim,subj,freq,signal)

# -- Plot -- #
# create figure
f,ax = pyplot.subplots(1,1)
f.set_figheight(3.6/3) # convert size to inches
f.set_figwidth(6.2/3)
f.set_dpi(1000)

# define labels and variables
labels = {'ylabel':'1/f Corrected Power\n(Enc. > Ret.; a.u.)','xlabel':'Frequency (Hz.)'}
variables = {'x':'freq','y':'signal'}
limits = {'xline':False,'yline':True,'xlim':[30,100],'ylim':[-.06,.06],
          'xticks':np.arange(40,100,10),'yticks':[-.06,0,.06],'xprecision':3,'yprecision':3}

# plot
eeg_timeseriesplot(data,variables,limits,labels,colpal['Red'],ax,[[40,50,0.02],[60,70,0.001],[70,80,0.02]])
pyplot.savefig(wdir + "figures/fig2b.tif",bbox_inches='tight',transparent=True,dpi='figure')
del(data,labels,variables,limits)


# %% -- Figure 2c -- Raw Traces --------------------------------------------- #
# -- Prep Data -- #
# load data
dataSlow = pandas.read_csv(wdir + "data/fig2_45HzRawTrace.csv", delimiter=",",header=None)
dataFast = pandas.read_csv(wdir + "data/fig2_70HzRawTrace.csv", delimiter=",",header=None)
ddim = np.size(dataSlow)

# get raw data
timeSlow = np.linspace(21.7,21.95,ddim)
timeFast = np.linspace(32.0,32.25,ddim)

# reformat data as frame
dataSlow = pandas.DataFrame(data=np.rot90([np.ndarray.flatten(dataSlow.values),timeSlow]),columns=['signal','time']) # reshape signal
dataFast = pandas.DataFrame(data=np.rot90([np.ndarray.flatten(dataFast.values),timeFast]),columns=['signal','time']) # reshape signal
del(timeSlow,timeFast,ddim)

# define labels and variables
labels = {'ylabel':'Amp. (uV)','xlabel':'','title':'Raw Slow Gamma Trace'}
variables = {'x':'time','y':'signal'}
limits = {'xline':False,'yline':False,'xlim':[21.7,21.95],'ylim':[-100,100],
          'xticks':[],'yticks':[-100,0,100],'xprecision':2,'yprecision':2}

# plot slow gamma
f,ax = pyplot.subplots(1,1)
f.set_figheight(1.8/5) # convert size to inches
f.set_figwidth(6/5)
f.set_dpi(1000)
eeg_timeseriesplot(dataSlow,variables,limits,labels,colpal['Purple'],ax,[21.8,21.85,1])
pyplot.savefig(wdir + "figures/fig2ci.tif",bbox_inches='tight',transparent=True,dpi='figure')

# plot fast gamma
f,ax = pyplot.subplots(1,1)
f.set_figheight(1.8/5) # convert size to inches
f.set_figwidth(6/5)
f.set_dpi(1000)
labels = {'ylabel':'Amp. (uV)','xlabel':'50ms','title':'Raw Fast Gamma Trace'}
limits = {'xline':False,'yline':False,'xlim':[32.0,32.25],'ylim':[-100,100],
          'xticks':[],'yticks':[-100,0,100],'xprecision':2,'yprecision':2}
eeg_timeseriesplot(dataFast,variables,limits,labels,colpal['Red'],ax,[32.1,32.15,1])
pyplot.savefig(wdir + "figures/fig2cii.tif",bbox_inches='tight',transparent=True,dpi='figure')
del(dataFast,dataSlow,variables,limits,labels)


# %% -- Figure 2d -- Peak Locked Averages  ---------------------------------- #
# -- Prep Data -- #
# load data
dataSlow = pandas.read_csv(wdir + "data/fig2_45HzPeakAvg.csv", delimiter=",",header=None)
dataFast = pandas.read_csv(wdir + "data/fig2_70HzPeakAvg.csv", delimiter=",",header=None)
ddim = np.shape(dataSlow)

# normalise signal
signalSlow = np.reshape(dataSlow.values,np.prod(ddim))
signalFast = np.reshape(dataFast.values,np.prod(ddim))

# get raw data
subj = np.repeat(np.arange(0,ddim[0]),[ddim[1]])
tlim = np.floor(ddim[1]/2)/500
time = np.tile(np.linspace(-tlim,tlim,ddim[1]),[ddim[0]])

# reformat data as frame
dataSlow = pandas.DataFrame(data=np.rot90([signalSlow,subj,time]),columns=['signal','subj','time']) # reshape signal
dataFast = pandas.DataFrame(data=np.rot90([signalFast,subj,time]),columns=['signal','subj','time']) # reshape signal
del(ddim,subj,time,tlim,signalSlow,signalFast)

# define labels and variables
labels = {'ylabel':'Amp. (uV)','xlabel':'','title':'Slow Gamma'}
variables = {'x':'time','y':'signal'}
limits = {'xline':False,'yline':False,'xlim':[-0.1,0.1],'ylim':[-10,10],
          'xticks':[],'yticks':[-10,0,10],'xprecision':2,'yprecision':2}

# plot slow gamma
f,ax = pyplot.subplots(1,1)
f.set_figheight(1.9/5) # convert size to inches
f.set_figwidth(3.3/5)
f.set_dpi(1000)
eeg_timeseriesplot(dataSlow,variables,limits,labels,colpal['Purple'],ax)
pyplot.savefig(wdir + "figures/fig2di.tif",bbox_inches='tight',transparent=True,dpi='figure')

# plot fast gamma
f,ax = pyplot.subplots(1,1)
f.set_figheight(1.9/5) # convert size to inches
f.set_figwidth(3.3/5)
f.set_dpi(1000)
labels = {'ylabel':'Amp. (uV)','xlabel':'Time (s)','title':'Fast Gamma'}
limits = {'xline':False,'yline':False,'xlim':[-0.1,0.1],'ylim':[-10,10],
          'xticks':[-0.1,0,0.1],'yticks':[-10,0,10],'xprecision':2,'yprecision':2}
eeg_timeseriesplot(dataFast,variables,limits,labels,colpal['Red'],ax)
pyplot.savefig(wdir + "figures/fig2dii.tif",bbox_inches='tight',transparent=True,dpi='figure')


# %% -- Figure 2e. -- Gamma Raincloud Contrast ------------------------------ #
# load raincloud data
rain_45Hz = pandas.read_csv(wdir + "data/fig2_45HzGammaRain.csv", delimiter=",",header=None,names=['encoding','retrieval'])
rain_70Hz = pandas.read_csv(wdir + "data/fig2_70HzGammaRain.csv", delimiter=",",header=None,names=['encoding','retrieval'])
condition = np.append(np.zeros(7),np.ones(5))

# create figure
f,ax = pyplot.subplots(1,1)
f.set_figheight(4.5/2.54) # 4inches 
f.set_figwidth(4.0/2.54) # 12inches
f.set_dpi(1000)

# define labels
labels = {'title':'Slow Gamma (40-50Hz)',
          'ylabel':'Gamma Power (a.u.)',
          'xticklabel':['Encoding','Retrieval'],
          'yticks':[1.15,1.25,1.35,1.45],
          'yticklabel':['1.15','1.25','1.35','1.45']}

# plot and save
eeg_raincomparison(rain_45Hz,[colpal['Gray'],colpal['Purple']],ax,labels,[1.15,1.45],condition)
pyplot.savefig(wdir + "figures/fig2ei.tif",bbox_inches='tight',transparent=True,dpi='figure')

# create figure
f,ax = pyplot.subplots(1,1)
f.set_figheight(4.5/2.54) # 4inches 
f.set_figwidth(4.0/2.54) # 12inches
f.set_dpi(1000)

# define labels
labels = {'title':'Fast Gamma (60-80Hz)',
          'ylabel':'Gamma Power (a.u.)',
          'xticklabel':['Encoding','Retrieval'],
          'yticks':[1.2,1.3,1.4,1.5,1.6,1.7,1.8],
          'yticklabel':['1.25','1.35','1.45','1.55','1.65','1.75','1.85']}

# plot and save
eeg_raincomparison(rain_70Hz,[colpal['Red'],colpal['Gray']],ax,labels,[1.25,1.75],condition)
pyplot.savefig(wdir + "figures/fig2eii.tif",bbox_inches='tight',transparent=True,dpi='figure')


# %% -- Figure 2f. -- Gamma TimeSeries -------------------------------------- #
# -- Prep Data -- #
# load raincloud data
dataEnc45Hz = pandas.read_csv(wdir + "data/fig2_hippEnc45HzTime.csv", delimiter=",",header=None)
dataRet45Hz = pandas.read_csv(wdir + "data/fig2_hippRet45HzTime.csv", delimiter=",",header=None)
dataEnc70Hz = pandas.read_csv(wdir + "data/fig2_hippEnc70HzTime.csv", delimiter=",",header=None)
dataRet70Hz = pandas.read_csv(wdir + "data/fig2_hippRet70HzTime.csv", delimiter=",",header=None)
ddim = np.shape(dataEnc45Hz)

# normalise signal
signalEnc45Hz = np.reshape(dataEnc45Hz.values,np.prod(ddim))
signalRet45Hz = np.reshape(dataRet45Hz.values,np.prod(ddim))
signalEnc70Hz = np.reshape(dataEnc70Hz.values,np.prod(ddim))
signalRet70Hz = np.reshape(dataRet70Hz.values,np.prod(ddim))
signalEnc = np.concatenate((signalEnc45Hz,signalEnc70Hz))
signalRet = np.concatenate((signalRet45Hz,signalRet70Hz))

# get raw data
subj = np.tile(np.repeat(np.arange(0,ddim[0]),[ddim[1]]),2)
time = np.tile(np.linspace(0,1.5,ddim[1]),ddim[0]*2)
condition = np.concatenate((np.zeros(np.size(signalEnc45Hz)),np.ones(np.size(signalEnc45Hz))))

# reformat data as frame
dataEnc = pandas.DataFrame(data=np.rot90([signalEnc,subj,time,condition]),columns=['signal','subj','time','condition']) # reshape signal
dataRet = pandas.DataFrame(data=np.rot90([signalRet,subj,time,condition]),columns=['signal','subj','time','condition']) # reshape signal
del(signalEnc45Hz,signalRet45Hz,subj,time,ddim)

# -- Plot -- #
# define labels and variables
labels = {'ylabel':'Gamma Power\n(Rem. > Forg.; z)','xlabel':'','title':'Encoding'}
variables = {'x':'time','y':'signal','hue':'condition'}
limits = {'xline':False,'yline':True,'xlim':[0,1.5],'ylim':[-.6,.6],
          'xticks':[0,0.5,1,1.5],'yticks':[-.6,0,.6],'xprecision':1,'yprecision':1}

# plot encoding
f,ax = pyplot.subplots(1,1)
f.set_figheight(1.4/2.54) # convert size to inches
f.set_figwidth(2.8/2.54)
f.set_dpi(1000)
eeg_timeseriesplot(dataEnc,variables,limits,labels,[colpal['Purple'],colpal['Red']],ax)
pyplot.savefig(wdir + "figures/fig2fi.tif",bbox_inches='tight',transparent=True,dpi='figure')

# plot retrieval
f,ax = pyplot.subplots(1,1)
f.set_figheight(1.4/2.54) # convert size to inches
f.set_figwidth(2.8/2.54)
f.set_dpi(1000)
labels = {'ylabel':'Gamma Power\n(Rem. > Forg.; z)','xlabel':'Time (s)','title':'Retrieval'}
eeg_timeseriesplot(dataRet,variables,limits,labels,[colpal['Purple'],colpal['Red']],ax)
pyplot.savefig(wdir + "figures/fig2fii.tif",bbox_inches='tight',transparent=True,dpi='figure')
del(labels,limits,variables,dataEnc,dataRet)


# %% -- Figure 2h. -- SME Gamma ANOVA --------------------------------------- %
# load raincloud data
data = pandas.read_csv(wdir + "data/fig2_hippGammaANOVA.csv", delimiter=",",header=None)

condition = np.append(np.zeros(7),np.ones(5))
colour = [colpal['Purple'],colpal['Red'],colpal['Purple'],colpal['Red']]

# create figure
f,ax = pyplot.subplots(1,1)
f.set_figheight(5/2.54) # 4inches 
f.set_figwidth(4.4/2.54) # 12inches
f.set_dpi(1000)

# define labels
labels = {'title':'',
          'ylabel':'Gamma Power (Hits > Misses; z)',
          'xticklabel':['Encoding','Retrieval'],
          'yticks':[-0.7,-0.35,0,0.35,0.7],
          'yticklabel':['-0.7','','0','','0.7'],
          'legend':['Slow Gamma','Fast Gamma']}

# plot and save
eeg_anovaboxplot(data,colour,condition,ax,labels,[-0.7,0.7],0.003)
pyplot.savefig(wdir + "figures/fig2g.tif",bbox_inches='tight',transparent=True,dpi='figure')


# %% -- Figure 3a. -- Neocortical TimeSeries -------------------------------- #
# -- Prep Data -- #
# load raincloud data
dataEnc = pandas.read_csv(wdir + "data/fig2_atlEncTime.csv", delimiter=",",header=None)
dataRet = pandas.read_csv(wdir + "data/fig2_atlRetTime.csv", delimiter=",",header=None)
ddim = np.shape(dataEnc)

# normalise signal
signalEnc = np.reshape(dataEnc.values,np.prod(ddim))
signalRet = np.reshape(dataRet.values,np.prod(ddim))

# get raw data
subj = np.repeat(np.arange(0,ddim[0]),[ddim[1]])
time = np.tile(np.linspace(0,1.5,ddim[1]),ddim[0])

# reformat data as frame
dataEnc = pandas.DataFrame(data=np.rot90([signalEnc,subj,time]),columns=['signal','subj','time']) # reshape signal
dataRet = pandas.DataFrame(data=np.rot90([signalRet,subj,time]),columns=['signal','subj','time']) # reshape signal
del(signalEnc,signalRet,subj,time,ddim)

# -- Plot -- #
# define labels and variables
labels = {'ylabel':'Alpha/Beta Power\n(Rem. > Forg.; z)','xlabel':'Time (s)','title':'Encoding'}
variables = {'x':'time','y':'signal'}
limits = {'xline':False,'yline':True,'xlim':[0,1.5],'ylim':[-.4,.4],
          'xticks':[0,0.5,1,1.5],'yticks':[-.4,0,.4],'xprecision':1,'yprecision':1}

# plot encoding
f,ax = pyplot.subplots(1,1)
f.set_figheight(2/2.54) # convert size to inches
f.set_figwidth(3.4/2.54)
f.set_dpi(1000)
eeg_timeseriesplot(dataEnc,variables,limits,labels,colpal['Blue'],ax,[0.4,0.6,0.03])
pyplot.savefig(wdir + "figures/fig3ai.tif",bbox_inches='tight',transparent=True,dpi='figure')

# plot retrieval
f,ax = pyplot.subplots(1,1)
f.set_figheight(2/2.54) # convert size to inches
f.set_figwidth(3.4/2.54)
f.set_dpi(1000)
labels = {'ylabel':'','xlabel':'Time (s)','title':'Retrieval'}
eeg_timeseriesplot(dataRet,variables,limits,labels,colpal['Blue'],ax,[[0.8,1.0,0.03],[1.0,1.2,0.03]])
pyplot.savefig(wdir + "figures/fig3aii.tif",bbox_inches='tight',transparent=True,dpi='figure')
del(labels,limits,variables,dataEnc,dataRet)


# %% -- Figure 3b. -- Neocortical SME Raincloud Contrast -------------------- %
# load raincloud data
rain_Enc = pandas.read_csv(wdir + "data/fig2_atlEncRain.csv", delimiter=",",header=None,names=['rem','forg'])
rain_Ret = pandas.read_csv(wdir + "data/fig2_atlRetRain.csv", delimiter=",",header=None,names=['rem','forg'])
condition = np.append(np.zeros(7),np.ones(5))

# create figure
f,ax = pyplot.subplots(1,1)
f.set_figheight(4.2/2.54) # 4inches 
f.set_figwidth(3.5/2.54) # 12inches
f.set_dpi(1000)

# define labels
labels = {'title':'Encoding (0.4-0.6s)',
          'ylabel':'Alpha/Beta Power (z)',
          'xticklabel':['Recalled','Forgotten'],
          'yticks':[-0.6,-0.3,0,0.3,0.6],
          'yticklabel':['-0.6','','0','','0.6']}

# plot and save
eeg_raincomparison(rain_Enc,[colpal['Blue'],colpal['Gray']],ax,labels,[-0.6,0.6],condition)
pyplot.savefig(wdir + "figures/fig3bi.tif",bbox_inches='tight',transparent=True,dpi='figure')

# create figure
f,ax = pyplot.subplots(1,1)
f.set_figheight(4.2/2.54) # 4inches 
f.set_figwidth(3.5/2.54) # 12inches
f.set_dpi(1000)

# define labels
labels = {'title':'Retrieval (0.8-1.2s)',
          'ylabel':'',
          'xticklabel':['Recalled','Forgotten'],
          'yticks':[-0.6,-0.3,0,0.3,0.6],
          'yticklabel':['-0.6','','0','','0.6']}

# plot and save
eeg_raincomparison(rain_Ret,[colpal['Blue'],colpal['Gray']],ax,labels,[-0.6,0.6],condition)
pyplot.savefig(wdir + "figures/fig3bii.tif",bbox_inches='tight',transparent=True,dpi='figure')


# %% -- Figure 4a. -- Xcorr TimeSeries -------------------------------- #
# -- Prep Data -- #
# load raincloud data
dataEnc = pandas.read_csv(wdir + "data/fig4_encXcSeries.csv", delimiter=",",header=None)
dataRet = pandas.read_csv(wdir + "data/fig4_retXcSeries.csv", delimiter=",",header=None)
dataDiff = pandas.read_csv(wdir + "data/fig4_diffXcSeries.csv", delimiter=",",header=None)
ddim = np.shape(dataEnc)

# normalise signal
signalEnc = np.reshape(dataEnc.values,np.prod(ddim))
signalRet = np.reshape(dataRet.values,np.prod(ddim))
signalDiff = np.reshape(dataDiff.values,np.prod(ddim))

# get raw data
subj = np.repeat(np.arange(0,ddim[0]),[ddim[1]])
time = np.tile(np.linspace(-0.3,0.3,ddim[1]),ddim[0])

# reformat data as frame
dataEnc = pandas.DataFrame(data=np.rot90([signalEnc,subj,time]),columns=['signal','subj','time']) # reshape signal
dataRet = pandas.DataFrame(data=np.rot90([signalRet,subj,time]),columns=['signal','subj','time']) # reshape signal
dataDiff = pandas.DataFrame(data=np.rot90([signalDiff,subj,time]),columns=['signal','subj','time']) # reshape signal
del(signalEnc,signalRet,signalDiff,subj,time,ddim)

# -- Plot -- #
# define labels and variables
labels = {'ylabel':'Cross-Correlation\n(Rem. > Forg.; Fisher z)','xlabel':'Time (s)','title':'Encoding'}
variables = {'x':'time','y':'signal'}
limits = {'xline':True,'yline':True,'xlim':[-0.3,0.3],'ylim':[-.04,.04],
          'xticks':[-0.3,-0.15,0,0.15,0.3],'yticks':[-.04,0,.04],'xprecision':2,'yprecision':2}

# plot encoding
f,ax = pyplot.subplots(1,1)
f.set_figheight(3/2.54) # convert size to inches
f.set_figwidth(3.8/2.54)
f.set_dpi(1000)
eeg_timeseriesplot(dataEnc,variables,limits,labels,colpal['Blue'],ax,[-0.2,-0.1,0.006])
pyplot.savefig(wdir + "figures/fig4ai.tif",bbox_inches='tight',transparent=True,dpi='figure')

# plot retrieval
f,ax = pyplot.subplots(1,1)
f.set_figheight(3/2.54) # convert size to inches
f.set_figwidth(3.8/2.54)
f.set_dpi(1000)
labels = {'ylabel':'Cross-Correlation\n(Rem. > Forg.; Fisher z)','xlabel':'Time (s)','title':'Retrieval'}
eeg_timeseriesplot(dataRet,variables,limits,labels,colpal['Blue'],ax,[0.2,0.3,0.031])
pyplot.savefig(wdir + "figures/fig4aii.tif",bbox_inches='tight',transparent=True,dpi='figure')

# plot diff
f,ax = pyplot.subplots(1,1)
f.set_figheight(3/2.54) # convert size to inches
f.set_figwidth(3.8/2.54)
f.set_dpi(1000)
labels = {'ylabel':'Cross-Correlation\n(Rem. > Forg.; Fisher z)','xlabel':'Time (s)','title':'Enc. > Ret.'}
limits = {'xline':True,'yline':True,'xlim':[-0.3,0.3],'ylim':[-.07,.07],
          'xticks':[-0.3,-0.15,0,0.15,0.3],'yticks':[-.07,0,.07],'xprecision':2,'yprecision':2}
eeg_timeseriesplot(dataDiff,variables,limits,labels,colpal['Blue'],ax,[[-0.2,-0.1,0.006],[0.2,0.3,0.031]])
pyplot.savefig(wdir + "figures/fig4aiii.tif",bbox_inches='tight',transparent=True,dpi='figure')

del(labels,limits,variables,dataEnc,dataRet)


# %% -- Figure 4b. -- Xcorr ANOVA -------------------------------- #
# load raincloud data
data = pandas.read_csv(wdir + "data/fig4_XcANOVA.csv", delimiter=",",header=None)

condition = np.append(np.zeros(7),np.ones(5))
colour = [colpal['Purple'],colpal['Red'],colpal['Purple'],colpal['Red']]

# create figure
f,ax = pyplot.subplots(1,1)
f.set_figheight(3.8/2.54) # 4inches 
f.set_figwidth(3.2/2.54) # 12inches
f.set_dpi(1000)

# define labels
labels = {'title':'',
          'ylabel':'Cross-Correlation\n(Rem. > Forg.; Fisher z)',
          'xticklabel':['Enc.','Ret.'],
          'yticks':[-0.08,-0.04,0,0.04],
          'yticklabel':['-0.08','-0.04','0','0.04'],
          'legend':''}

# plot and save
eeg_anovaboxplot(data,colour,condition,ax,labels,[-0.08,0.04],0.003)
pyplot.savefig(wdir + "figures/fig4b.tif",bbox_inches='tight',transparent=True,dpi='figure')


# %% -- Figure 4c. -- Xcorr Raw -------------------------------- #
# load raincloud data
data = pandas.read_csv(wdir + "data/fig4_encRaw.csv", delimiter=",",header=None)
ddim = np.shape(data)

# normalise signal
signal = data.values
env1 = signal[:,1]
env1 = (env1 - min(env1)) / (max(env1)-min(env1))
env2 = signal[:,4]
env2 = (env2 - min(env2)) / (max(env2)-min(env2))
p1 = np.concatenate((signal[:,0],signal[:,1]))
p2 = np.concatenate((signal[:,3],signal[:,4]))
p3 = np.concatenate((env1,env2))

# get raw data
time = np.tile(np.linspace(0.5,1.1,ddim[0]),2)
condition = np.concatenate((np.zeros(ddim[0]),np.ones(ddim[0])))

# reformat data as frame
data1 = pandas.DataFrame(data=np.rot90([p1,time,condition]),columns=['signal','time','condition']) # reshape signal
data2 = pandas.DataFrame(data=np.rot90([p2,time,condition]),columns=['signal','time','condition']) # reshape signal
data3 = pandas.DataFrame(data=np.rot90([p3,time,condition]),columns=['signal','time','condition']) # reshape signal
del(p1,time)

# -- Plot -- #
# define labels and variables
labels = {'ylabel':'Alpha/Beta\nAmp. (a.u.)','xlabel':'','title':'Encoding'}
variables = {'x':'time','y':'signal','hue':'condition'}
limits = {'xline':False,'yline':False,'xlim':[0.5,1.1],'ylim':[-3,3],
          'xticks':[0.5,0.7,0.9,1.1],'yticks':[-3,3],'xprecision':2,'yprecision':2}

# plot encoding
f,ax = pyplot.subplots(1,1)
f.set_figheight(0.8/2.54) # convert size to inches
f.set_figwidth(4/2.54)
f.set_dpi(1000)
eeg_timeseriesplot(data1,variables,limits,labels,[colpal['Blue'],colpal['Blue']],ax,[0.6,0.85,1])
pyplot.savefig(wdir + "figures/fig4ci.tif",bbox_inches='tight',transparent=True,dpi='figure')

f,ax = pyplot.subplots(1,1)
f.set_figheight(0.8/2.54) # convert size to inches
f.set_figwidth(4/2.54)
f.set_dpi(1000)
labels = {'ylabel':'Gamma\nAmp. (a.u.)','xlabel':'','title':''}
eeg_timeseriesplot(data2,variables,limits,labels,[colpal['Red'],colpal['Red']],ax,[0.75,1,1])
pyplot.savefig(wdir + "figures/fig4cii.tif",bbox_inches='tight',transparent=True,dpi='figure')

f,ax = pyplot.subplots(1,1)
f.set_figheight(0.8/2.54) # convert size to inches
f.set_figwidth(4/2.54)
f.set_dpi(1000)
labels = {'ylabel':'Amp.\n(norm.)','xlabel':'Time (s)','title':''}
limits = {'xline':False,'yline':False,'xlim':[0.5,1.1],'ylim':[-0.05,1.05],
          'xticks':[0.5,0.7,0.9,1.1],'yticks':[0,1],'xprecision':2,'yprecision':2}
eeg_timeseriesplot(data3,variables,limits,labels,[colpal['Gray'],colpal['Gray']],ax)
time = np.linspace(0.5,1.1,ddim[0])
idx1 = (time > 0.6) & (time < 0.85)
idx2 = (time > 0.75) & (time < 1)
pyplot.plot(time[idx1],env1[idx1],figure=f,c=colpal['Blue'])
pyplot.plot(time[idx2],env2[idx2],figure=f,c=colpal['Red'])
pyplot.savefig(wdir + "figures/fig4ciii.tif",bbox_inches='tight',transparent=True,dpi='figure')


# -- RETRIEVAL -- #
# load raincloud data
data = pandas.read_csv(wdir + "data/fig4_retRaw.csv", delimiter=",",header=None)
ddim = np.shape(data)

# normalise signal
signal = data.values
env1 = signal[:,1]
env1 = (env1 - min(env1)) / (max(env1)-min(env1))
env2 = signal[:,4]
env2 = (env2 - min(env2)) / (max(env2)-min(env2))
p1 = np.concatenate((signal[:,0],signal[:,1]))
p2 = np.concatenate((signal[:,3],signal[:,4]))
p3 = np.concatenate((env1,env2))

# get raw data
time = np.tile(np.linspace(0.75,1.35,ddim[0]),2)
condition = np.concatenate((np.zeros(ddim[0]),np.ones(ddim[0])))

# reformat data as frame
data1 = pandas.DataFrame(data=np.rot90([p1,time,condition]),columns=['signal','time','condition']) # reshape signal
data2 = pandas.DataFrame(data=np.rot90([p2,time,condition]),columns=['signal','time','condition']) # reshape signal
data3 = pandas.DataFrame(data=np.rot90([p3,time,condition]),columns=['signal','time','condition']) # reshape signal
del(p1,time)

# -- Plot -- #
# define labels and variables
labels = {'ylabel':'','xlabel':'','title':'Retrieval'}
variables = {'x':'time','y':'signal','hue':'condition'}
limits = {'xline':False,'yline':False,'xlim':[0.75,1.35],'ylim':[-5,5],
          'xticks':[0.75,0.95,1.15,1.35],'yticks':[-5,5],'xprecision':2,'yprecision':2}

# plot encoding
f,ax = pyplot.subplots(1,1)
f.set_figheight(0.8/2.54) # convert size to inches
f.set_figwidth(4/2.54)
f.set_dpi(1000)
eeg_timeseriesplot(data1,variables,limits,labels,[colpal['Blue'],colpal['Blue']],ax,[0.9,1.15,1])
pyplot.savefig(wdir + "figures/fig4cvi.tif",bbox_inches='tight',transparent=True,dpi='figure')

f,ax = pyplot.subplots(1,1)
f.set_figheight(0.8/2.54) # convert size to inches
f.set_figwidth(4/2.54)
f.set_dpi(1000)
labels = {'ylabel':'','xlabel':'','title':''}
eeg_timeseriesplot(data2,variables,limits,labels,[colpal['Purple'],colpal['Purple']],ax,[0.85,1,1])
pyplot.savefig(wdir + "figures/fig4cvii.tif",bbox_inches='tight',transparent=True,dpi='figure')

f,ax = pyplot.subplots(1,1)
f.set_figheight(0.8/2.54) # convert size to inches
f.set_figwidth(4/2.54)
f.set_dpi(1000)
labels = {'ylabel':'','xlabel':'Time (s)','title':''}
limits = {'xline':False,'yline':False,'xlim':[0.75,1.35],'ylim':[-0.05,1.05],
          'xticks':[0.75,0.95,1.15,1.35],'yticks':[0,1],'xprecision':2,'yprecision':2}
eeg_timeseriesplot(data3,variables,limits,labels,[colpal['Gray'],colpal['Gray']],ax)
time = np.linspace(0.75,1.35,ddim[0])
idx1 = (time > 0.9) & (time < 1.15)
idx2 = (time > 0.85) & (time < 1)
pyplot.plot(time[idx1],env1[idx1],figure=f,c=colpal['Blue'])
pyplot.plot(time[idx2],env2[idx2],figure=f,c=colpal['Purple'])
pyplot.savefig(wdir + "figures/fig4cviii.tif",bbox_inches='tight',transparent=True,dpi='figure')


# %% -- Figure 4d. -- Neocortical SME Raincloud Contrast -------------------- %
# load raincloud data
rain_Enc = pandas.read_csv(wdir + "data/fig4_EncRain.csv", delimiter=",",header=None,names=['rem','forg'])
rain_Ret = pandas.read_csv(wdir + "data/fig4_RetRain.csv", delimiter=",",header=None,names=['rem','forg'])
condition = np.append(np.zeros(7),np.ones(5))

# create figure
f,ax = pyplot.subplots(1,1)
f.set_figheight(3.2/2.54) # 4inches 
f.set_figwidth(3/2.54) # 12inches
f.set_dpi(1000)

# define labels
labels = {'title':'Encoding (-0.2s to -0.1s)',
          'ylabel':'Cross-Correlation (Fisher z)',
          'xticklabel':['Recalled','Forgotten'],
          'yticks':[-0.05,-0.025,0,0.025,0.05],
          'yticklabel':['-0.05','','0','','0.05']}

# plot and save
eeg_raincomparison(rain_Enc,[colpal['Blue'],colpal['Gray']],ax,labels,[-0.05,0.05],condition)
pyplot.savefig(wdir + "figures/fig4di.tif",bbox_inches='tight',transparent=True,dpi='figure')

# create figure
f,ax = pyplot.subplots(1,1)
f.set_figheight(3.2/2.54) # 4inches 
f.set_figwidth(3/2.54) # 12inches
f.set_dpi(1000)

# define labels
labels = {'title':'Retrieval (0.2s to 0.3s)',
          'ylabel':'Cross-Correlation (Fisher z)',
          'xticklabel':['Recalled','Forgotten'],
          'yticks':[-0.05,-0.025,0,0.025,0.05],
          'yticklabel':['-0.05','','0','','0.05']}

# plot and save
eeg_raincomparison(rain_Ret,[colpal['Blue'],colpal['Gray']],ax,labels,[-0.05,0.05],condition)
pyplot.savefig(wdir + "figures/fig4dii.tif",bbox_inches='tight',transparent=True,dpi='figure')


# %% -- Supplementary 1a ---------------------------------------------------- %
# -- Prep Data -- #
# load raincloud data
data = pandas.read_csv(wdir + "data/fig2_encRetGammaSpecHits.csv", delimiter=",",header=None)
ddim = np.shape(data)

# normalise signal
signal = np.reshape(data.values,np.prod(ddim))

# get raw data
subj = np.repeat(np.arange(0,ddim[0]),[ddim[1]])
freq = np.tile(np.linspace(30,100,ddim[1]),ddim[0])

# reformat data as frame
data = pandas.DataFrame(data=np.rot90([signal,subj,freq]),columns=['signal','subj','freq']) # reshape signal
del(ddim,subj,freq,signal)

# -- Plot -- #
# create figure
f,ax = pyplot.subplots(1,1)
f.set_figheight(3/3) # convert size to inches
f.set_figwidth(4.4/3)
f.set_dpi(1000)

# define labels and variables
labels = {'ylabel':'1/f Corrected Power\n(Enc. > Ret.; a.u.)','xlabel':'Frequency (Hz.)'}
variables = {'x':'freq','y':'signal'}
limits = {'xline':False,'yline':True,'xlim':[30,100],'ylim':[-.06,.06],
          'xticks':np.arange(40,100,10),'yticks':[-.06,0,.06],'xprecision':3,'yprecision':3}

# plot
eeg_timeseriesplot(data,variables,limits,labels,colpal['Red'],ax,[[40,50,0.02],[60,70,0.001],[70,80,0.02]])
pyplot.savefig(wdir + "figures/sup1a.tif",bbox_inches='tight',transparent=True,dpi='figure')
del(data,labels,variables,limits)


# %% -- Supplementary 1b ---------------------------------------------------- %
# load raincloud data
data = pandas.read_csv(wdir + "data/fig2_encRetGammaSpecMisses.csv", delimiter=",",header=None)
ddim = np.shape(data)

# normalise signal
signal = np.reshape(data.values,np.prod(ddim))

# get raw data
subj = np.repeat(np.arange(0,ddim[0]),[ddim[1]])
freq = np.tile(np.linspace(30,100,ddim[1]),ddim[0])

# reformat data as frame
data = pandas.DataFrame(data=np.rot90([signal,subj,freq]),columns=['signal','subj','freq']) # reshape signal
del(ddim,subj,freq,signal)

# -- Plot -- #
# create figure
f,ax = pyplot.subplots(1,1)
f.set_figheight(3/3) # convert size to inches
f.set_figwidth(4.4/3)
f.set_dpi(1000)

# define labels and variables
labels = {'ylabel':'','xlabel':'Frequency (Hz.)'}
variables = {'x':'freq','y':'signal'}
limits = {'xline':False,'yline':True,'xlim':[30,100],'ylim':[-.06,.06],
          'xticks':np.arange(40,100,10),'yticks':[-.06,0,.06],'xprecision':3,'yprecision':3}

# plot
eeg_timeseriesplot(data,variables,limits,labels,colpal['Gray'],ax)
pyplot.savefig(wdir + "figures/sup1b.tif",bbox_inches='tight',transparent=True,dpi='figure')
del(data,labels,variables,limits)


# %% -- Supplementary 1c ---------------------------------------------------- %# load raincloud data
data1 = pandas.read_csv(wdir + "data/fig2_encGammaSpecHits.csv", delimiter=",",header=None)
data2 = pandas.read_csv(wdir + "data/fig2_retGammaSpecHits.csv", delimiter=",",header=None)
ddim = np.shape(data1)

# normalise signal
signal1 = np.reshape(data1.values,np.prod(ddim))
signal2 = np.reshape(data2.values,np.prod(ddim))
signal  = np.concatenate((signal1,signal2))

# get raw data
subj = np.tile(np.repeat(np.arange(0,ddim[0]),[ddim[1]]),2)
freq = np.tile(np.linspace(30,100,ddim[1]),ddim[0]*2)
cond = np.concatenate((np.zeros(np.prod(ddim)),np.ones(np.prod(ddim))))

# reformat data as frame
data = pandas.DataFrame(data=np.rot90([signal,subj,freq,cond]),columns=['signal','subj','freq','cond']) # reshape signal
del(ddim,subj,freq,signal,signal1,signal2)

# -- Plot -- #
# create figure
f,ax = pyplot.subplots(1,1)
f.set_figheight(3/3) # convert size to inches
f.set_figwidth(4.4/3)
f.set_dpi(1000)

# define labels and variables
labels = {'ylabel':'1/f Corrected Power (a.u.)','xlabel':'Frequency (Hz.)'}
variables = {'x':'freq','y':'signal','hue':'cond'}
limits = {'xline':False,'yline':True,'xlim':[30,100],'ylim':[1,1],
          'xticks':np.arange(40,100,10),'yticks':[1,1.5,2],'xprecision':3,'yprecision':3}

# plot
eeg_timeseriesplot(data,variables,limits,labels,[colpal['Red'],colpal['Gray']],ax)
pyplot.savefig(wdir + "figures/sup1c.tif",bbox_inches='tight',transparent=True,dpi='figure')
del(data,labels,variables,limits)


# %% -- Supplementary 1d ---------------------------------------------------- %# load raincloud data
data1 = pandas.read_csv(wdir + "data/fig2_encGammaSpecMisses.csv", delimiter=",",header=None)
data2 = pandas.read_csv(wdir + "data/fig2_retGammaSpecMisses.csv", delimiter=",",header=None)
ddim = np.shape(data1)

# normalise signal
signal1 = np.reshape(data1.values,np.prod(ddim))
signal2 = np.reshape(data2.values,np.prod(ddim))
signal  = np.concatenate((signal1,signal2))

# get raw data
subj = np.tile(np.repeat(np.arange(0,ddim[0]),[ddim[1]]),2)
freq = np.tile(np.linspace(30,100,ddim[1]),ddim[0]*2)
cond = np.concatenate((np.zeros(np.prod(ddim)),np.ones(np.prod(ddim))))

# reformat data as frame
data = pandas.DataFrame(data=np.rot90([signal,subj,freq,cond]),columns=['signal','subj','freq','cond']) # reshape signal
del(ddim,subj,freq,signal,signal1,signal2)

# -- Plot -- #
# create figure
f,ax = pyplot.subplots(1,1)
f.set_figheight(3/3) # convert size to inches
f.set_figwidth(4.4/3)
f.set_dpi(1000)

# define labels and variables
labels = {'ylabel':'','xlabel':'Frequency (Hz.)'}
variables = {'x':'freq','y':'signal','hue':'cond'}
limits = {'xline':False,'yline':True,'xlim':[30,100],'ylim':[1,1],
          'xticks':np.arange(40,100,10),'yticks':[1,1.5,2],'xprecision':3,'yprecision':3}

# plot
eeg_timeseriesplot(data,variables,limits,labels,[colpal['Red'],colpal['Gray']],ax)
pyplot.savefig(wdir + "figures/sup1d.tif",bbox_inches='tight',transparent=True,dpi='figure')
del(data,labels,variables,limits)


# %% -- Supplementary 1e ---------------------------------------------------- %# load raincloud data
data1 = pandas.read_csv(wdir + "data/fig2_encGammaSpecHits.csv", delimiter=",",header=None)
data2 = pandas.read_csv(wdir + "data/fig2_retGammaSpecHits.csv", delimiter=",",header=None)
ddim = np.shape(data1)

# normalise signal
signal1 = np.reshape(data1.values,np.prod(ddim))
signal2 = np.reshape(data2.values,np.prod(ddim))
signal  = np.concatenate((signal1,signal2))

# get raw data
subj = np.tile(np.repeat(np.arange(0,ddim[0]),[ddim[1]]),2)
freq = np.tile(np.linspace(30,100,ddim[1]),ddim[0]*2)
cond = np.concatenate((np.zeros(np.prod(ddim)),np.ones(np.prod(ddim))))

# reformat data as frame
data = pandas.DataFrame(data=np.rot90([signal,subj,freq,cond]),columns=['signal','subj','freq','cond']) # reshape signal
del(ddim,subj,freq,signal,signal1,signal2)

# -- Plot -- #
# create figure
f,ax = pyplot.subplots(1,1)
f.set_figheight(3/3) # convert size to inches
f.set_figwidth(4.4/3)
f.set_dpi(1000)

# define labels and variables
labels = {'ylabel':'1/f Corrected Power (a.u.)','xlabel':'Frequency (Hz.)'}
variables = {'x':'freq','y':'signal','hue':'cond'}
limits = {'xline':False,'yline':True,'xlim':[38,72],'ylim':[1.2,2],
          'xticks':np.arange(40,80,10),'yticks':[1.2,1.6,2],'xprecision':3,'yprecision':3}

# plot
eeg_timeseriesplot(data,variables,limits,labels,[colpal['Red'],colpal['Gray']],ax)
pyplot.savefig(wdir + "figures/sup1e.tif",bbox_inches='tight',transparent=True,dpi='figure')
del(data,labels,variables,limits)
