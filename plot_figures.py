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
from matplotlib import transforms
import matplotlib.font_manager as font_manager

# %% -- define functions ---------------------------------------------------- #
# plot raincloud
def custom_rainplot(data,colour,axes,fontname,labels,ylim,offset,pvalue):
    
    # get transform data
    trans   = axes.transData
    offset  = transforms.ScaledTranslation(offset,0,f.dpi_scale_trans)
    
    # plot violin
    axes=pt.half_violinplot(data = data,bw = "scott",inner = None, scale = "count",
                          width = 0.5, linewidth = 1, cut = 1, palette = colour,
                          ax = axes, edgecolor = [0,0,0])
    
    # plot single points
    axes=sns.swarmplot(data = data, edgecolor =[0,0,0], size = 1.5, 
                       transform = trans + offset, palette = colour,
                       ax = axes)
    
    # plot mean and confidence intervals
    axes=sns.boxplot(data = data, palette = colour, width = 0.1, ax = axes, linewidth = 1, fliersize = 1)
    
    # plot significance    
    sig_offset = ylim[1]-(ylim[1]*0.05)
    for i in np.arange(0,np.size(pvalue)):
        if pvalue[i] < 0.001:
            pyplot.scatter(np.array([-0.05,0,0.05])+i,[sig_offset,sig_offset,sig_offset],s=1,c=[0.8,0.8,0.8],marker='*',edgecolors=None)
        elif pvalue[i] < 0.01:
            pyplot.scatter(np.array([-0.025,0.025])+i,[sig_offset,sig_offset],s=1,c=[0.8,0.8,0.8],marker='*',edgecolors=None)
        elif pvalue[i] < 0.05:
            pyplot.scatter(np.array([0])+i,[sig_offset],s=2,c='black',marker='*',edgecolors=None)
    
    # add horizontal line
    axes.axhline(y=0, xmin=-1, xmax=3, color=[0,0,0], linewidth = 1)
    
    # aesthetics
    axes.set_ylabel(labels['ylabel'],fontname=fontname,fontsize=7,labelpad=5,fontweight='light')   # add Y axis label
    axes.set_ylim(ylim)                  # set Y axis range to 0 -> 1
    axes.set_xlim(-0.65,-0.5+len(data.columns))                  # set Y axis range to 0 -> 1
    axes.tick_params(axis='x',          # change X tick parameters
                   which='both',          # affect both major and minor ticks
                   bottom=False,          # turn off bottom ticks
                   labelbottom=True,  # keep bottom labels
                   pad=2.5,
                   width=1)   
    axes.tick_params(axis='y',          # change X tick parameters
                       pad=3,
                       width=1,
                       length=2.5)
    axes.set_yticks(labels['yticks'])
    axes.set_xticklabels(labels['xticklabel'],fontname=fontname,fontweight='light',fontsize=6)
    axes.set_yticklabels(labels['yticklabel'],fontname=fontname,fontweight='light',fontsize=6)

    # change axes
    axes.spines['top'].set_visible(False)
    axes.spines['right'].set_visible(False)
    axes.spines['bottom'].set_visible(False)
    axes.spines['left'].set_linewidth(1)
    
               
def eeg_timeseriesplot(data,variables,limits,labels,colour,ax):
 
   # plot multicondition timeseries
   if ('hue' in variables):
      sns.lineplot(x=variables['x'],
                   y=variables['y'],
                   data=data,
                   hue=variables['hue'],
                   hue_order=[0,1],
                   palette=colour,
                   ax = ax,
                   linewidth=1)
      
      # add leened
      font = font_manager.FontProperties(family='Arial',weight='light',size=7)
      ax.legend(prop=font)
      ax.legend(labels['legend'],frameon=False,fontsize=7,bbox_to_anchor=(0., 1.02, 1., .102), 
                  loc=3, ncol=2, mode="expand", borderaxespad=0.,prop=font)
      
   # plot single condition timeseries   
   else:
      sns.lineplot(x=variables['x'],
                   y=variables['y'],
                   data=data,
                   palette=colour,
                   ax = ax,
                   linewidth=1)
      #ax.get_legend().remove()
               
   # add horizontal line
   if limits['xline']:
      ax.axvline(x=0, linewidth = 1, color = [0,0,0], linestyle='--')
   if limits['yline']:
      ax.axhline(y=0, linewidth = 1, color = [0,0,0], linestyle='-')
   
   # set x scale limits
   if ('xlim' in limits):   
      
      # set limits
      ax.set_xlim(limits['xlim'])
         
   # set y scale limits
   if ('ylim' in limits):  
      
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
      ax.set_xticklabels(limits['xticks'],fontname='Arial',fontweight='light',fontsize=7)
         
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
      ax.set_yticklabels(limits['yticks'],fontname='Arial',fontweight='light',fontsize=7)
         
   
   # aesthetics
   ax.set_ylabel(labels['ylabel'],fontname='Arial',fontsize=7,labelpad=2,fontweight='light')   # add Y axis label
   ax.set_xlabel(labels['xlabel'],fontname='Arial',fontsize=7,labelpad=2,fontweight='light')   # add Y axis label
   ax.tick_params(axis='both',          # change X tick parameters
                      pad=3,
                      length=2.5)
   
   # change axes
   ax.spines['top'].set_visible(False)
   ax.spines['right'].set_visible(False)

# %% -- define key parameters ----------------------------------------------- #
# get current directory
wdir = 'C:/Users/bengr/Documents/Github/hippcortex-interactions/'

# set context
sns.set_context("paper")

# set plot defaults
mpl.rcParams['xtick.major.width'] = 1
mpl.rcParams['ytick.major.width'] = 1
mpl.rcParams['xtick.color'] = [0,0,0]
mpl.rcParams['ytick.color'] = [0,0,0]
mpl.rcParams['lines.linewidth'] = 1
mpl.rcParams['axes.linewidth'] = 1

# get colour dict
RdBu = sns.color_palette("RdBu_r", 10)
colour = {'Red':RdBu[8],
          'Blue':RdBu[1]}

# %% -- Figure 2a ----------------------------------------------------------- #
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

# -- Plot -- #
# create figure
f,ax = pyplot.subplots(1,1)
f.set_figheight(3.33/2.54) # convert size to inches
f.set_figwidth(5.4/2.54)
f.set_dpi(1000)

# define labels and variables
labels = {'legend':['ATL','Hippocampus'],'ylabel':'1/f Corrected Power (z)','xlabel':'Frequency (Hz.)'}
variables = {'x':'freq','y':'signal','hue':'condition'}
limits = {'xline':False,'yline':False,'xlim':[1.5,64],'ylim':[-1,4],
          'xticks':[4,8,16,32,64],'yticks':[-2,0,2,4],'xprecision':1,'yprecision':1}

# plot
eeg_timeseriesplot(data,variables,limits,labels,[colour['Red'],colour['Blue']],ax)

# %% -- Figure 2b ----------------------------------------------------------- #
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

# -- Plot -- #
# create figure
f,ax = pyplot.subplots(1,1)
f.set_figheight(3.33/2.54) # convert size to inches
f.set_figwidth(5.4/2.54)
f.set_dpi(1000)

# define labels and variables
labels = {'ylabel':'1/f Corrected Power\n(Enc. > Ret.; a.u.)','xlabel':'Frequency (Hz.)'}
variables = {'x':'freq','y':'signal'}
limits = {'xline':False,'yline':True,'xlim':[32,98],'ylim':[-.05,.05],
          'xticks':np.arange(40,100,10),'yticks':[-.05,0,.05],'xprecision':3,'yprecision':3}

# plot
eeg_timeseriesplot(data,variables,limits,labels,colour['Red'],ax)

# %% 
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

# -- Plot -- #
# define labels and variables
labels = {'ylabel':'1/f Corrected Power\n(Rem. > Rog.; z)','xlabel':'Time (s)'}
variables = {'x':'time','y':'signal'}
limits = {'xline':False,'yline':True,'xlim':[0,1.5],'ylim':[-.4,.4],
          'xticks':[0,0.5,1,1.5],'yticks':[-.4,0,.4],'xprecision':1,'yprecision':1}

# plot encoding
f,ax = pyplot.subplots(1,1)
f.set_figheight(3.33/2.54) # convert size to inches
f.set_figwidth(5.4/2.54)
f.set_dpi(1000)
eeg_timeseriesplot(dataEnc,variables,limits,labels,colour['Blue'],ax)

# plot retrieval
f,ax = pyplot.subplots(1,1)
f.set_figheight(3.33/2.54) # convert size to inches
f.set_figwidth(5.4/2.54)
f.set_dpi(1000)
eeg_timeseriesplot(dataRet,variables,limits,labels,colour['Blue'],ax)

































# %% -- prepare data -------------------------------------------- #
# define raincloud data filename
data_fname = wdir + "data/fig2_45HzGammaRain.csv"

# load raincloud data
rain_45Hz = pandas.read_csv(wdir + "data/fig2_45HzGammaRain.csv", delimiter=",",header=None,names=['encoding','retrieval'])
rain_70Hz = pandas.read_csv(wdir + "data/fig2_70HzGammaRain.csv", delimiter=",",header=None,names=['encoding','retrieval'])

# create figure
f,ax = pyplot.subplots(1,1)
f.set_figheight(6.2/2.54) # 4inches 
f.set_figwidth(5.2/2.54) # 12inches
f.set_dpi(200)

# define colour scheme
colour = sns.color_palette("Blues",n_colors=7)
colour = [colour[3]]

# define labels
labels = {'title':'',
          'ylabel':'"Slow" Gamma Power (a.u.)',
          'xticklabel':[''],
          'yticks':[1.1,1.2,1.3,1.4,1.5],
          'yticklabel':['1.1','1.2','1.3','1.4','1.5']}

# plot raincloud
custom_rainplot(rain_45Hz,colour,ax,'Calibri',labels,[1.11,1.5],0.15,[1])
   
# create figure
f,ax = pyplot.subplots(1,1)
f.set_figheight(6.2/2.54) # 4inches 
f.set_figwidth(5.2/2.54) # 12inches
f.set_dpi(200)

# define colour scheme
colour = sns.color_palette("Blues",n_colors=7)
colour = [colour[3]]

# define labels
labels = {'title':'',
          'ylabel':'"Fast" Gamma Power (a.u.)',
          'xticklabel':[''],
          'yticks':[1.2,1.3,1.4,1.5,1.6,1.7],
          'yticklabel':['1.2','1.3','1.4','1.5','1.6','1.7']}

# plot raincloud
custom_rainplot(rain_70Hz,colour,ax,'Calibri',labels,[1.25,1.65],0.15,[1])





# -- prep frequency
# load frequency data
datatmp = pandas.read_csv(wdir + "data/fig2_data/group_task-memory_eeg-freqseries.csv",
                                 delimiter=',',
                                 header=None)

# create new structure for frequency data
data_frequency = pandas.DataFrame(data=np.reshape(datatmp.values,[datatmp.size]),columns=['signal'])

# create subject number array
data_frequency = data_frequency.assign(subj=pandas.Series(np.repeat(np.arange(0,21),[150])).values)
    
# create condition array
data_frequency = data_frequency.assign(condition=pandas.Series(np.tile(np.append(np.zeros(75),np.ones(75)),[21])).values)

# create frequency array
data_frequency = data_frequency.assign(frequency=pandas.Series(np.tile(np.linspace(3,40,75),[42])).values)


# -- prep frequency diff
# get frames for hits and misses
data_freqA = data_frequency[data_frequency['condition']==1];
data_freqB = data_frequency[data_frequency['condition']==0];
data_freqA = data_freqA.reset_index()
data_freqB = data_freqB.reset_index()

# create new frame
data_freqdiff = data_freqA
data_freqdiff['signal'] = data_freqB['signal']-data_freqA['signal']


# -- prep time
# load frequency data
datatmp = pandas.read_csv(wdir + "data/fig2_data/group_task-memory_eeg-timeseries.csv",
                                 delimiter=',',
                                 header=None)

# create new structure for frequency data
data_time = pandas.DataFrame(data=np.reshape(datatmp.values,[datatmp.size]),columns=['signal'])

# create subject number array
data_time = data_time.assign(subj=pandas.Series(np.repeat(np.arange(0,21),[122])).values)
    
# create condition array
data_time = data_time.assign(condition=pandas.Series(np.tile(np.append(np.zeros(61),np.ones(61)),[21])).values)

# create frequency array
data_time = data_time.assign(frequency=pandas.Series(np.tile(np.linspace(-1,2,61),[42])).values)


# -- prep frequency diff
# get frames for hits and misses
data_timeA = data_time[data_time['condition']==1];
data_timeB = data_time[data_time['condition']==0];
data_timeA = data_timeA.reset_index()
data_timeB = data_timeB.reset_index()

# create new frame
data_timediff = data_timeA
data_timediff['signal'] = data_timeB['signal']-data_timeA['signal']


# %% ----- Raincloud Plot ----- # 
# create figure
f,ax = pyplot.subplots(1,1)
f.set_figheight(6.2/2.54) # 4inches 
f.set_figwidth(3.1/2.54) # 12inches
f.set_dpi(1000)

# define colour scheme
colour = sns.color_palette("Blues",n_colors=7)
colour = [colour[3]]

# define labels
labels = {'title':'',
          'ylabel':'Alpha/Beta Power (Rem. > Forgot.; z)',
          'xticklabel':[''],
          'yticks':[-0.5,-0.25,0,0.25],
          'yticklabel':['-0.5','-0.25','0','0.25']}

# plot raincloud
custom_rainplot(data_raincloud,colour,ax,'Calibri',labels,[-0.5,0.25],0.15,[1])
   
# save image
pyplot.savefig(wdir + "/figures/fig2a.tif",bbox_inches='tight',transparent=True,dpi='figure')
  
# %% ----- Frequency Individual Series ----- # 
# create figure
f,ax = pyplot.subplots(1,1)
f.set_figheight(2.5/2.54) # 4inches 
f.set_figwidth(3.8/2.54) # 12inches
f.set_dpi(1000)

# define colour scheme
colour = sns.color_palette("Blues",n_colors=7)
colour = [(0.7,0.7,0.7),colour[5]]

# define labels and variables
labels = {'legend':['Forgot.','Rem.'],
          'ylabel':'Power (z)',
          'xlabel':'Frequency (Hz.)'}

variables = {'x':'frequency',
             'y':'signal',
             'hue':'condition'}

# plot frequency series
custom_timeseriesplot(data_frequency,variables,ax,colour,labels,[3,40],[-0.25,0.25],[5,10,15,20,25,30,35,40],True,False)

# save image
pyplot.savefig(wdir + "/figures/fig2b.tif",bbox_inches='tight',transparent=True,dpi='figure')

# %% ----- Frequency Difference Series ----- # 
# create figure
f,ax = pyplot.subplots(1,1)
f.set_figheight(2.5/2.54) # 4inches 
f.set_figwidth(3.8/2.54) # 12inches
f.set_dpi(1000)

# define colour scheme
colour = sns.color_palette("Blues",n_colors=7)
colour = [colour[5],colour[5]]

# define labels and variables
labels = {'legend':[''],
          'ylabel':'Power\n(Rem. > Forg.; z)',
          'xlabel':'Frequency (Hz.)'}

variables = {'x':'frequency',
             'y':'signal',
             'hue':'condition'}

# plot frequency series
custom_timeseriesplot(data_freqdiff,variables,ax,colour,labels,[3,40],[-0.25,0.1],[5,10,15,20,25,30,35,40],False,True,True)

# save image
pyplot.savefig(wdir + "/figures/fig2c.tif",bbox_inches='tight',transparent=True,dpi='figure')

# %% ----- Time Series ----- # 
# create figure
f,ax = pyplot.subplots(1,1)
f.set_figheight(2.5/2.54) # 4inches 
f.set_figwidth(3.8/2.54) # 12inches
f.set_dpi(1000)

# define colour scheme
colour = sns.color_palette("Blues",n_colors=7)
colour = [(0.7,0.7,0.7),colour[5]]

# define labels and variables
labels = {'legend':['Forgot.','Rem.'],
          'ylabel':'Power (z)',
          'xlabel':'Time (s)'}

variables = {'x':'frequency',
             'y':'signal',
             'hue':'condition'}

# plot frequency series
custom_timeseriesplot(data_time,variables,ax,colour,labels,[-0.5,2],[-0.25,0.25],[-0.5,0,0.5,1,1.5,2],True,True)

# save image
pyplot.savefig(wdir + "/figures/fig2d.tif",bbox_inches='tight',transparent=True,dpi='figure')

# %% ----- Time Series Individual ----- # 
# create figure
f,ax = pyplot.subplots(1,1)
f.set_figheight(2.5/2.54) # 4inches 
f.set_figwidth(3.8/2.54) # 12inches
f.set_dpi(1000)

# define colour scheme
colour = sns.color_palette("Blues",n_colors=7)
colour = [colour[5],colour[5]]

# define labels and variables
labels = {'legend':[''],
          'ylabel':'Power\n(Rem. > Forg.; z)',
          'xlabel':'Time (s)'}

variables = {'x':'frequency',
             'y':'signal',
             'hue':'condition'}

# plot frequency series
custom_timeseriesplot(data_timediff,variables,ax,colour,labels,[-0.5,2],[-0.3,0.2],[-0.5,0,0.5,1,1.5,2],False,True,True)

# save image
pyplot.savefig(wdir + "/figures/fig2e.tif",bbox_inches='tight',transparent=True,dpi='figure')