

from matplotlib.ticker import ScalarFormatter
import math
import cPickle as pickle
import runner 
import matplotlib.pylab as plt
import numpy as np
import downSamplePickles as dsp
from runShannonTest import *
import json
import taufitting as tf
import scipy.fftpack as fftp
runner.init()
import os.path

class empty:pass
mM_to_uM = 1e3
ms_to_s = 1e-3

plotlyAuth=None

# For plotly support 
def PlotViaPlotly(casesSubset,state):
  downsample=10
  import plotly
  if plotlyAuth==None:
    with open('/net/share/pmke226/PLOTLY', 'r') as f:
        first_line = f.readline()
    plotlyKey = first_line
    plotly.tools.set_credentials_file(
      username='huskeypm', api_key='3x5rz5d19r') #plotlyKey)
  import plotly.tools as tls
  import plotly.plotly as py

  fig = plt.figure()
  ax = plt.subplot(111)
  case1=casesSubset[0]
  case2=casesSubset[1]
  title ="%s_%s_%s"%(state,case1.label,case2.label)
  ax.set_title("%s: %s,%s"%(state,case1.label,case2.label))
  for i,case in enumerate(casesSubset):
    pkg = GetData(case.data,state)
    ax.plot(pkg.t[::downsample],pkg.valsIdx[::downsample], label=case.label)
  plotly_fig = tls.mpl_to_plotly( fig )

  # Adding custom attributes to legend
  plotly_fig['layout']['showlegend'] = True
#layout = go.Layout(
#    xaxis=dict(


#        range=[2, 5]
#    )
#)
#fig = go.Figure(data=data, layout=layout)

  plot_url = py.iplot(plotly_fig, filename = title)
  print plot_url.resource

### 
### I/O 
###








#-----------------------------------------------Here begins relevant stuff-------------------------------------------
def makePackage(p,p_idx,s,s_idx,j,j_idx,t):

  return {'p':p,'s':s,'t':t,'j':np.asarray(j),\
           'p_idx':p_idx,'s_idx':s_idx,'j_idx':j_idx}
  
def writePickle(name,p,p_idx,s,s_idx,j,j_idx,t):
  # store to pickle
  # using 'asarray' since my 'j' was getting stored as its transpose 
  #data1 = {'p':p,'s':s,'t':t,'j':np.asarray(j),\
  #         'p_idx':p_idx,'s_idx':s_idx,'j_idx':j_idx}
  
  data1 = makePackage(p,p_idx,s,s_idx,j,j_idx,t)

  #print "j again: ", len(j) 
  #print "j_idx: ",np.shape(j_idx)
  print "name: ", name
  if ".pkl" not in name[-4:]:
    name += ".pkl"
  output = open(name, 'wb')
  pickle.dump(data1, output)
  output.close()
  print "SUCCESS! Wrote output to", name

def readPickle(name = "PCa0.75kss0.25.pkl",verbose=True):          
  if verbose: 
    print "Reading " + name  
  pkl_file = open(name, 'rb')
  data1 = pickle.load(pkl_file)  
  pkl_file.close()

  return data1  

def LoadPickles(caseDict,noOverwrite=False,verbose=True):
  for key,case in caseDict.iteritems():
    if verbose:
      print "# ", key
      print "Loading"  , case.name

    if hasattr(case,'data') and noOverwrite==True:
