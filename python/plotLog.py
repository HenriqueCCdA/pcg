import matplotlib.pyplot as plt 
# **********************************************************************
def plot(name):
  """
**********************************************************************
  """
# ... leitura 
  with open(name,'r') as f:
    lines = f.read().split('\n')
# ......................................................................

# ... retira a ultima linha
  lines = lines[:-1]
# ......................................................................

# ...
  it = []
  r  = []
  for line in lines:
    it.append(int(line.split(',')[0]))
    r.append(float(line.split(',')[1]))  
# ......................................................................
  
# ...
  plt.semilogy(it,r,label=name.split('.')[0])
# ......................................................................
# **********************************************************************
#
# **********************************************************************
def initPlot(iCod=0):
  """
**********************************************************************
  """  
  fig = plt.figure(iCod,figsize=(6,6))
  ax  = fig.add_subplot(111)
# **********************************************************************
#
# **********************************************************************
def show():
  """
**********************************************************************
  """
  plt.xlim(0)
  plt.ylim(1.e-16,10.0)
  plt.legend(bbox_to_anchor=(0.75, 0.95), loc=2, borderaxespad=0.
            ,prop={'size':8})
  plt.ylabel(r'${r/r_0}$')
  plt.xlabel(r'$Iteration$')
  plt.grid(color='gray', linestyle='--', linewidth=0.25)
  plt.show()
# **********************************************************************
