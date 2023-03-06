import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams['contour.negative_linestyle'] = 'solid'

def plot_streamlines(simulation_name,Grid,v,PSI,psi_min,psi_max,label):

    # x velocity internpolation from faces to cell centers
    Xp,Yp = np.meshgrid(Grid.Vx.xc,Grid.Vy.yc)
    fig = plt.figure(figsize=(15,15) , dpi=100)
    CS = plt.contour(Xp, Yp, PSI,100,colors='k')
    
    if label=='label_yes':
        plt.clabel(CS, fontsize=9, inline=1)
    
    plt.xlabel('x')
    plt.ylabel('y')
    plt.xlim([Grid.p.xmin,Grid.p.xmax])
    plt.ylim([Grid.p.ymin,Grid.p.ymax])
    plt.axis('scaled')    
    plt.savefig(f'{simulation_name}.pdf')  


