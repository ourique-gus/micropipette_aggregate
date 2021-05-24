
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import matplotlib.ticker as ticker
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm
import os


def plot_SimDescription(path, name):
    ax.cla()
    
    #ax.xaxis.set_major_locator(ticker.MaxNLocator(4))
    #ax.yaxis.set_major_locator(ticker.MaxNLocator(4))

    x_minorLocator = ticker.AutoMinorLocator(4)
    y_minorLocator = ticker.AutoMinorLocator(4)
    ax.xaxis.set_minor_locator(x_minorLocator)
    ax.yaxis.set_minor_locator(y_minorLocator)

    ax.tick_params(axis='both', which='major', direction='in', top = True, right=True)
    ax.tick_params(axis='both', which='minor', direction='in', top=True,right=True)


    #ax.set_xlabel('$x$')
    #ax.set_ylabel('$y$')

    majorFormatter = FormatStrFormatter('')
    ax.yaxis.set_major_formatter(majorFormatter)
    ax.xaxis.set_major_formatter(majorFormatter)


    
    dlim_var=0
    dlim_var_y=15
    

    c_cell='C3'
    c_wall='C0'
    c_pres='black'
    c_gray='gray'
    



    read=np.load(path, allow_pickle=True)        
    config=read['configuration'].item()
    particles=read['particles']
    wall=read['wall']
    forces=read['forces'].item()
    
    wc=np.mean(wall[:,1])
    ywd=125
    
    print(wc)
    print(wall[:,1])
    
    ax.plot(wall[:int(len(wall)/2),0],wall[:int(len(wall)/2),1],'-',c=c_wall,lw=2)
    ax.plot(wall[int(len(wall)/2):,0],wall[int(len(wall)/2):,1],'-',c=c_wall,lw=2)
    
    cx=np.mean(particles[:-1,0,:])
    cy=np.mean(particles[:-1,1,:])
    
    dy=100
    
    #ax.set_xlim([0,120])
    #ax.set_ylim([0,60])
    
    ax.set_xlim([30,70])
    ax.set_ylim([20,40])
    
    #ax.set_xlim([config['wall']['xl'], config['wall']['xl']+50])
    #ax.set_ylim([config['wall']['yb'], config['wall']['yt']])
    
    #ax.set_xticks(np.linspace(50,500,5))
    #ax.set_yticks(np.linspace(wc-ywd,wc+ywd,5))
    
    ks=config['membrane']['ks']
   

    for j in range(config['system']['Npart']):
        klx=np.hstack([particles[:-1,0,j], particles[0,0,j]])
        kly=np.hstack([particles[:-1,1,j], particles[0,1,j]])
        plt.plot(klx, kly, '-', c='C3')    
    

        
    

        q_nskip=5
        q_len=10
        
        
        #ax.quiver(    particles[:-1,0,j][::q_nskip]-1.05*q_len*forces['sAf'][:,0,j][::q_nskip],
        #                particles[:-1,1,j][::q_nskip]-1.05*q_len*forces['sAf'][:,1,j][::q_nskip],
        #                q_len*forces['sAf'][:,0,j][::q_nskip], q_len*forces['sAf'][:,1,j][::q_nskip],
        #                angles='xy', scale_units='xy', scale=1, width=0.001, zorder=-98, color="gray", lw=0.1 ) ##Sugando
                        
        #ax.quiver(    particles[:-1,0,j][::q_nskip]-1.05*q_len*forces['mA'][:,0,j][::q_nskip],
        #                particles[:-1,1,j][::q_nskip]-1.05*q_len*forces['mA'][:,1,j][::q_nskip],
        #                q_len*forces['mA'][:,0,j][::q_nskip], q_len*forces['mA'][:,1,j][::q_nskip],
        #                angles='xy', scale_units='xy', scale=1, width=0.001, zorder=-98, color="blue", lw=0.1 ) ##Sugando
                        
        #ax.quiver(    particles[:-1,0,j][::q_nskip],
        #                particles[:-1,1,j][::q_nskip],
        #                10*np.cos(particles[:-1,4,j][::q_nskip]), 10*np.sin(particles[:-1,4,j][::q_nskip]),
        #                angles='xy', scale_units='xy', scale=1, width=0.001, zorder=-98, color="orange", lw=0.05 ) ##Sugando


    fig.savefig("Graph/"+name+'.png',bbox_inches='tight')


    
    return None
    
    
#path_init='/base/gustavo/BioPhy/Trabalhos/2020-04-18-MicropipetteGraph/BaseFile/Save/DC-00100000-00000100-01500000-00001000-00006000-Init.npz'
#path_final='/base/gustavo/BioPhy/Trabalhos/2020-04-18-MicropipetteGraph/BaseFile/Save/DC-00100000-00000100-01500000-00001000-00006000-Final.npz'

#fig=plt.figure(figsize=(14,14/(450./250)))
fig=plt.figure(figsize=(24,24*60/120))
ax=fig.add_axes([0.15,0.15,0.6,0.6])

flist=np.sort(np.array([i for i in os.listdir('Save/')]))
glist=np.sort(np.array([i for i in os.listdir('Graph')]))

aux_f=[i.split(".")[0] for i in flist]
aux_g=[i.split(".")[0] for i in glist]

index=np.array([False if i in aux_g else True for i in aux_f])

flist=flist[index]

#flist=flist[-50:]

for i in flist[::-1]:
    name=i.split('.')[0]
    path='Save/'+i
    plot_SimDescription(path,name)
