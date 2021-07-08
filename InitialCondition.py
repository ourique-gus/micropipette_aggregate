import numpy as np
import time
import CalcForce as CalcForce
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import os


###########################


def plot_data(name, particles, wall, config):

	fig=plt.figure(figsize=(8*config['system']['Lx']/config['system']['Ly'],8))
	
	ax=fig.add_axes([0.1,0.1,0.6,0.6])
	
	
	ax.xaxis.set_major_locator(ticker.MaxNLocator(4))
	ax.yaxis.set_major_locator(ticker.MaxNLocator(4))

	x_minorLocator = ticker.AutoMinorLocator(4)
	y_minorLocator = ticker.AutoMinorLocator(4)
	ax.xaxis.set_minor_locator(x_minorLocator)
	ax.yaxis.set_minor_locator(y_minorLocator)

	ax.tick_params(axis='both', which='major', direction='in', top = True, right=True)
	ax.tick_params(axis='both', which='minor', direction='in', top=True,right=True)
	ax.set_ylabel('y')
	ax.set_xlabel('x')
	
	
	dlim_var=0*10
	
	c_cell='C3'
	c_wall='C0'
	c_pres='black'

	
	#ax.set_xlim([0,config['system']['Lx']])
	#ax.set_ylim([0,config['system']['Ly']])
	ax.set_xlim([0,config['system']['Lx']])
	ax.set_ylim([0,config['system']['Ly']])
	
	#ax.set_xlim([40,100])
	#ax.set_ylim([45,75])
	
	ax.plot(wall[:int(len(wall)/2),0],wall[:int(len(wall)/2),1],'-',c=c_wall,lw=2)
	ax.plot(wall[int(len(wall)/2):,0],wall[int(len(wall)/2):,1],'-',c=c_wall,lw=2)
	

	for j in range(config['system']['Npart']):
	
		#ax.set_xlim([0,config['system']['Lx']])
		#ax.set_ylim([0,config['system']['Ly']])
		
	
		
		klx=np.hstack([particles[:-1,0,j], particles[0,0,j]])
		kly=np.hstack([particles[:-1,1,j], particles[0,1,j]])

		
		ax.plot(klx,kly,'-', c=c_cell, lw=1, zorder=-96)
		
		

	fig.savefig("Graph/"+name+'.pdf',bbox_inches='tight')
	plt.close(fig)

###########################

#ax=fig.add_axes([0.15,0.15,0.5,0.5])

np.set_printoptions(precision=4, formatter={'all': lambda x: '{:10.6f}'.format(x)}, linewidth=240)
print()
print('Starting...')
print()


config={'membrane':{}, 'nucleus':{}, 'cell':{}, 'active':{}, 'system':{}, 'wall':{}}

#Area=2000

config['membrane']['Nm']=50
config['membrane']['kb']=0.01
config['membrane']['kA']=10
config['membrane']['ks']=100
#config['membrane']['r0']=np.sqrt(Area/np.pi)*2*np.sin(np.pi/config['membrane']['Nm'])
config['membrane']['d0']=1
config['membrane']['r0']=2*np.pi*config['membrane']['d0']/config['membrane']['Nm']
config['membrane']['Fr']=1
config['membrane']['re']=config['membrane']['r0']
config['membrane']['dT']=2*np.pi/config['membrane']['Nm']
#config['membrane']['d0']=config['membrane']['r0']/(2*np.sin(config['membrane']['dT']/2))
config['membrane']['A0']=np.pi*config['membrane']['d0']**2


config['nucleus']['ks']=0.00000
config['nucleus']['r0']=config['membrane']['d0']
config['nucleus']['re']=config['nucleus']['r0']/(10**(1./3))
config['nucleus']['Fr']=0
config['nucleus']['ms']=1
config['nucleus']['dp']=0.1


config['cell']['Fr']=100
config['cell']['Fa']=0.1
config['cell']['re']=1.25*config['membrane']['r0']
config['cell']['rl']=2*config['membrane']['r0']
print(config['membrane']['r0'])


config['active']['v0']=0
config['active']['Dt']=1.0
config['active']['Dr']=0.01
config['active']['tau']=1

config['system']['dt']=0.001
config['system']['Lx']=200
config['system']['Ly']=60
config['system']['dx']=0.3
config['system']['dy']=0.3
config['system']['Ncalc']=1
config['system']['sqdt']=np.sqrt(config['system']['dt'])
config['system']['sq3']=np.sqrt(3)
config['system']['Acte']=np.sqrt(2*config['active']['Dr']*config['system']['dt'])
config['system']['mu0']=1
config['system']['Ndt']=10000
config['system']['Nret']=10
config['system']['Af']=0.0
config['system']['tv']=0
config['system']['Nit']=0

#######################################
########### Criando paredes ###########

chn_width=10
chn_length=150
node_spc=config['membrane']['r0']/2.

node_x=[0, -chn_length, -chn_length, 0]
node_y=[chn_width/2, chn_width/2, 0, 0]

xr=np.hstack([np.linspace(node_x[i], node_x[i+1],100000) for i in range(len(node_x)-1)])
yr=np.hstack([np.linspace(node_y[i], node_y[i+1],100000) for i in range(len(node_y)-1)])

wx=[]
wy=[]

wx.append(xr[0])
wy.append(yr[0])

pos=0

while pos < len(xr)-1:
	while (np.sqrt((wx[-1]-xr[pos])**2+(wy[-1]-yr[pos])**2) < node_spc) & (pos < len(xr)-1):
		#print( wx[-1],xr[pos],wy[-1],yr[pos], ((wx[-1]-xr[pos])**2+(wy[-1]-yr[pos])**2) , node_spc*node_spc)
		pos+=1
	wx.append(xr[pos])
	wy.append(yr[pos])
	
wx=np.array(wx)
wy=np.array(wy)

n_smooth=100
f_smooth=0.5

for ii in range(n_smooth):
	wx[1:-1]=(1-f_smooth)*wx[1:-1]+0.5*f_smooth*(wx[:-2]+wx[2:])
	wy[1:-1]=(1-f_smooth)*wy[1:-1]+0.5*f_smooth*(wy[:-2]+wy[2:])
	
wx=np.hstack([wx,wx])+config['system']['Lx']
wy=np.hstack([wy+chn_width/2,-wy-chn_width/2])+config['system']['Ly']/2

	
	
wall=np.vstack([wx,wy, 0*wx, 0*wy]).T

#plt.plot(wall[:,0],wall[:,1])
#plt.show()
#exit()

print (np.min(wall[:,0]), np.min(wall[0,1]), np.min(wall[-1,1]), np.min(wall[int(len(wall)/2)-1,1]))

config['wall']['Nw']=len(wall)
config['wall']['Fr']=1000
config['wall']['re']=1.75*config['membrane']['r0']
#wall=np.array([[-9999,-9998],[-9999,-9998],[0,0],[0,0]]).T


config['wall']['xl']=np.min(wall[:,0])
config['wall']['xr']=config['system']['Lx']
config['wall']['yb']=np.mean(wall[int(len(wall)/2):,1])
config['wall']['yt']=np.mean(wall[:int(len(wall)/2),1])

print(config['wall']['xl'], config['wall']['xr'], config['wall']['yb'], config['wall']['yt'])

print('Cell radius: ~{:7.3f}'.format(config['membrane']['d0']))
print('Best box width: {:7.3}'.format(np.max([config['wall']['re'], config['cell']['re'], config['membrane']['re'], config['cell']['rl']])))
print(config['wall']['re'], config['cell']['re'], config['membrane']['re'], config['cell']['rl'])


#######################################
########### Criando células ###########

config['system']['Npart']=500
Ncells=0
csm=0.7
cells=[]
agg_iradius=21

xc_v=[]
yc_v=[]

while len(cells) < config['system']['Npart']:
    xc=np.random.uniform(-agg_iradius,agg_iradius)
    yc=np.random.uniform(-agg_iradius,agg_iradius)
    while np.sqrt(xc*xc+yc*yc) > agg_iradius:
        xc=np.random.uniform(-agg_iradius,agg_iradius)
        yc=np.random.uniform(-agg_iradius,agg_iradius)
    dist=[np.sqrt((xc-xc_v[i])*(xc-xc_v[i])+(yc-yc_v[i])*(yc-yc_v[i])) for i in range(len(cells))]
    if (len(dist) and min(dist) or 2*config['membrane']['d0']) > csm*2*config['membrane']['d0']:
        print(len(cells))
        xc_v.append(xc)
        yc_v.append(yc)
        
        cell_x=np.hstack([csm*config['membrane']['d0']*np.cos((i)/config['membrane']['Nm']*2*np.pi)+xc for i in range(config['membrane']['Nm'])]+ [xc])
        cell_y=np.hstack([csm*config['membrane']['d0']*np.sin((i)/config['membrane']['Nm']*2*np.pi)+yc for i in range(config['membrane']['Nm'])]+ [yc])
	        
        rAngle=0*np.random.uniform(0,2*np.pi, config['membrane']['Nm']+1)
	
        cell=np.vstack([cell_x, cell_y, 
	        config['active']['v0']*np.cos(rAngle),
	        config['active']['v0']*np.sin(rAngle),
	        0*np.random.uniform(0,2*np.pi, config['membrane']['Nm']+1)]).T
        cells.append(cell)
	
particles=np.dstack(cells)

particles[:,1,:]+=-np.mean(particles[:,1,:])+config['system']['Ly']/2
particles[:,0,:]+=-np.max(particles[:,0,:])+np.min(wall[:,0])-0.2

plot_data('temp', particles, wall, config)


#np.savez('BaseParticles',particles=particles, wall=wall, configuration=config, k=0)
#config['system']['Npart']=500
#particles=np.load('BaseParticles.npz')['particles']


#plot_data('F00000000', particles, wall, 0, config)


s=time.time()
for i in np.linspace(0.1,1,10):
    out=out=CalcForce.integrate(particles, wall, config['membrane']['Nm'], config['system']['Npart'],
							    config['wall']['Nw'], config['system']['Lx'], config['system']['Ly'], config['system']['dx'],
							    config['system']['dy'],

							    config['membrane']['kb'], i*config['membrane']['kA'], config['membrane']['ks'],
							    config['membrane']['A0'], config['membrane']['r0'],
							    config['membrane']['Fr'], config['membrane']['re'],

							    config['nucleus']['ks'], config['nucleus']['r0'],
							    config['nucleus']['Fr'], config['nucleus']['re'],
							    config['nucleus']['ms'], config['nucleus']['dp'],

							    config['cell']['Fr'], config['cell']['Fa'],
							    config['cell']['re'], config['cell']['rl'],

							    config['wall']['Fr'], config['wall']['re'], 
							    config['wall']['xl'], config['wall']['xr'],
							    config['wall']['yb'], config['wall']['yt'],

							    config['active']['v0'],config['active']['Dt'], config['active']['Dr'],
							    config['active']['tau'], config['system']['mu0'],  np.zeros(config['system']['Npart'])*config['system']['Af'],
							    config['system']['dt'], 5000, 1)

    print('time',time.time()-s)

    particles=out[0][:,:,:,-1]

    f_mb, f_mA, f_ms, f_mF, f_ns, f_nF, f_NNs, \
    f_NNF, f_cF, f_sAf, f_wF = out[2:-1]


    forces={'mb':f_mb[:,:,:,-1], 'mA':f_mA[:,:,:,-1], 'ms':f_ms[:,:,:,-1], 
	    'mF':f_mF[:,:,:,-1], 'ns':f_ns[:,:,:,-1], 'nF':f_nF[:,:,:,-1],
	    'NNs':f_NNs[:,:,-1], 'NNF':f_NNF[:,:,-1], 'cF':f_cF[:,:,:,-1],
	    'sAf':f_sAf[:,:,:,-1], 'swF':f_wF[:,:,:,-1]}
	
"""
dx=0.001
while np.sum(np.abs(forces['swF'])) < 1E-5:
    particles[:,0,:]+=dx
    print(np.max(particles[:,0,:]))
    
    out=CalcForce.integrate(particles, wall, config['membrane']['Nm'], config['system']['Npart'],
				config['wall']['Nw'], config['system']['Lx'], config['system']['Ly'],
				config['system']['dx'], config['system']['dy'],

				config['membrane']['kb'], config['membrane']['kA'], config['membrane']['ks'],
				config['membrane']['A0'], config['membrane']['r0'],
				config['membrane']['Fr'], config['membrane']['re'],

				config['nucleus']['ks'], config['nucleus']['r0'],
				config['nucleus']['Fr'], config['nucleus']['re'],
				config['nucleus']['ms'], config['nucleus']['dp'],

				config['cell']['Fr'], config['cell']['Fa'],
				config['cell']['re'], config['cell']['rl'],

				config['wall']['Fr'], config['wall']['re'], 
				config['wall']['xl'], config['wall']['xr'],
				config['wall']['yb'], config['wall']['yt'],

				config['active']['v0'],config['active']['Dt'], config['active']['Dr'],
				config['active']['tau'], config['system']['mu0'], 0*config['system']['Af'],
				config['system']['dt'], 1, 1)


    f_mb, f_mA, f_ms, f_mF, f_ns, f_nF, f_NNs, \
    f_NNF, f_cF, f_sAf, f_wF, f_ms_up, f_ms_dw = out[2:-1]


    forces={'mb':f_mb[:,:,:,-1], 'mA':f_mA[:,:,:,-1], 'ms':f_ms[:,:,:,-1], 
	    'mF':f_mF[:,:,:,-1], 'ns':f_ns[:,:,:,-1], 'nF':f_nF[:,:,:,-1],
	    'NNs':f_NNs[:,:,-1], 'NNF':f_NNF[:,:,-1], 'cF':f_cF[:,:,:,-1],
	    'sAf':f_sAf[:,:,:,-1], 'swF':f_wF[:,:,:,-1],
	    'ms_up':f_ms_up[:,:,:,-1], 'ms_dw':f_ms_dw[:,:,:,-1]}
    

"""
## Ajustando configurações 

## Iniciando a figura

for i in os.listdir('InitialCondition'):
	os.remove('InitialCondition/'+i)
	
for i in os.listdir('Save'):
	os.remove('Save/'+i)
	
plot_data('var', particles, wall, config)
					
np.savez('InitialCondition/IC',particles=particles, wall=wall, forces=forces , configuration=config, k=0)
				
		


