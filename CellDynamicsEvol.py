import numpy as np
import time
import CalcForce as CalcForce
import os, sys
from scipy.spatial import Voronoi
from scipy.spatial import ConvexHull
from matplotlib.patches import Polygon

def GetHull(particles,config):
    nucleus=particles[-1,0:2,:].T
    idvar=np.arange(len(nucleus)).astype(int)
    index=[]
    index.append(idvar[nucleus[:,0] > config['wall']['xl']])
    #mx, my = np.mean(nucleus[:,0]), np.mean(nucleus[:,1])
    
    for i in range(3):
    
        hull=ConvexHull(nucleus)
        index.append(idvar[hull.vertices])
        aux=np.ones(len(idvar)).astype(bool)
        aux[hull.vertices]=False
        idvar=idvar[aux]
        nucleus=nucleus[aux]

    index=np.hstack(index)
    #pol=Polygon(var)
    #index=pol.contains_points(particles[-1,0:2,:].T)
    
    var=np.ones(len(particles[-1,0,:])).astype(bool)
    var[index]=False
    
    index=(particles[-1,0,:] > config['wall']['xl'])
    var[index]=False
    
    index=   (particles[-1,0,:] > config['wall']['xl'])\
            *(particles[-1,0,:] < config['wall']['xr'])\
            *(particles[-1,1,:] > config['wall']['yb'])\
            *(particles[-1,1,:] < config['wall']['yt'])
            
    var[index]=True

    return var
    
def GetDensity(particles,config):
    #nucleus=particles[-1,0:2,:].T
    nucleus=np.vstack([
        [np.mean(particles[:-1,0,i].T) for i in range(config['system']['Npart'])],
        [np.mean(particles[:-1,1,i].T) for i in range(config['system']['Npart'])]
        ]).T
    nn=len(nucleus)
    nd=np.zeros((nn,nn))
    for ii in range(nn-1):
        for jj in range(ii+1,nn):
            nd[ii,jj]=1/((nucleus[ii,0]-nucleus[jj,0])*(nucleus[ii,0]-nucleus[jj,0])+ \
                        (nucleus[ii,1]-nucleus[jj,1])*(nucleus[ii,1]-nucleus[jj,1]))
    nd+=nd.T
    
    rho=np.sum(nd,axis=1)
    
    index=   (nucleus[:,0] > config['wall']['xl']-2)\
            *(nucleus[:,0] < config['wall']['xr'])\
            *(nucleus[:,1] > config['wall']['yb']-2)\
            *(nucleus[:,1] < config['wall']['yt']+2)
    
    mrho=np.mean(rho[~index])
    var=rho < 1*mrho 
    var[index]=False
    
    return var
    

from optparse import OptionParser

parser = OptionParser()
parser.add_option("-f", "--file", dest="filename", default="IC.npz")
parser.add_option("-s", "--ks", dest="ks", type='float', default=10)
parser.add_option("-A", "--kA", dest="kA", type='float', default=1)
parser.add_option("-P", "--dP", dest="dP", type='float', default=1)
parser.add_option("-a", "--fa", dest="fa", type='float', default=0.1)
parser.add_option("-b", "--kb", dest="kb", type='float', default=0.01)

#parser.add_option("-s", "--ks", dest="ks", type='float', default=100)
#parser.add_option("-A", "--kA", dest="kA", type='float', default=100)
#parser.add_option("-P", "--dP", dest="dP", type='float', default=1)
#parser.add_option("-a", "--fa", dest="fa", type='float', default=0.1)
#parser.add_option("-b", "--kb", dest="kb", type='float', default=0.01)

(options, args) = parser.parse_args()
options=vars(options)

read=np.load(options['filename'],allow_pickle=True)
ks=options['ks']
kA=options['kA']
dP=options['dP']
fa=options['fa']
kb=options['kb']

particles=read['particles']


wall=read['wall']
if wall.shape[1] < 3:
    wall=np.vstack([wall.T,0*wall.T]).T
config=read['configuration'].item()

max_NaN=10000
num_tt=5000
num_tvar=1

LpZero=57.

#particles[:,0,0]-=35


kNaN=0
kevolpos=read['k']

#config['membrane']['kb']=0.1
#config['membrane']['kA']=100
#config['membrane']['ks']=100
#config['membrane']['Fr']=10
#config['nucleus']['ks']=0.001
#config['nucleus']['Fr']=10
#config['nucleus']['ms']=1

config['membrane']['kb']=kb
config['membrane']['kA']=kA
config['membrane']['ks']=ks
config['cell']['Fa']=fa
Af0=dP

config['system']['dt']=0.001
config['system']['Ndt']=10000
config['system']['Nret']=1

config['wall']['xl']=config['wall']['xl']-5
config['wall']['yb']=config['wall']['yb']
config['wall']['yt']=config['wall']['yt']
config['wall']['re']=1.75*config['membrane']['r0']

config['system']['dx']=0.30
config['system']['dy']=0.30
config['cell']['rl']=2*config['membrane']['r0']

config['wall']['Fr']=500

#Af0=0.00500
#Af0=0.00250
#Af0=0.00100
#Af0=0.00050


mkb=int(config['membrane']['kb']*1E4)
mkA=int(config['membrane']['kA']*1E4)
mks=int(config['membrane']['ks']*1E4)
nks=int(config['nucleus']['ks']*1E4)
cfa=int(config['cell']['Fa']*1E4)
sAf=int(Af0*1E4)

base_name='D-{:08d}-{:08d}-{:08d}-{:08d}-{:08d}-{:08d}'.format(mkb,mkA,mks,nks,sAf,cfa)

#   if not os.path.exists("Save/"+base_name):
#    os.makedirs("Save/"+base_name)



#par_dict[name]={'particles':particles, 'wall':wall , 'configuration':config}


#var=config['system']['Af']
#config['system']['Af']=0

######################################################################
########################## Condição Inicial ##########################
######################################################################

name_Af='{:08d}'.format(int(Af0*1E4))
name_it='{:08d}'.format(int(kevolpos))

#config['system']['Ndt'], config['system']['Nret']

######################################################################
######################################################################
######################################################################

for tt in range(num_tt):
	s=time.time()
	for tvar in range(num_tvar):

		kevolpos+=1
		name_it='{:08d}'.format(int(kevolpos))
		
		base=np.zeros(config['system']['Npart'])
		base[GetDensity(particles,config)]=Af0
        
		config['system']['Af']=base
		
		pv=np.dstack([particles[:,:,i].T for i in range(config['system']['Npart']) ])
		wv=wall.T
		#print(pv.shape, particles.shape)
		out=CalcForce.integrate(particles, wall, config['membrane']['Nm'], config['system']['Npart'],
							config['wall']['Nw'], config['system']['Lx'], config['system']['Ly'], config['system']['dx'],
							config['system']['dy'],

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
							config['active']['tau'], config['system']['mu0'],  config['system']['Af'],
							config['system']['dt'], config['system']['Ndt'], config['system']['Nret'])
		nanvalue=out[-1]
		while nanvalue==True:
			print ('Warning: NaN',kNaN)
			kNaN+=1
			out=CalcForce.integrate(particles, wall, config['membrane']['Nm'], config['system']['Npart'],
							config['wall']['Nw'], config['system']['Lx'], config['system']['Ly'], config['system']['dx'],
							config['system']['dy'],

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
							config['active']['tau'], config['system']['mu0'],  config['system']['Af'],
							config['system']['dt'], config['system']['Ndt'], config['system']['Nret'])
			nanvalue=out[-1]

		particles=out[0][:,:,:,-1]
		wall=out[1][:,:,-1]

		
		f_mb, f_mA, f_ms, f_mF, f_ns, f_nF, f_NNs, \
		f_NNF, f_cF, f_sAf, f_wF = out[2:-1]

		forces={'mb':f_mb[:,:,:,-1], 'mA':f_mA[:,:,:,-1], 'ms':f_ms[:,:,:,-1], 
			'mF':f_mF[:,:,:,-1], 'ns':f_ns[:,:,:,-1], 'nF':f_nF[:,:,:,-1],
			'NNs':f_NNs[:,:,-1], 'NNF':f_NNF[:,:,-1], 'cF':f_cF[:,:,:,-1],
			'sAf':f_sAf[:,:,:,-1], 'swF':f_wF[:,:,:,-1]}

	
	if kevolpos/1.==int(kevolpos/1):
	    #np.savez_compressed("Save/"+base_name+'/S'+name_it,particles=particles, wall=wall, forces=forces , configuration=config, k=kevolpos)
	    np.savez_compressed("Save/"+'S'+name_it,particles=particles, wall=wall, forces=forces , configuration=config, k=kevolpos)
	    print(kevolpos,time.time()-s)
	
	#for ii in range(config['system']['Npart']):
	#    if particles[-1,0,ii] > 95:
	#        particles[:,0,ii]+=-particles[-1,0,ii]+np.min(particles[-1,0,:])-5.0*config['membrane']['d0']
	#        particles[:,1,ii]+=-particles[-1,1,ii]+config['system']['Ly']/2+ np.random.uniform(-13,13)
	
