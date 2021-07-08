import numpy as np
import matplotlib.pyplot as plt
import os, time

class micropipette_aggregate():
    def __init__(self, comp=False):
        self.path_self=os.path.join(os.path.dirname(__file__))
        if comp:
            os.system("make -C " + self.path_self)
        import micropipette_aggregate.CalcForce as CalcForce
        self.CalcForce = CalcForce
        
        
    def load_ic(self,path=""):
        if path=="":
            self.path_ic=self.path_self+"/InitialCondition/IC.npz"
        else:
            self.path_ic = path
        self.ic = np.load(self.path_ic,allow_pickle=True)
        self.data_ic = {
            'particles':self.ic['particles'],
            'wall':self.ic['wall'],
            'configuration':self.ic['configuration'].item(),
            'forces':self.ic['forces'],
            'instant':self.ic.get('k') and self.ic['k'] or self.ic['instant']
        }
        self.data_cc = {i:self.data_ic[i].copy() for i in self.data_ic}
        
    def print_configuration(self):
        for par in self.data_cc['configuration']:
            print("\n"+par+":")
            for val in self.data_cc['configuration'][par]:
                print("{:>6s}:".format(val),self.data_cc['configuration'][par][val])
            
    def set_configuration(self,config, option=None, val=None):
        if option != None:
            self.data_cc['configuration'][config][option]=val
        else:
            for par in config:
                self.data_cc['configuration'][par[0]][par[1]]=par[2]
        
    def get_cell_mean_xy(self):
        cxy=np.vstack([
            [np.mean(self.data_cc['particles'][:-1,0,i].T) for i in range(self.data_cc['configuration']['system']['Npart'])],
            [np.mean(self.data_cc['particles'][:-1,1,i].T) for i in range(self.data_cc['configuration']['system']['Npart'])]
            ]).T
        return cxy
        
        
    def get_name_string(self):
        config = self.data_cc['configuration']
        strs=["{:08d}".format(int(10000*i)) for i in [
            config['membrane']['ks'],
            config['membrane']['kb'],
            config['membrane']['kA'],
            config['cell']['Fa'],
            config['cell']['rl'],
            config['system']['Af'],
            config['system']['dP']
            ]]
        return "-".join(strs)
        
    def get_Nit_string(self):
        return "{:08d}".format(self.data_cc['configuration']['system']['Nit'])
        
    def get_af_val(self):
        config = self.data_cc['configuration']
        Af=config['system']['Af']
        dP=config['system']['dP']
        particles=self.data_cc['particles']
        center=np.vstack([
            [np.mean(particles[:-1,0,i].T) for i in range(config['system']['Npart'])],
            [np.mean(particles[:-1,1,i].T) for i in range(config['system']['Npart'])]
            ]).T
        nn=len(center)
        nd=np.zeros((nn,nn))
        for ii in range(nn-1):
            for jj in range(ii+1,nn):
                nd[ii,jj]=1/((center[ii,0]-center[jj,0])*(center[ii,0]-center[jj,0])+ \
                            (center[ii,1]-center[jj,1])*(center[ii,1]-center[jj,1]))
        nd+=nd.T
        
        rho=np.sum(nd,axis=1)
        
        var=np.zeros(len(center))
        
        ## outside micropipette
        
        index=   (center[:,0] > config['wall']['xl']-2)\
                *(center[:,0] < config['wall']['xr'])\
                *(center[:,1] > config['wall']['yb']-2)\
                *(center[:,1] < config['wall']['yt']+2)
        
        mrho=np.mean((rho*(center[:,0] < config['wall']['xl']))[~index])
        var[rho < mrho]=Af
        
        index=   (center[:,0] > config['wall']['xl']-6)\
                *(center[:,0] < config['wall']['xl']+3)\
                *(center[:,1] > config['wall']['yt'])\
                *(center[:,1] < config['wall']['yt']+8)
        
        var[index]=Af
        
        index=   (center[:,0] > config['wall']['xl']-6)\
                *(center[:,0] < config['wall']['xl']+3)\
                *(center[:,1] < config['wall']['yb'])\
                *(center[:,1] > config['wall']['yt']-8)
                
        var[index]=Af
        ## inside micropipette
        
        index=   (center[:,0] > config['wall']['xl']-10)\
        *(center[:,0] < config['wall']['xr'])\
        *(center[:,1] > config['wall']['yb']-1.5)\
        *(center[:,1] < config['wall']['yt']+1.5)
        
        mrho=np.mean(rho[index])
        var[index*(rho < mrho) ]=Af-dP
        
        return var
            
    def integrate(self,num_cores):
        config = self.data_cc['configuration'].copy()
        af_val=self.get_af_val()
        out=self.CalcForce.integrate(self.data_cc['particles'], self.data_cc['wall'], config['membrane']['Nm'], config['system']['Npart'],
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
                    config['active']['tau'], config['system']['mu0'],  af_val,
                    config['system']['dt'], config['system']['Ndt'], config['system']['Nret'], num_cores)
                    
        state=out[-1]
        particles=out[0][:,:,:,-1]
        wall=out[1][:,:,-1]
        
        f_mb, f_mA, f_ms, f_mF, f_ns, f_nF, f_NNs, \
        f_NNF, f_cF, f_sAf, f_wF = out[2:-1]

        forces={'mb':f_mb[:,:,:,-1], 'mA':f_mA[:,:,:,-1], 'ms':f_ms[:,:,:,-1], 
            'mF':f_mF[:,:,:,-1], 'ns':f_ns[:,:,:,-1], 'nF':f_nF[:,:,:,-1],
            'NNs':f_NNs[:,:,-1], 'NNF':f_NNF[:,:,-1], 'cF':f_cF[:,:,:,-1],
            'sAf':f_sAf[:,:,:,-1], 'swF':f_wF[:,:,:,-1]}
            
        
        config['system']['tv']+=config['system']['dt']*config['system']['Ndt']
        config['system']['Nit']+=config['system']['Ndt']
        
        instant=self.data_cc['instant']+1
        
        data = {
            'particles':particles,
            'wall':wall,
            'configuration':config,
            'forces':forces,
            'instant':instant,
        }
        
        return data, state

    def evolve(self,num_cores=1):
        state=True
        while state:
            data,state=self.integrate(num_cores)
        self.data_cc=data
        
    def save_data(self,path):
        np.savez_compressed(path,
            particles=self.data_cc['particles'],
            wall=self.data_cc['wall'],
            forces=self.data_cc['forces'],
            configuration=self.data_cc['configuration'],
            instant=self.data_cc['instant']
        )
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
