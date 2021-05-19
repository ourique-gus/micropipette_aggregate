import micropipette_aggregate.CalcForce as CalcForce
import micropipette_aggregate.CalcForceOMP as CalcForceOMP
import numpy as np
import matplotlib.pyplot as plt
import os

class micropipette_aggregate():
    def __init__(self):
        self.CalcForce = CalcForce
        self.CalcForceOMP = CalcForceOMP
        self.path_self=os.path.join(os.path.dirname(__file__))
        
    def load_ic(self,path=""):
        if path=="":
            self.path_ic=self.path_self+"/InitialCondition/IC.npz"
        else:
            self.path_ic = path
        self.ic = np.load(self.path_ic,allow_pickle=True)
        self.data_ic = {
            'particles':self.ic['particles'],
            'walls':self.ic['wall'],
            'configuration':self.ic['configuration'].item(),
            'forces':self.ic['forces'],
            'instant':self.ic['k']
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
        
            
    def integrate(self,omp=False):
        func = omp and self.CalcForceOMP or self.CalcForce
        config = self.data_cc['configuration'].copy()
        out=func.integrate(self.data_cc['particles'], self.data_cc['walls'], config['membrane']['Nm'], config['system']['Npart'],
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
            'walls':wall,
            'configuration':config,
            'forces':forces,
            'instant':instant,
        }
        
        return data, state

    def evolve(self,omp=False):
        state=True
        while state:
            data,state=self.integrate(omp)
        self.data_cc=data
        
            
