import os
import numpy as np
import sys
from dataclasses import dataclass

@dataclass
class Molecule():
    name: str
    weight: float
    magic_a: float
    magic_b: float
    mole_fraction: float
    
class mrkHolloway():
    """_summary_
    """
    def __init__(self):
        self.output_file = input('Enter name of output file \n')
        starting_vals = input('Enter T (deg C) and Molar Volume (cc/mole)\n')
        temp, volume = starting_vals.split()
        self.tstart = float(temp)
        if self.tstart < 0.01:
            sys.exit()
        self.vstart = float(volume)
        self.molecules = [
            Molecule(name='CO2', weight=44.0, magic_a=46.e6, magic_b=2.97e1, mole_fraction=input(f'Enter mole fraction CO2\n')),
            Molecule(name='CO', weight=28.0, magic_a=16.98e6, magic_b=2.738e1, mole_fraction=input('Enter mole fraction CO \n')),
            Molecule(name='CH4', weight=16.0, magic_a=31.59e6, magic_b=2.9703e1, mole_fraction=input('Enter mole fraction CH4 \n')),
            Molecule(name='H2', weight=2.0, magic_a=3.56e6, magic_b=1.515e1, mole_fraction=input('Enter mole fraction H2\n')),
            Molecule(name='H2O', weight=18.0, magic_a=35.e6, magic_b=1.46e1, mole_fraction=input('Enter mole fraction H20\n')),
            Molecule(name='H2S', weight=34.0, magic_a=87.9e6, magic_b=2.0e1, mole_fraction=input('Enter mole fraction H2S \n')),
            Molecule(name='SO2', weight=64.06, magic_a=142.6e6, magic_b=3.94e1, mole_fraction=input('Enter mole fraction SO2 \n')),
            Molecule(name='N2', weight=28.0, magic_a=15.382e6, magic_b=2.68e1, mole_fraction=input('Enter mole fraction N2 \n'))]
        
        # for molecule in self.molecules:
        #     molecule.mole_fraction = float(input(f'Enter mole fraction {molecule}\n'))
        # self.names = ['CO2', 'CO', 'CH4', 'H2', 'H2O', 'H2S', 'SO2', 'N2']
        # self.formwt = [44.0, 28.0, 16.0, 2.0, 18.0, 34.0, 64.06, 28.0]
        #self.mixnum = 8
        
            
     
    
    def check_totals(self):
        if not sum(Molecule.mole_fraction) == 1.0:
            print('The mole fractions do not add up to 1.0! Start again')
            self.obtain_input()
        else:
            pass
        
        
    def calc_form_weight(self):
        self.form = 0.0
        for i in range(0,self.mixnum):
            self.form += self.input_mixture[i] * self.formwt[i]
            
            
    def calc_isochore_volume(self):
        self.voli = [self.vstart + 10*i - 10 for i in range(1,6,1)]
        self.density = [self.form / x for x in self.voli]
            
    
    def temp_pressure_calc(self, tk):
        rbar = 83.117
        for i in range(1,12, 1):
            t = tk - 1e2 + 1e2 * i
            tc = t - 273.15
            #bsum = bterm for each gas * mole fraction
            asum, bsum = self.mrkmix(tk)
            pout = []
            for j in range(5):
                
                vol = self.voli[j]
                p = 0.0
                if not (vol < bsum - 1.0):
                    # aterm = a_molecule
                    # bterm = b_moledcule
                    aterm = ((rbar * t) / (vol - bsum))
                    bterm = asum / ((np.sqrt(t)) * ((vol**2) + (bsum * vol)))
                    p = aterm - bterm
                pout.append(p)
                
            # Print temp/pressure outputs
            print(f'{tc:10.0f} {"":<10s} ' + 
                  ''.join(f'{x:10.0f}' for x in pout))
            
            
    def additional_runs(self):
        rerun = input('Do you want to do another composition? No=0, Yes=1 \n')
        if rerun == '1':
            self.run()
        else:
            sys.exit(0)

    def mrkmix(self, tk):
        # finish translating mrkmix
        """_summary_

        Args:
            tk (_type_): _description_

        Returns:
            _type_: _description_
        """
        a = [46.e6, 16.98e6, 31.59e6, 3.56e6, 35.e6, 87.9e6, 142.6e6, 15.382e6]
        b = [2.97e1, 2.738e1, 2.9703e1, 1.515e1, 1.46e1, 2.0e1, 3.94e1, 2.68e1]
        r = 82.05
        if (tk < 1e-4):
            t = 1.0
        else:
            t = tk
            
        tcelsius = t - 273.15
        r2t = r**2 * t**2.5
        rt = r*t**1.5
        ah2om = 166.8 - .19308 * tcelsius + .1864e-3 * tcelsius**2 - .71288e-7 * tcelsius**3
        if (tcelsius < 6e2):
            ah2om = 4.221e3 - 3.1227e1 * tcelsius + 8.7485e-2 * tcelsius**2 - 1.07295e-4 \
                * tcelsius**3 + 4.86111e-8 * tcelsius**4
        elif (tcelsius > 1200):
            ah2om = 140.0 - 0.05 * tcelsius
            
        ah2om *= 10e5
        aco2m = 73.03 - 0.0714 * tcelsius + 2.157e-5 * tcelsius**2
        aco2m *= 10e5
        xk = np.exp(-11.071 + (5953./t) - (2.746e6/(t**2)) + (4.646e8/t**3))
        co2h2o = xk * 0.5 * r2t
        co2h2o += np.sqrt(a[0] * a[4])
        asum = 0.0
        bsum = 0.0
        for i in range(self.mixnum):
            bsum += b[i] * self.input_mixture[i]
            for j in range(self.mixnum):
                if (i == j):
                    if ((i != 4) and (j != 0)):
                        asum += self.input_mixture[i] * self.input_mixture[j] * a[i]
                    elif (i == 4):
                        asum += self.input_mixture[i] * self.input_mixture[j] * ah2om
                    elif (i == 0):
                        asum += self.input_mixture[i] * self.input_mixture[j] * aco2m
                elif (( i == 4) and (j == 0) or (i == 0) and (j==4)):
                    asum += self.input_mixture[i] * self.input_mixture[j] * co2h2o
                else:
                    asum += self.input_mixture[i] * self.input_mixture[j] * np.sqrt(a[i]*a[j])
        asum /= 1.013
        return asum, bsum
    
    #### ADD A DATACLASS FOR THIS  ####
    ### number of temp/pressure could be variablen###
    ### how many increments of what temp ####
    ### max temp to go to ###
    ### add plots ###
    ### add write to output ###
    ### at end ask if new composition?
    ### if yes ask if new file
    ### if same file re-enter starting info
    # if bad composition save name and just start from fractions again
    
    
    
    
    
    def run(self):
       
        self.check_totals()
        self.calc_form_weight()
        print('The mole fractions are:')
        formstr = ' '.join(f"{x:<9}" for x in self.names)
        mix_str = ' '.join(f"{x:<9.3f}" for x in self.input_mixture)
        
        # Print input mixture
        print(formstr)
        print(mix_str)
        print(f'The starting temperature = {self.tstart:3f} deg. c')
        print(f'The starting molar volume = {self.vstart:3f} cm^3/mol')
        
        # Start Calculation Isochore
        tk = self.tstart + 273.15 # convert to kelvin
        self.calc_isochore_volume()
        # Print volume and density
        print(f'{"":<30s} Isochore pressure in bars')
        print(f'{"T deg. c mol vol =":<30s}' + 
              ' '.join(f"{x:<10.3f}" for x in self.voli))
        
        print(f'{"density =":<30s}' + 
              ' '.join(f'{x:<10.4f}' for x in self.density))
        
        # Calc temp pressure
        self.temp_pressure_calc(tk)
        self.additional_runs()
        
                                  

if __name__ == '__main__':
    go = mrkHolloway()
    go.run()
    
    
    
        
        
        