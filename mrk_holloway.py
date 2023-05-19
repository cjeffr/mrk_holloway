import os
import numpy as np
import sys
import pandas as pd
from dataclasses import dataclass

@dataclass
class Molecule():
    name: str
    weight: float
    mrk_molecular_interaction: float
    mrk_molecular_volume: float
    mole_fraction: float
    
class mrkHolloway():
    """_summary_.
    """
    def __init__(self, output_file, data):
        self.output_file = output_file
        self.tstart = data['tstart']
        self.vstart = data['vstart']
        self.t_increment = data['t_increment']
        self.max_t = data['max_t']
        self.molecules = data['molecules']
        self.temps_k = self.calculate_temp_increments()
        self.mixnum = len(self.molecules)
        self.formula_weight = sum(gas.weight * gas.mole_fraction for gas in self.molecules)
        
    
    def check_totals(self):
        if not sum(gas.mole_fraction for gas in self.molecules) == 1.0:
            print('The mole fractions do not add up to 1.0! Start again')
            self.obtain_input()
        else:
            pass
        
    def calculate_temp_increments(self):
        temp_list = []
        temp_list.append(self.tstart)
        nearest_temp = ((self.tstart+ (self.t_increment -1)) // self.t_increment) * self.t_increment
        next_t = nearest_temp
        temp_list.append(next_t)
        while next_t < self.max_t:
            next_t += self.t_increment
            temp_list.append(next_t)
        temp_list = [x + 273.15 for x in temp_list]
        return temp_list
              
            
    def calc_isochore_volume(self):
        self.density = self.formula_weight / self.vstart
            
    
    def temp_pressure_calc(self, tk):
        rbar = 83.117
        self.pressure_output = []
        for t in self.temps_k:
            tc = t - 273.15
            asum, bsum = self.mrkmix(tk)
            
            vol = self.vstart
            p = 0.0
            if not (vol < bsum - 1.0):
                aterm = ((rbar * t) / (vol - bsum))
                bterm = asum / ((np.sqrt(t)) * ((vol**2) + (bsum * vol)))
                p = aterm - bterm
            pout = p
                
            self.pressure_output.append(pout)
        # Print temp/pressure outputs
            print(f'{tc:10.0f} {"":<10s} ' + 
                    ''.join(f'{pout:10.0f}'))
            
            
    def additional_runs(self):
        rerun = input('Do you want to do another composition? (y/n) \n')
        if rerun.lower() == 'y':
            print(f'does {self.output_file} exist?', os.path.exists(self.output_file))
            while os.path.exists(self.output_file):
                print('made it past the os.path.exists')
                choice = input(f"{self.output_file} already exists. Do you want to overwrite it? (y/n): ")
                if choice.lower() == 'n':
                    output_file = input('Enter a new output file name: ')
                    data = read_input_from_user()
                    go = mrkHolloway(output_file, data)
                    go.run()
                    go.additional_runs()
                else:
                    data = read_input_from_user()
                    go = mrkHolloway(self.output_file, data)
                    go.run()
                    go.additional_runs()
        else:
            sys.exit(0)


    def mrkmix(self, tk):
        """_summary_

        Args:
            tk (_type_): _description_

        Returns:
            _type_: _description_
        """
        
        r = 82.05
        if (tk < 1e-4):
            t = 1.0
        else:
            t = tk
            5
            
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
        co2h2o += np.sqrt(sum([mol.mrk_molecular_interaction
                               for mol in self.molecules if mol.name in ('CO2', 'H20')]))
        asum = 0.0
        bsum = 0.0
        for gas_i in self.molecules:
            bsum += gas_i.mrk_molecular_volume * gas_i.mole_fraction
            for gas_j in self.molecules:
                if gas_i == gas_j:
                    if gas_i.name != 'H20' and gas_j.name != 'CO2':
                        asum += gas_i.mole_fraction * gas_j.mole_fraction * gas_i.mrk_molecular_interaction
                    elif gas_i.name == 'H20':
                        asum += gas_i.mole_fraction * gas_j.mole_fraction * ah2om
                    elif gas_i == 'CO2':
                        asum += gas_i.mole_fraction * gas_j.mole_fraction * aco2m
                elif ((gas_i == 'H20') and (gas_j == 'CO2') or (gas_i == 'CO2') and (gas_j == 'H20')):
                    asum += gas_i.mole_fraction * gas_j.mole_fraction * co2h2o
                else:
                    asum += gas_i.mole_fraction * gas_j.mole_fraction * \
                    np.sqrt(gas_i.mrk_molecular_interaction * gas_j.mrk_molecular_interaction)
          
        asum /= 1.013
        return asum, bsum
    
    
    def isochore_plot(self, num_output):
        import matplotlib.pyplot as plt
        plt.style.use('./mystyle.mplstyle')
        
        # # Sort the data by density
        # self.data_by_density = {}
        # for idx, d in enumerate(self.density):
        #     self.data_by_density[d] = {}
        #     for temp_idx, temp in enumerate(self.temps_k):
        #         self.data_by_density[d][temp-273.15] = self.temp_p_data[temp-273.15][idx]
                
        # split the dictionary to be plotting temp v. pressure      
        # for density in self.data_by_density:
        #     data = self.data_by_density[density]
        #     lists = sorted(data.items()) # sorted by key, return a list of tuples

        #     temps, pressures = zip(*lists) # unpack a list of pairs into two tuples
        temps = [x - 273.15 for x in self.temps_k]
        plt.plot(temps, self.pressure_output, color='black')
        
        # Set the axis labels
        plt.xlabel('Temperature (c)')
        plt.ylabel('Pressure (bar)')
        
        # Update x and y limits
        plt.xlim(0, max(temps))
        plt.ylim(0, max(self.pressure_output)) # enter max pressure you want to see 
        # Show the plot
        
        plt.savefig(f'{self.output_file}_figure_{num_output}.png')
        plt.show()
        

    def write_output(self):
        self.output_file += '_output.txt'
        with open(self.output_file, 'w') as f:
            f.write(f'The starting temperature = {self.tstart:3f} deg. c\n')
            f.write(f'The Molar volume = {self.vstart:3f} cm^3/mol\n')
            f.write('The mole fractions are:\n')
            formstr = ' '.join(f"{x.name:<9}" for x in self.molecules)
            mix_str = ' '.join(f"{x.mole_fraction:<9.3f}" for x in self.molecules)
            f.write(formstr + '\n')
            f.write(mix_str + '\n')
            f.write(f'{"":<30s} Isochore pressure in bars\n')
            f.write(f'{"T deg. c  mol vol =":<30s}' + 
                    ' '.join(f"{self.vstart:<10.3f}") + '\n')
            f.write(f'{"density =":<30s}' + 
                    ' '.join(f'{self.density:<10.4f}') + '\n')
            
            # Print temp/pressure outputs
            for t in self.temps_k:
                f.write(f'{t - 273.15:10.0f} {"":<13s} ' + 
                    ''.join(f'{x:10.0f}' for x in self.pressure_output) + '\n')
            f.close()
       
            
    def run(self, num_output=1):
        """Run everything calculate isochore and plot and print output to screen
        """
        self.check_totals()
        print('The mole fractions are:')
        formstr = ' '.join(f"{x.name:<9}" for x in self.molecules)
        mix_str = ' '.join(f"{x.mole_fraction:<9.3f}" for x in self.molecules)
        
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
              ' '.join(f"{self.vstart:<10.3f}" ))
        
        print(f'{"density =":<30s}' + 
              ' '.join(f'{self.density:<10.4f}'))
        
        # Calc temp pressure
        self.temp_pressure_calc(tk)
        self.isochore_plot(num_output)
        self.write_output()
        
        
def create_molecules(row):
    """Takes input from file and creates the Molecule dataclass for the gases

    Args:
        row (_type_): _description_

    Returns:
        _type_: _description_
    """
    molecules = [
    Molecule(name='CO2', weight=44.0, mrk_molecular_interaction=46.e6, mrk_molecular_volume=2.97e1, mole_fraction=row['CO2']),
    Molecule(name='CO', weight=28.0, mrk_molecular_interaction=16.98e6, mrk_molecular_volume=2.738e1, mole_fraction=row['CO']),
    Molecule(name='CH4', weight=16.0, mrk_molecular_interaction=31.59e6, mrk_molecular_volume=2.9703e1, mole_fraction=row['CH4']),
    Molecule(name='H2', weight=2.0, mrk_molecular_interaction=3.56e6, mrk_molecular_volume=1.515e1, mole_fraction=row['H2']),
    Molecule(name='H2O', weight=18.0, mrk_molecular_interaction=35.e6, mrk_molecular_volume=1.46e1, mole_fraction=row['H2O']),
    Molecule(name='H2S', weight=34.0, mrk_molecular_interaction=87.9e6, mrk_molecular_volume=2.0e1, mole_fraction=row['H2S']),
    Molecule(name='SO2', weight=64.06, mrk_molecular_interaction=142.6e6, mrk_molecular_volume=3.94e1, mole_fraction=row['SO2']),
    Molecule(name='N2', weight=28.0, mrk_molecular_interaction=15.382e6, mrk_molecular_volume=2.68e1, mole_fraction=row['N2'])]
    return molecules


def read_input_from_user():
    """Takes input from terminal if user isn't supplying an input file

    Returns:
        _type_: _description_
    """
    
    data = {}
    temp = input('Enter T (deg C) \n') 
    volume = input("Enter molar volume \n")
    data['tstart'] = float(temp)
    data['max_t'] = float(input('What is the maximum temperature to calculate?\n'))
    data['t_increment'] = int(input('What increment of degrees?\n'))
    
    if data['tstart'] < 0.01:
        sys.exit()
    data['vstart'] = float(volume)
    data['molecules'] = [
        Molecule(name='CO2', weight=44.0, mrk_molecular_interaction=46.e6, mrk_molecular_volume=2.97e1, mole_fraction=float(input(f'Enter mole fraction CO2\n'))),
        Molecule(name='CO', weight=28.0, mrk_molecular_interaction=16.98e6, mrk_molecular_volume=2.738e1, mole_fraction=float(input('Enter mole fraction CO \n'))),
        Molecule(name='CH4', weight=16.0, mrk_molecular_interaction=31.59e6, mrk_molecular_volume=2.9703e1, mole_fraction=float(input('Enter mole fraction CH4 \n'))),
        Molecule(name='H2', weight=2.0, mrk_molecular_interaction=3.56e6, mrk_molecular_volume=1.515e1, mole_fraction=float(input('Enter mole fraction H2\n'))),
        Molecule(name='H2O', weight=18.0, mrk_molecular_interaction=35.e6, mrk_molecular_volume=1.46e1, mole_fraction=float(input('Enter mole fraction H20\n'))),
        Molecule(name='H2S', weight=34.0, mrk_molecular_interaction=87.9e6, mrk_molecular_volume=2.0e1, mole_fraction=float(input('Enter mole fraction H2S \n'))),
        Molecule(name='SO2', weight=64.06, mrk_molecular_interaction=142.6e6, mrk_molecular_volume=3.94e1, mole_fraction=float(input('Enter mole fraction SO2 \n'))),
        Molecule(name='N2', weight=28.0, mrk_molecular_interaction=15.382e6, mrk_molecular_volume=2.68e1, mole_fraction=float(input('Enter mole fraction N2 \n')))]
     
    return data
     
        
def read_input_from_file(filename):
    """If at file run a file is offered read all of the inputs from the file except the name of the output file

    Args:
        filename (_type_): _description_

    Returns:
        _type_: _description_
    """
    df = pd.read_csv(filename)
    data = []
    for idx, row in df.iterrows():
        tstart = row['starting_temp']
        vstart = row['starting_volume']
        max_t = row['max_temp']
        max_v = row['max_volume']
        num_temps = int(row['temp_increments'])
        num_volume = int(row['volume_increments'])
        if tstart < 0.01:
            sys.exit()
        molecules = [
            Molecule(name='CO2', weight=44.0, mrk_molecular_interaction=46.e6, mrk_molecular_volume=2.97e1, mole_fraction=row['CO2']),
            Molecule(name='CO', weight=28.0, mrk_molecular_interaction=16.98e6, mrk_molecular_volume=2.738e1, mole_fraction=row['CO']),
            Molecule(name='CH4', weight=16.0, mrk_molecular_interaction=31.59e6, mrk_molecular_volume=2.9703e1, mole_fraction=row['CH4']),
            Molecule(name='H2', weight=2.0, mrk_molecular_interaction=3.56e6, mrk_molecular_volume=1.515e1, mole_fraction=row['H2']),
            Molecule(name='H2O', weight=18.0, mrk_molecular_interaction=35.e6, mrk_molecular_volume=1.46e1, mole_fraction=row['H2O']),
            Molecule(name='H2S', weight=34.0, mrk_molecular_interaction=87.9e6, mrk_molecular_volume=2.0e1, mole_fraction=row['H2S']),
            Molecule(name='SO2', weight=64.06, mrk_molecular_interaction=142.6e6, mrk_molecular_volume=3.94e1, mole_fraction=row['SO2']),
            Molecule(name='N2', weight=28.0, mrk_molecular_interaction=15.382e6, mrk_molecular_volume=2.68e1, mole_fraction=row['N2'])]
        mixnum = len(molecules)
        formula_weight = sum(gas.weight * gas.mole_fraction for gas in molecules)
        values = {'tstart': tstart,
                'vstart': vstart,
                'max_t' : max_t,
                'max_v' : max_v,
                'num_temps': num_temps,
                'num_volume': num_volume,
                'molecules': molecules,
                'mixnum': mixnum}
        data = {idx: values}
    return data
                                 

if __name__ == '__main__':
    if len(sys.argv) <= 1:
        output_file = input('Enter output filename \n')
        data = read_input_from_user()
        go = mrkHolloway(output_file, data)
        go.run()
        go.additional_runs()
        
    else:
        output_file_name = input('Enter output filename \n')
        df = pd.read_csv('isochore_data.csv')
        df['molecules'] = df.apply(create_molecules, axis=1)
        for idx, row in df.iterrows():
            data = {
                'tstart': row['starting_temp'],
                'vstart': row['starting_volume'],
                't_increment': row['temp_increments'],
                'max_t': row['max_temp'],
                'molecules': row['molecules']
            }
            output_file = f'{output_file_name}_{idx}'
            num_output = len(df)
            go = mrkHolloway(output_file, data)
            go.run(num_output)
       
        
        