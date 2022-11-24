import numpy as np
import os.path
import matplotlib.pyplot as plt
from map_tools import *

# plot parameters
plt.rcParams["font.size"]=10
plt.rcParams["lines.markersize"]=3
plt.rcParams["figure.figsize"]=[12,8]
plt.rcParams["figure.dpi"]=200

class Temperature:
    
    def __init__(self, run, cyc, calib=False):
        self.path = "../maps_bin/{}_{}_000_temperature.EDMdat".format(str(run).zfill(6), str(cyc).zfill(6)) 
        if os.path.exists(self.path):
            self.getData()
            self.file_exists=True
        else:
            self.file_exists=False
            return None
        if calib==True:
            self.calibrate()
            self.calibrated=True
        else:
            self.calibrated=False
        self.temps = np.array([self.mid_over, self.front_bot, self.back_top, self.front_top, self.mid_bot, self.mid_top, self.back_bot, self.mid_under])
        self.avg = np.mean(self.temps, axis=1)
        self.std = np.std(self.temps, axis=1)
    
    def getData(self, dtformat='v2', ncol=50):
            # define data type according to the header. there is an extra column ('other') not specified in the header.
            self.temp_keys = ["temp {}".format(i) for i in range(1, ncol-1)]
            if dtformat=='v2':
                dtype = np.dtype([("time", np.uint64)] 
                                 + [(key, np.double) for key in self.temp_keys] 
                                 + [("digital inputs", np.uint64)] 
                                 + [("other", np.uint64)])

            with open(self.path, "rb") as f:
                data = np.fromfile(f, dtype=dtype)
                
            self.data = data
            self.t = 1e-9*np.array(data["time"] - data["time"][0])
            self.mid_over = data["temp 9"]
            self.front_bot = data["temp 10"]
            self.back_top = data["temp 29"]
            self.front_top = data["temp 30"]
            self.mid_bot = data["temp 33"]
            self.mid_top = data["temp 34"]
            self.back_bot = data["temp 35"]
            self.mid_under = data["temp 36"]
            
    def calibrate(self, offsets=[-0.019, -0.005, 0.001]):
        self.mid_bot = self.mid_bot ## calibration reference
        self.mid_top = self.mid_top + offsets[0]
        self.back_bot = self.back_bot + offsets[1]
        self.mid_under = self.mid_under + offsets[2]
            
class TemperatureRun:
    
    def __init__(self, run, cycrange, calib=True):
        cycles = {}
        avgs = []
        stds = []
        temps = {}
        time = []
        temps_init = 0
        for cyc in cycrange:
            cycles[cyc] = Temperature(run, cyc, calib=calib)
            cycle = cycles[cyc]
            if cycle.file_exists:
                avgs.append(cycle.avg)
                stds.append(cycle.std)
                time.append(cycle.data["time"])
                # initialize temps dic with first cycle keys
                if temps_init==0:
                    keys = cycle.temp_keys
                    for k in keys:
                        temps[k] = []
                        temps_init = 1
                for k in keys:
                    temps[k].append(cycle.data[k])
        # flatten temperature data dic entry by entry
        for k in keys:
            temps[k] = np.concatenate(temps[k])
        self.temp_keys = keys
        self.time = 1e-9*np.concatenate(time)
        self.cycles = cycles
        self.temps = temps
        self.cycAvgs = np.transpose(avgs)
        self.cycStds = np.transpose(stds)
        self.avgs = np.mean(avgs, axis=0)
        self.stds = np.mean(stds, axis=0)
        
class TemperatureRunSet:
    
    def __init__(self, runrange, cycrange, calib=True):
        runs = {}
        avgs = []
        stds = []
        temps = {}
        time = []
        temps_init = 0
        for r in runrange:
            run = TemperatureRun(r, cycrange, calib=calib)
            if len(run.avgs)>0:
                runs[r] = run
                avgs.append(run.avgs)
                stds.append(run.stds)
                time.append(run.time)
                # initialize temps dic with first run
                if temps_init==0:
                    keys = run.temp_keys
                    for k in keys:
                        temps[k] = []
                        temps_init = 1
                for k in keys:
                    temps[k].append(run.temps[k])
        # flatten temperature data dic entry by entry
        for k in keys:
            temps[k] = np.concatenate(temps[k])
        self.temp_keys = keys
        self.time = np.concatenate(time)
        self.temps = temps
        self.runs = runs
        self.runAvgs = np.transpose(avgs)
        self.runStds = np.transpose(stds)
            