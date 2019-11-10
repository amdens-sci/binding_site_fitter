import numpy as np, scipy, pandas as pd
from scipy import optimize
from scipy.optimize import minimize
import progress_Bar

import ctypes, os

import numpy as np

from numpy.ctypeslib import ndpointer


lib = np.ctypeslib.load_library('mylib.so', os.getcwd())
rootfind = lib.find_single_positive_root
rootfind.restype = None
rootfind.argtypes = [ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                ctypes.c_int, ctypes.c_int,
                ndpointer(ctypes.c_double, flags="C_CONTIGUOUS")]

#This is a generic model which may belong to one of several types (single binding site, two binding site or three binding
#site). If it is single or two binding site only some of the parameters in the parameter dictionary will be used.
#Model type corresponds to 0 = single binding site, 1 = two binding site and 2 = three binding site.
#Weight type corresponds to 0 = none, 1 = 1/C and 2 = 1/C**2.
#Error type corresponds to 0 = fisher information, 1 = bootstrap.
class protbind_model:
        def __init__(self, model_type = 0, weight_type = 0, error_type = 0):
                self.current_dataset = pd.DataFrame()
                self.params = {'kd1':0.1, 'kd2':10, 'kd3':100,
                                'p1':10, 'p2':700, 'p3':700}
                self.param_errors = {'kd1':0.0, 'kd2':0.0, 'kd3':0.0,
                                'p1':0.0, 'p2':0.0, 'p3':0.0}
                self.model_type = model_type
                self.weight_type = weight_type
                self.error_type = error_type
                self.associated_data = None
                self.predicted_free = None
                self.x_values_spanning_range = None
                self.lowpreds = None
                self.highpreds = None
                self.num_bootstrap = 1000
                self.stoch_search_num_iters = 40


###These functions are for the single-site binding model

        def onebind_resid(self, params, x, actual):
                residuals = (self.onebind(x, params) - actual)
                if self.weight_type == 0:
                        return np.dot(residuals, residuals)
                elif self.weight_type == 1:
                        return np.sum(residuals * residuals * 1/x)
                elif self.weight_type == 2:
                        return np.sum(residuals * residuals * 1/(x**2))

        def onebind(self, x, params):
                kd = params[0]
                ptot = params[1]
                predicted = -0.5*(kd + ptot - x) + 0.5*np.sqrt( (kd + ptot - x)**2 + 
                                                                4*kd*x)
                return predicted


###These functions are for the two-site binding model

        def twobind_resid(self, params, x, actual):
                residuals = (self.twobind(x, params) - actual)
                if self.weight_type == 0:
                        return np.dot(residuals, residuals)
                elif self.weight_type == 1:
                        return np.sum(residuals * residuals * 1/x)
                elif self.weight_type == 2:
                        return np.sum(residuals * residuals * 1/(x**2))
                

        def twobind(self, x, params):
                kd1, kd2, p1, p2 = params[0], params[1], params[2], params[3]
                coefficients = np.zeros((x.shape[0], 4))
                coefficients[:,-1] = 1.0
                coefficients[:,-2] = p1 + p2 + kd1 + kd2 - x
                coefficients[:,-3] = p1*kd2 + p2*kd1 + kd1*kd2 - kd2*x - kd1*x
                coefficients[:,-4] = -kd1*kd2*x
                predicted = np.zeros((x.shape[0]))
                rootfind(coefficients, x.shape[0], 4, predicted)
                return predicted



        ##This function is to fit a new model. Depending on the user's selected option for model fitting, we may
        #call one of several functions specific to certain model types.
        def model_fit(self):
                try:
                        x = np.asarray([float(z) for z in self.associated_data['total'].values])
                        y = np.asarray([float(z) for z in self.associated_data['free'].values])
                except:
                        return 'could not convert data to numeric form! Are you sure your input protein binding dataset is only non-numeric characters?'
                self.x_values_spanning_range = np.exp(np.arange(np.log(np.min(x)), np.log(np.max(x)), 0.1) )
                if self.model_type == 0:
                        self.onebindsite_fitter(x,y)
                elif self.model_type == 1:
                        self.twobindsite_fitter(x,y)
                elif self.model_type == 2:
                        self.threebindsite_fitter(x,y)
                return '0'


        def onebindsite_fitter(self, x, y):
                start_values = np.zeros((self.stoch_search_num_iters,2))
                start_values[:,0] = np.exp(np.log(0.1) + np.log(500)*np.random.rand(self.stoch_search_num_iters))
                start_values[:,1] = np.exp(np.log(0.5) + np.log(1000)*np.random.rand(self.stoch_search_num_iters))
                best_params = np.zeros((2))
                best_fun = -1
                for j in range(0, self.stoch_search_num_iters):
                        current_starting_params = start_values[j,:]
                        ssbm = minimize(self.onebind_resid, x0=current_starting_params, args=(x, y), 
                                                                method='L-BFGS-B', bounds=[(1e-9,None), (1e-9,None)])
                        if ssbm.success == True:
                                if best_fun < 0 or ssbm.fun < best_fun:
                                        best_fun = ssbm.fun
                                        best_params = ssbm.x
                self.params['kd1'] = best_params[0]
                self.params['p1'] = best_params[1]
                self.predicted_free = self.onebind(self.x_values_spanning_range, best_params)
                print('found optimal solution...now assessing error.')
                if self.error_type == 1 or self.error_type == 0:
                        predictions_stack, parameter_stack = [], []
                        for i in range(0, self.num_bootstrap):
                                current_cases = np.random.choice(np.arange(0,x.shape[0]), x.shape[0], replace=True)
                                ssbm = minimize(self.onebind_resid, x0=best_params, args=(x[current_cases], y[current_cases]), 
                                                                method='L-BFGS-B', bounds=[(1e-9,None), (1e-9,None)])
                                parameter_stack.append(ssbm.x)
                                predictions_stack.append(self.onebind(self.x_values_spanning_range, ssbm.x))
                                if i % 100 == 0:
                                        print(i)
                        predictions_stack = np.sort(np.vstack(predictions_stack), axis=0)
                        parameter_stack = np.vstack(parameter_stack)
                        self.param_errors['kd1'] = np.std(parameter_stack[:,0])
                        self.param_errors['p1'] = np.std(parameter_stack[:,1])
                        self.lowpreds = predictions_stack[25,:]
                        self.highpreds = predictions_stack[975,:]
                        
                                        
        def twobindsite_fitter(self, x, y):
                start_values = np.zeros((self.stoch_search_num_iters,4))
                start_values[:,0] = np.exp(np.log(0.1) + np.log(500)*np.random.rand(self.stoch_search_num_iters))
                start_values[:,1] = np.exp(np.log(0.1) + np.log(500)*np.random.rand(self.stoch_search_num_iters))
                start_values[:,2] = np.exp(np.log(0.5) + np.log(1000)*np.random.rand(self.stoch_search_num_iters))
                start_values[:,3] = np.exp(np.log(0.5) + np.log(1000)*np.random.rand(self.stoch_search_num_iters))
                best_params = np.zeros((4))
                best_fun = -1
                for j in range(0, self.stoch_search_num_iters):
                        current_starting_params = start_values[j,:]
                        ssbm = minimize(self.twobind_resid, x0=current_starting_params, args=(x, y), 
                                                                method='L-BFGS-B', bounds=[(1e-10,None), (1e-10,None),
                                                                                           (1e-10,None), (1e-10,None)])
                        if ssbm.success == True:
                                if best_fun < 0 or ssbm.fun < best_fun:
                                        best_fun = ssbm.fun
                                        best_params = ssbm.x
                        if j % 5 == 0:
                                print(j)
                best_params[0:2] = np.sort(best_params[0:2])
                best_params[2:] = np.sort(best_params[2:])
                self.params['kd1'] = best_params[0]
                self.params['kd2'] = best_params[1]
                self.params['p1'] = best_params[2]
                self.params['p2'] = best_params[3]
                self.predicted_free = self.twobind(self.x_values_spanning_range, best_params)
                print('found optimal solution...now assessing error.')

                if self.error_type == 1 or self.error_type == 0:
                        predictions_stack, parameter_stack = [], []
                        for i in range(0, self.num_bootstrap):
                                current_cases = np.random.choice(np.arange(0,x.shape[0]), x.shape[0], replace=True)
                                ssbm = minimize(self.twobind_resid, x0=best_params, args=(x[current_cases], y[current_cases]), 
                                                                method='L-BFGS-B', bounds=[(1e-10,None), (1e-10,None),
                                                                                           (1e-10,None), (1e-10,None)])
                                parameter_stack.append(ssbm.x)
                                predictions_stack.append(self.twobind(self.x_values_spanning_range, ssbm.x))
                                if i % 100 == 0:
                                        print(i)
                        predictions_stack = np.sort(np.vstack(predictions_stack), axis=0)
                        parameter_stack = np.vstack(parameter_stack)
                        parameter_stack[:,0:2] = np.sort(parameter_stack[:,0:2], axis=1)
                        parameter_stack[:,2:] = np.sort(parameter_stack[:,2:], axis=1)
                        self.param_errors['kd1'] = np.std(parameter_stack[:,0])
                        self.param_errors['kd2'] = np.std(parameter_stack[:,1])
                        self.param_errors['p1'] = np.std(parameter_stack[:,2])
                        self.param_errors['p2'] = np.std(parameter_stack[:,3])
                        self.lowpreds = predictions_stack[25,:]
                        self.highpreds = predictions_stack[975,:]
                        print(self.param_errors)
                        print(self.params)
                
        def threebindsite_fitter(self, x, y):
                return 'This method is not offered yet!!'



