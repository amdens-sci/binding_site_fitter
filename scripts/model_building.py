import numpy as np, scipy, pandas as pd, multiprocessing as mp, numdifftools
from scipy import optimize
from scipy.optimize import minimize
from sklearn.model_selection import ParameterGrid
import ctypes, os
from numpy.ctypeslib import ndpointer
from PyQt5.QtWidgets import QMessageBox

#We need to load the dll / so and set up two specific C++ functions for use in this
#program. We do this using ctypes. The C++ functions will take 4 arguments: pointer
#to the input numpy array, dimensions of the input array and pointer to the output numpy array
#which the C++ functions will assume is correctly sized. THe output numpy array will be modified
#in place so nothing needs to be returned. Extremely important the numpy arrays are C-contiguous!
lib = np.ctypeslib.load_library('rootfinder.so', os.path.join('..','lib'))
rootfind2 = lib.single_pos_special_cubic
rootfind2.restype = None
rootfind2.argtypes = [ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                ctypes.c_int, ctypes.c_int,
                ndpointer(ctypes.c_double, flags="C_CONTIGUOUS")]

rootfind3 = lib.single_pos_special_quartic
rootfind3.restype = None
rootfind3.argtypes = [ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                ctypes.c_int, ctypes.c_int,
                ndpointer(ctypes.c_double, flags="C_CONTIGUOUS")]

#This is a generic model which may belong to one of several types (single binding site, two binding site or three binding
#site). If it is single or two binding site only some of the parameters in the parameter dictionary will be used.
#Model type corresponds to 0 = single binding site, 1 = two binding site and 2 = three binding site.
#Weight type corresponds to 0 = none, 1 = 1/C and 2 = 1/C**2.
class protbind_model:
        def __init__(self, model_type = 0, weight_type = 0, error_type = 0):
                self.current_dataset = pd.DataFrame()
                self.param_ids = {'kd1':0, 'p1':1, 'kd2':2,
                                  'p2':3}
                self.params = np.zeros((4))
                self.param_errors = np.zeros((4))
                self.model_type = model_type
                self.model_associated_params = {0:['kd1', 'p1'],
                                                1:['kd1', 'p1', 'kd2', 'p2']}
                self.weight_type = weight_type
                self.error_type = error_type
                self.associated_data = None
                self.predicted_free = None
                self.x_values_spanning_range = None
                self.AIC = 0
                self.plot_color = 'b'
                self.plot_title = None


###This function enables a fitted model to make predictions on demand.
        def make_preds(self, xvalues):
                if self.associated_data is not None and self.x_values_spanning_range is not None:
                        if self.model_type == 0:
                                return (self.onebind(xvalues, self.params[0:2])), '0'
                        if self.model_type == 1:
                                return (self.twobind(xvalues, self.params)), '0'
                else:
                        return 0, 'err'


###These functions are for the single-site binding model

        def onebind_resid(self, params, x, actual):
                residuals = (self.onebind(x, params) - actual)
                if self.weight_type == 0:
                        return np.dot(residuals, residuals)
                elif self.weight_type == 1:
                        return np.sum(residuals * residuals * 1/x)
                elif self.weight_type == 2:
                        return np.sum(residuals * residuals * 1/(x**2))
                elif self.weight_type == 3:
                        return np.sum(residuals * residuals * 1/(x**3))

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
                elif self.weight_type == 3:
                        return np.sum(residuals * residuals * 1/(x**3))
                

        def twobind(self, x, params):
                kd1, kd2 = params[self.param_ids['kd1']], params[self.param_ids['kd2']]
                p1, p2 = params[self.param_ids['p1']], params[self.param_ids['p2']]
                coefficients = np.zeros((x.shape[0], 4))
                coefficients[:,-1] = 1.0
                coefficients[:,-2] = p1 + p2 + kd1 + kd2 - x
                coefficients[:,-3] = p1*kd2 + p2*kd1 + kd1*kd2 - kd2*x - kd1*x
                coefficients[:,-4] = -kd1*kd2*x
                predicted = np.zeros((x.shape[0]))
                rootfind2(coefficients, x.shape[0], 4, predicted)
                return predicted




        ##This function is to fit a new model. Depending on the user's selected option for model fitting, we may
        #call one of several functions specific to certain model types.
        def model_fit(self):
                try:
                        x = np.asarray([float(z) for z in self.associated_data['total'].values])
                        y = np.asarray([float(z) for z in self.associated_data['free'].values])
                except:
                        return ('could not convert data to numeric form! Are you sure your '
                                'input protein binding dataset is only non-numeric characters?')
                self.x_values_spanning_range = np.exp(np.arange(np.log(np.min(x)), np.log(np.max(x)), 0.1) )
                residual_function_dict = {0:self.onebind_resid, 1:self.twobind_resid}
                preds_function_dict = {0:self.onebind, 1:self.twobind}
                
                self.binding_fitter(x, y, residual_function_dict[self.model_type],
                                    preds_function_dict[self.model_type])
                self.actual_free = np.copy(y)
                self.input_x_values = np.copy(x)
                error_code = self.calc_error(x, y, residual_function_dict[self.model_type])
                if error_code == '0':
                        return '0'
                else:
                        return error_code
                

        def binding_fitter(self, x, y, residual_function, pred_function):
                print('Beginning fit')
                self.grid_search(x,y,residual_function)
                self.predicted_free, _ = self.make_preds(self.x_values_spanning_range)
                print('Fit complete.')

        def grid_search(self, x, y, residual_function):
                self.params = np.zeros((4))
                if self.model_type == 0:
                        param_grid = {'kd1':[1, 10.0, 100], 'p1':[5.0, 50, 500]}
                if self.model_type == 1:
                        param_grid = {'kd1':[1, 10.0, 100], 'p1':[1.0, 10, 100],
                                      'kd2':[10, 100, 1000], 'p2':[10.0, 100, 1000]}
                        
                start_parameters = ParameterGrid(param_grid)
                bounds = [(1e-10, None) for j in range(0, len(self.model_associated_params[self.model_type]))]
                best_params = np.zeros((len(bounds)))
                best_fun = -1
                for j in range(0, len(start_parameters)):
                        current_starting_params = [start_parameters[j][param] for param in
                                                   self.model_associated_params[self.model_type]]
                        ssbm = minimize(residual_function, x0=current_starting_params, args=(x, y), 
                                                                method='L-BFGS-B', bounds=bounds)
                        if ssbm.success == True:
                                if best_fun < 0 or ssbm.fun < best_fun:
                                        best_fun = ssbm.fun
                                        best_params = ssbm.x
                        if j % 25 == 0:
                                print('%s iterations'%j)
                item_counter = 0
                for j, current_param in enumerate(self.model_associated_params[self.model_type]):
                        self.params[self.param_ids[current_param]] = best_params[j]
                print(self.params)


        def calc_error(self, x, y, residual_fun):
                net_residuals = residual_fun(self.params, x, y)
                num_params = len(self.model_associated_params[self.model_type])
                degrees_of_freedom = x.shape[0] - num_params
                hess_calculator = numdifftools.Hessian(residual_fun)
                hessian = hess_calculator(self.params[0:num_params], x, y)
                hessian_eigenvalues = np.linalg.eig(hessian)[0]
                cov_mat = np.linalg.inv(hessian)
                if np.min(cov_mat) < 0:
                        return ("There was an error inverting the Hessian; the problem is likely ill-posed."
                                      "Try using a less flexible model with fewer free-to-vary parameters.")
                for i in range(0, num_params):
                        self.param_errors[i] = np.sqrt((net_residuals / degrees_of_freedom)
                                                               * cov_mat[i,i])
                self.AIC = net_residuals + 2*num_params + 2*num_params*(num_params+1)/(x.shape[0]-num_params-1)
                if (np.max(hessian_eigenvalues) / np.min(hessian_eigenvalues)) > 1e8:
                        alert = QMessageBox()
                        alert.setText("The Hessian is ill-conditioned; error calculations may not be reliable."
                                      "The problem is likely ill-posed. Do not use this result.")
                        alert.setWindowTitle('Protein Fitter')
                        alert.exec_()
                return '0'
