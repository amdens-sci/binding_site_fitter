import numpy as np, matplotlib
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
from matplotlib.pyplot import cm
import matplotlib.colors as colors
import matplotlib.pyplot as plt


def cell_text(qtapp):
  num_params = qtapp.current_model.model_type * 2 + 2
  labels = ['kd1', 'p1', 'kd2', 'p2'][0:num_params]
  celltext = [[labels[i], calc_sig_figs(qtapp.current_model.params[i]),
               calc_sig_figs(qtapp.current_model.param_errors[i])]
              for i in range(0, num_params)]
  celltext.append(['','',''])
  celltext.append(['AIC', calc_sig_figs(qtapp.current_model.AIC), ''])
  return celltext
  
def calc_sig_figs(parameter):
  num_sig_figs = 3
  return '{:g}'.format(float('{:.{p}g}'.format(parameter, p=num_sig_figs)))
  

def gen_residual_plot(qtapp):
  popup_plot = plt.figure()
  popup_ax = popup_plot.add_subplot(111)
  predicted_free, _ = qtapp.current_model.make_preds(qtapp.current_model.associated_data['total'].values)
  actual_free = qtapp.current_model.associated_data['free'].values
  popup_ax.scatter(qtapp.current_model.input_x_values, 100*(1 - predicted_free / actual_free), marker = 'o',
                     facecolors='None', edgecolors='red', s=8)
  horizontal_array = np.arange(qtapp.current_model.input_x_values[0],qtapp.current_model.input_x_values[-1],1)
  popup_ax.plot(horizontal_array, np.zeros((horizontal_array.shape[0])), color='black')
  popup_ax.set_xscale('log')
  popup_ax.set_title('Residuals as percent of actual for the current fit')
  popup_ax.set_xlabel('Total drug concentration from model-associated dataset (ug/mL)')
  popup_ax.set_ylabel('% error on predicted free drug')
  popup_plot.show()


def gen_plot(qtapp, generate_residual_popup_plot = False):
  qtapp.central_plot.clf()
  ax = qtapp.central_plot.add_subplot(121)
  ax2 = qtapp.central_plot.add_subplot(122)
  ax2.clear()
  ax.clear()

  ax.set_ylabel('Free drug concentration (ug/mL)')
  ax.set_xlabel('Total drug concentration (ug/mL)')
  ax.set_yscale('log')
  ax.set_xscale('log')

  xactual = qtapp.current_model.associated_data['total'].values
  yactual = qtapp.current_model.associated_data['free'].values
  ypredicted = qtapp.current_model.predicted_free

  ax.scatter(xactual, yactual, color='k', s=1.0)
  ax.plot(qtapp.current_model.x_values_spanning_range, qtapp.current_model.predicted_free, 
      color=qtapp.current_model.plot_color, linewidth=0.5)
  ax2.axis('off')
  ax2.table(cellText = cell_text(qtapp), cellLoc='center', loc='center right',
            colLabels = ['Parameter', 'Value', 'Uncertainty'],
            colWidths=[0.33, 0.33, 0.33])
  ax2.set_title('Parameter estimates, associated error and other\nstatistics for current model')
  if qtapp.current_model.plot_title is not None:
                ax.set_title(qtapp.current_model.plot_title)
  qtapp.canvas.draw()
  return '0'
