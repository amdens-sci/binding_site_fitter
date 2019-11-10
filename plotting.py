import numpy as np, matplotlib
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
from matplotlib.pyplot import cm
import matplotlib.colors as colors
	



def gen_plot(qtapp):
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
			color=qtapp.color_choice, linewidth=0.5)
	if qtapp.current_model.lowpreds is not None and qtapp.current_model.highpreds is not None:
                ax.fill_between(qtapp.current_model.x_values_spanning_range, qtapp.current_model.lowpreds,
				qtapp.current_model.highpreds, color = qtapp.color_choice,
				alpha=0.4)
	qtapp.canvas.draw()
	return '0'
