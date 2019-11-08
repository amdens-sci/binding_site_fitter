from PyQt5 import QtCore
from PyQt5.QtWidgets import QApplication, QLabel, QWidget, QPushButton, QVBoxLayout, QMainWindow, QMessageBox, QFileDialog, QLineEdit, QHBoxLayout, QCheckBox, QComboBox
from PyQt5.QtGui import QPixmap
import pandas as pd, numpy as np, scipy, model_building, plotting
from scipy import optimize
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure

class protfitter(QMainWindow):

	def __init__(self):
		super().__init__()
		self.current_model = model_building.protbind_model()

		self.weight_options = {'No weighting':0, '1/C weighting':1, '1/(C^2) weighting':2}
		self.regression_options = {'One-binding site model':0, 'Two-binding site model':1,
						'Three-binding site model':2}
		self.color_options = {'Blue':'b', 'Green':'r', 'Red':'g'}
		self.color_choice = 'b'
		self.error_options = {'Fisher information':0, '1000 iteration bootstrap':1}

		self.central_plot = Figure(figsize=(8,6))
		self.canvas = FigureCanvas(self.central_plot)
		self.toolbar = NavigationToolbar(self.canvas, self)
		
		self.setWindowTitle('Protein Binding Fitter')
		self.central_widget = QWidget()
		self.setCentralWidget(self.central_widget)
		mainlayout = QVBoxLayout(self.central_widget)
		horiz_layout_top = QHBoxLayout()
		horiz_layout_bottom = QHBoxLayout()

		mainlayout.addWidget(self.toolbar)
		mainlayout.addWidget(self.canvas)


		import_button = QPushButton('Import new protein binding dataset')
		mainlayout.addWidget(import_button)
		import_button.clicked.connect(self.load_pb_file)

		fit_button = QPushButton('Fit/Plot protein binding data')
		fit_button.resize(100, fit_button.height())
		mainlayout.addWidget(fit_button)
		fit_button.clicked.connect(self.fit_pb_data)


		self.regression_type = QComboBox()
		self.regression_type.addItem("One-binding site model")
		self.regression_type.addItem("Two-binding site model")
		self.regression_type.addItem("Three-binding site model")
		self.regression_type.activated[str].connect(self.change_regression_type)
		horiz_layout_top.addWidget(self.regression_type)

		self.weight_type = QComboBox()
		self.weight_type.addItem("No weighting")
		self.weight_type.addItem("1/C weighting")
		self.weight_type.addItem("1/(C^2) weighting")
		self.weight_type.activated[str].connect(self.change_weight_type)
		horiz_layout_top.addWidget(self.weight_type)

		self.color_palette = QComboBox()
		self.color_palette.addItem("Blue")
		self.color_palette.addItem("Green")
		self.color_palette.addItem("Red")
		self.color_palette.activated[str].connect(self.change_color_palette)
		horiz_layout_bottom.addWidget(self.color_palette)

		self.error_type = QComboBox()
		self.error_type.addItem("Fisher information")
		self.error_type.addItem("1000 iteration bootstrap")
		self.error_type.activated[str].connect(self.change_error_type)
		horiz_layout_bottom.addWidget(self.error_type)

		mainlayout.addLayout(horiz_layout_top)
		mainlayout.addLayout(horiz_layout_bottom)
		self.show()






	def change_error_type(self, text):
		if text in self.error_options:
			self.current_model.error_type = self.error_options[text]
		else:
			self.sudden_death('Invalid option')

	def change_weight_type(self, text):
		if text in self.weight_options:
			self.current_model.weight_type = self.weight_options[text]
		else:
			self.sudden_death('Invalid option')

	def change_regression_type(self, text):
		if text in self.regression_options:
			self.current_model.model_type = self.regression_options[text]
		else:
			self.sudden_death('Invalid option')

	def change_color_palette(self, text):
		if text in self.color_options:
			self.color_choice = self.color_options[text]
		else:
			self.sudden_death('Invalid option')

	def load_pb_file(self):
		alert = QMessageBox()
		alert.setText("A quick reminder: When importing data, use csv files where the first column contains the total concentration "
			"and the second column contains free drug concentration, otherwise you'll get very strange results.")
		alert.setWindowTitle('Protein Fitter 2000')
		alert.exec_()
		options = QFileDialog.Options()
		filename, _ = QFileDialog.getOpenFileName(self,"Load new dataset",
						"","CSV Files (*.csv);;", options=options)
		if filename:
			try:
				self.current_model.associated_data = pd.read_csv(filename, header=None)
				self.current_model.associated_data.columns = ['total', 'free']
			except:
				self.sudden_death('There was an error opening the selected file! Clearly you have made a mistake. '
					'One reason why this may have occurred '
					'is if you selected a non-csv file or a file with more than two columns. '
					'Remember your instructions!')
			if len(self.current_model.associated_data.columns) != 2:
				self.sudden_death('There was an error opening the selected file! Clearly you have made a mistake. '
					'One reason why this may have occurred '
					'is if you selected a non-csv file or a file with more than two columns. '
					'Remember your instructions!')

	def fit_pb_data(self):
		if self.current_model.associated_data is None:
			self.sudden_death("You want to fit the data, but you haven't loaded any? Try loading some first. Now there's an idea!")
			return
		output_code = self.current_model.model_fit()
		if output_code != '0':
			self.sudden_death(output_code)
			return
		plotting.gen_plot(self)


	def save_pb_model(self):
		options = QFileDialog.Options()
		filename, _ = QFileDialog.getSaveFileName(self,"Save Model / Results",
						"",".mod Files (*.mod);;", options=options)
		if filename:
			error_code = data_export.export_results(self, filename)
			if error_code != '0':
				self.sudden_death(error_code)
				return
			alert = QMessageBox()
			alert.setText('Your model has been exported to a .mod file entitled "%s" . Results are now plotted to screen.'%filename)
			alert.setWindowTitle('Success!')
			alert.exec_()


	def sudden_death(self, error_message):
		alert = QMessageBox()
		alert.setText(error_message)
		alert.setWindowTitle('Sudden Death')
		alert.exec_()


if __name__ == '__main__':
	app = QApplication([])
	current_app = protfitter()
	app.exec_()
	exit()