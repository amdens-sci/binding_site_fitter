from PyQt5 import QtCore
from PyQt5.QtWidgets import QApplication, QLabel, QWidget, QPushButton, QVBoxLayout, QMainWindow, QMessageBox, QFileDialog
from PyQt5.QtWidgets import QLineEdit, QHBoxLayout, QCheckBox, QComboBox
from PyQt5.QtGui import QPixmap
import pandas as pd, numpy as np, scipy, model_building, plotting, pickle, os
from scipy import optimize
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure

class protfitter(QMainWindow):

  def __init__(self):
    super().__init__()
    self.current_model = model_building.protbind_model()
    self.current_data = pd.DataFrame()
    self.weight_options = {'No weighting':0, '1/C weighting':1, '1/(C^2) weighting':2,
                           '1/(C^3) weighting':3}
    self.regression_options = {'One-binding site model':0, 'Two-binding site model':1}
    
    self.current_model_selection = self.regression_options['One-binding site model']
    self.current_weight_selection = self.weight_options['No weighting']
    self.color_options = {'Blue':'b', 'Green':'g', 'Red':'r'}
    self.color_index = {'b':0, 'g':1, 'r':2}

    self.central_plot = Figure(figsize=(15,5))
    self.canvas = FigureCanvas(self.central_plot)
    self.toolbar = NavigationToolbar(self.canvas, self)
    
    self.setWindowTitle('Protein Binding Fitter')
    self.central_widget = QWidget()
    self.setCentralWidget(self.central_widget)
    mainlayout = QVBoxLayout(self.central_widget)
    self.central_widget.setFixedWidth(1200)
    top_panel = QHBoxLayout()
    left_panel = QVBoxLayout()
    right_panel = QVBoxLayout()
    left_panel.addStretch(1)

    right_panel.addWidget(self.toolbar)
    right_panel.addWidget(self.canvas)

    data_label = QLabel('Data options')
    left_panel.addWidget(data_label)

    import_button = QPushButton('Import new protein binding dataset')
    left_panel.addWidget(import_button)
    import_button.clicked.connect(self.load_pb_file)

    fit_button = QPushButton('Fit/Plot protein binding data')
    left_panel.addWidget(fit_button)
    fit_button.clicked.connect(self.fit_pb_data)

    model_label = QLabel('Model fitting options')
    left_panel.addWidget(model_label)
    
    self.regression_type = QComboBox()
    self.regression_type.addItem("One-binding site model")
    self.regression_type.addItem("Two-binding site model")
    self.regression_type.activated[str].connect(self.change_regression_type)
    left_panel.addWidget(self.regression_type)

    self.weight_type = QComboBox()
    self.weight_type.addItem("No weighting")
    self.weight_type.addItem("1/C weighting")
    self.weight_type.addItem("1/(C^2) weighting")
    self.weight_type.addItem("1/(C^3) weighting")
    self.weight_type.activated[str].connect(self.change_weight_type)
    left_panel.addWidget(self.weight_type)

    self.color_palette = QComboBox()
    self.color_palette.addItem("Blue")
    self.color_palette.addItem("Green")
    self.color_palette.addItem("Red")
    self.color_palette.activated[str].connect(self.change_color_palette)
    left_panel.addWidget(self.color_palette)

    save_label = QLabel('Model options')
    left_panel.addWidget(save_label)

    save_button = QPushButton('Save Model')
    left_panel.addWidget(save_button)
    save_button.clicked.connect(self.save_model)

    load_button = QPushButton('Load Saved Model')
    left_panel.addWidget(load_button)
    load_button.clicked.connect(self.load_saved_model)

    pk_calc_button = QPushButton('Run calculations for PK data')
    left_panel.addWidget(pk_calc_button)
    pk_calc_button.clicked.connect(self.run_pk_calcs)

    top_panel.addLayout(left_panel)
    top_panel.addLayout(right_panel)
    
    mainlayout.addLayout(top_panel)
    self.show()

  def run_pk_calcs(self):
    alert = QMessageBox()
    alert.setText("For PK data, your input must be a csv file with concentration in the first column. Your file may "
                  "contain other columns, but they will all be ignored. Do not use column headers.")
    alert.setWindowTitle('Protein Fitter')
    alert.exec_()
    options = QFileDialog.Options()
    filename, _ = QFileDialog.getOpenFileName(self,"Load a csv file containing PK data",
            "","csv files (*.csv);;", options=options)
    if filename:
      try:
        input_pk_data = pd.read_csv(filename, header=None)
        total_concentrations = input_pk_data.values[:,0].astype(float)
      except:
        self.sudden_death("There was an error reading your file. Are you sure there's no text in column 1? Are you "
                          "SURE it's a csv?")
      try:
        free_concentrations, error_code = self.current_model.make_preds(total_concentrations)
        output_dict = {'Total concentration':total_concentrations, 'Free concentration':free_concentrations}
        output_df = pd.DataFrame.from_dict(output_dict)
      except:
        self.sudden_death("There was an error generating the output. Data has not been saved. Make sure you have a model "
                          "fitted/loaded before attempting to run PK calculations.")
      options = QFileDialog.Options()
      filename, _ = QFileDialog.getSaveFileName(self,"Save calculated PK concentrations",
            "","csv files (*.csv);;", options=options)
      try:
        if '.csv' not in filename:
          filename = filename + '.csv'
        output_df.to_csv(filename)
      except:
        self.sudden_death("There was an error generating the output file. Data has not been saved.")
    alert = QMessageBox()
    alert.setText("Congratulations! Your results are saved to %s"%filename)
    alert.setIconPixmap(QPixmap(os.path.join('lib', 'success.jpg')))
    alert.setWindowTitle('Protein Fitter')
    alert.exec_()
        


  def save_model(self):
    if self.current_model.predicted_free is None:
        self.sudden_death('You must fit a model before you can save one. Nothing saved.')
        return
    self.current_model.plot_title = self.central_plot.axes[0].title._text
    options = QFileDialog.Options()
    filename, _ = QFileDialog.getSaveFileName(self,"Save model",
            "","MODEL Files (*.model);;", options=options)
    if filename:
      try:
        with open('%s.model'%filename, 'wb') as output_model:
          pickle.dump(self.current_model, output_model)
      except:
        self.sudden_death('You are trying to overwrite an existing file or '
                          'save in an invalid location. Model not saved.')
        return

  def load_saved_model(self):
    options = QFileDialog.Options()
    filename, _ = QFileDialog.getOpenFileName(self,"Load model",
            "","MODEL Files (*.model);;", options=options)
    if filename:
      try:
        with open(filename, 'rb') as input_model:
          self.current_model = pickle.load(input_model)
      except:
        self.sudden_death('There was an error loading the model you indicated. Try again, or check '
                          'to make sure you indicated a valid .mod file.')
        return
      plotting.gen_plot(self)
      self.regression_type.setCurrentIndex(self.current_model.model_type)
      self.weight_type.setCurrentIndex(self.current_model.weight_type)
      self.color_palette.setCurrentIndex(self.color_index[self.current_model.plot_color])

      
  def change_weight_type(self, text):
    if text in self.weight_options:
      self.current_weight_selection = self.weight_options[text]
    else:
      self.sudden_death('Invalid option')

  def change_regression_type(self, text):
    if text in self.regression_options:
      self.current_model_selection = self.regression_options[text]
    else:
      self.sudden_death('Invalid option')

  def change_color_palette(self, text):
    if text in self.color_options:
      self.current_model.plot_color = self.color_options[text]
      plotting.gen_plot(self)
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
    self.current_model.model_type = self.current_model_selection
    self.current_model.weight_type = self.current_weight_selection
    output_code = self.current_model.model_fit()
    if output_code != '0':
      self.sudden_death(output_code)
      return
    else:
      plotting.gen_plot(self)
      plotting.gen_residual_plot(self)



  def sudden_death(self, error_message):
    alert = QMessageBox()
    alert.setText(error_message)
    alert.setIconPixmap(QPixmap(os.path.join('lib', 'failure.jpg')))
    alert.setWindowTitle('Sudden Death')
    alert.exec_()


if __name__ == '__main__':
  app = QApplication([])
  current_app = protfitter()
  app.exec_()
  exit()
