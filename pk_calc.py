import pandas as pd, numpy as np, model_building
from PyQt5.QtWidgets import QFileDialog
import matplotlib.pyplot as plt

def pk_calculations(pb_model, file_to_load):
    try:
        raw_inputs = pd.read_csv(file_to_load, header=None).values
        x = np.asarray([float(z) for z in raw_inputs[:,0]])
        y = pb_model.make_preds(x)
        options = QFileDialog.Options()
        filename, _ = QFileDialog.getSaveFileName(self,"Save model",
                "","MODEL Files (*.model);;", options=options)
    except:
        return 'There was an error loading the selected file. Changes not saved.'
