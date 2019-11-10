from PyQt5.QtWidgets import QWidget, QProgressBar, QDialog
from PyQt5.QtCore import QBasicTimer

class ProgBar(QDialog):
    def __init__(self):
        super().__init__()
        self.central_widget = QWidget()
        self.setCentralWidget(self.central_widget)
        mainlayout = QVBoxLayout(self.central_widget)


        self.progress_indicator = QProgressBar(self)
        self.progress_indicator.setGeometry(30,40,300,50)
        self.setWindowTitle('Fitting progress')
        mainlayout.addWidget(self.progress_indicator)

        self.progress_indicator.setValue(0)
        self.progress_indicator.setMaximum(100)
        self.show()
