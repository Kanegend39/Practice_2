from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure


class MplCanvas_3D(FigureCanvas):
    def __init__(self, parent=None, width=10, height=8, dpi=100):
        fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = fig.add_subplot(111, projection='3d')
        self.axes.set_facecolor('white')
        self.axes.grid(color='black', linewidth=0.5)
        self.axes.grid(True)
        super(MplCanvas_3D, self).__init__(fig)
