from PyQt5 import QtWidgets
from PyQt5.Qt import Qt

import numpy as np
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.patches import Circle
from matplotlib.collections import PatchCollection
import matplotlib.colors as mplc
from scipy.special import comb


class MplWidget(QtWidgets.QWidget):
    def __init__(self, parent=None):
        super(MplWidget, self).__init__(parent)
        self.canvas = FigureCanvas(Figure())

        vertical_layout = QtWidgets.QVBoxLayout(self)
        vertical_layout.addWidget(self.canvas)

        self.canvas.axes = self.canvas.figure.add_subplot(111)
        self.ax0 = self.canvas.axes

        self.num_pts = 0
        self.bezpts = np.zeros((4,2))
        print("self.bezpts=",self.bezpts)
        self.num_eval = 10

        self.plot_xmin = -500
        self.plot_xmax = 500
        self.plot_ymin = -500
        self.plot_ymax = 500

        self.ax0.set_aspect(1.0)
        self.ax0.set_xlim(self.plot_xmin, self.plot_xmax)
        self.ax0.set_ylim(self.plot_ymin, self.plot_ymax)

        self.canvas.mpl_connect("button_press_event", self.button_press)
        # self.canvas.mpl_connect("button_release_event", self.on_release)
        # self.canvas.mpl_connect("motion_notify_event", self.on_move)

    def keyPressEvent(self, event):
        if event.key() == Qt.Key_Space:
            self.test_method()
        elif event.key() == Qt.Key_C:
            self.ax0.cla()
            self.ax0.set_aspect(1.0)
            self.ax0.set_xlim(self.plot_xmin, self.plot_xmax)
            self.ax0.set_ylim(self.plot_ymin, self.plot_ymax)
            self.canvas.update()
            self.canvas.draw()

    def test_method(self):
        print('Space key pressed')

    #---------------------------------------------------------------------------
    def circles(self, x, y, s, c='b', vmin=None, vmax=None, **kwargs):
        """
        See https://gist.github.com/syrte/592a062c562cd2a98a83 

        Make a scatter plot of circles. 
        Similar to plt.scatter, but the size of circles are in data scale.
        Parameters
        ----------
        x, y : scalar or array_like, shape (n, )
            Input data
        s : scalar or array_like, shape (n, ) 
            Radius of circles.
        c : color or sequence of color, optional, default : 'b'
            `c` can be a single color format string, or a sequence of color
            specifications of length `N`, or a sequence of `N` numbers to be
            mapped to colors using the `cmap` and `norm` specified via kwargs.
            Note that `c` should not be a single numeric RGB or RGBA sequence 
            because that is indistinguishable from an array of values
            to be colormapped. (If you insist, use `color` instead.)  
            `c` can be a 2-D array in which the rows are RGB or RGBA, however. 
        vmin, vmax : scalar, optional, default: None
            `vmin` and `vmax` are used in conjunction with `norm` to normalize
            luminance data.  If either are `None`, the min and max of the
            color array is used.
        kwargs : `~matplotlib.collections.Collection` properties
            Eg. alpha, edgecolor(ec), facecolor(fc), linewidth(lw), linestyle(ls), 
            norm, cmap, transform, etc.
        Returns
        -------
        paths : `~matplotlib.collections.PathCollection`
        Examples
        --------
        a = np.arange(11)
        circles(a, a, s=a*0.2, c=a, alpha=0.5, ec='none')
        plt.colorbar()
        License
        --------
        This code is under [The BSD 3-Clause License]
        (http://opensource.org/licenses/BSD-3-Clause)
        """

        if np.isscalar(c):
            kwargs.setdefault('color', c)
            c = None

        if 'fc' in kwargs:
            kwargs.setdefault('facecolor', kwargs.pop('fc'))
        if 'ec' in kwargs:
            kwargs.setdefault('edgecolor', kwargs.pop('ec'))
        if 'ls' in kwargs:
            kwargs.setdefault('linestyle', kwargs.pop('ls'))
        if 'lw' in kwargs:
            kwargs.setdefault('linewidth', kwargs.pop('lw'))
        # You can set `facecolor` with an array for each patch,
        # while you can only set `facecolors` with a value for all.

        zipped = np.broadcast(x, y, s)
        patches = [Circle((x_, y_), s_)
                for x_, y_, s_ in zipped]
        collection = PatchCollection(patches, **kwargs)
        if c is not None:
            c = np.broadcast_to(c, zipped.shape).ravel()
            collection.set_array(c)
            collection.set_clim(vmin, vmax)

        # ax = plt.gca()
        # ax.add_collection(collection)
        # ax.autoscale_view()
        self.ax0.add_collection(collection)
        self.ax0.autoscale_view()
        # plt.draw_if_interactive()
        if c is not None:
            # plt.sci(collection)
            self.ax0.sci(collection)
        # return collection

    #-----------------------------------------
    def bernstein_poly(self, i, n, t):
        """
        The Bernstein polynomial of n, i as a function of t
        """

        return comb(n, i) * ( t**(n-i) ) * (1 - t)**i


    def bezier_curve(self, points, nTimes=10):
        """
        Given a set of control points, return the
        bezier curve defined by the control points.

        points should be a list of lists, or list of tuples
        such as [ [1,1], 
                    [2,3], 
                    [4,5], ..[Xn, Yn] ]
            nTimes is the number of time steps, defaults to 1000

            See http://processingjs.nihongoresources.com/bezierinfo/
        """

        nPoints = len(points)
        xPoints = np.array([p[0] for p in points])
        yPoints = np.array([p[1] for p in points])

        t = np.linspace(0.0, 1.0, nTimes)

        polynomial_array = np.array([ self.bernstein_poly(i, nPoints-1, t) for i in range(0, nPoints)   ])

        xvals = np.dot(xPoints, polynomial_array)
        yvals = np.dot(yPoints, polynomial_array)

        return xvals, yvals

    #---------------------------------------------------------------------------
    def button_press(self, event):
        xval = event.xdata
        yval = event.ydata
        print("xval,yval=", xval,yval)
        self.bezpts[self.num_pts] = [xval,yval]
        print("self.bezpts=",self.bezpts)
        self.num_pts += 1

        rval = 8.0
        if self.num_pts == 4:
            rval = 12.0
            self.circles(xval,yval, s=rval, c='r', edgecolor='red')
            print("------- plot Bezier!")
            # xv, yv= bezier_curve(bezpts, nTimes=num_eval)  # eval Bezier
            # xv, yv= self.bezier_curve(self.bezpts, nTimes=self.num_eval)  # eval Bezier
            xv, yv= self.bezier_curve(self.bezpts,40)  # eval Bezier
            print("xv=",xv)
            print("yv=",yv)
            print("len(xv)=",len(xv))
            # self.csv_array = np.empty([1,4])  # should probably *just* np.delete, but meh
            # self.csv_array = np.delete(self.csv_array,0,0)
            for idx in range(len(xv)):
                print("idx=",idx)
            self.num_pts = 0

            # self.population_plot[self.discrete_scalar].ax0.plot(xv, yv, label=ctname, linewidth=lw, color=ctcolor)
            # self.ax0.plot(xv, yv, 'o')
            self.circles(xv,yv, s=rval, c='r', edgecolor='red')
            # self.ax0.set_xlabel('time (mins)')
            # self.ax0.set_ylabel('# of cells')
            # self.ax0.set_title("cell_type", fontsize=10)
            # self.ax0.legend(loc='center right', prop={'size': 8})
            self.canvas.update()
            self.canvas.draw()
            # self.population_plot[self.discrete_scalar].ax0.legend(loc='center right', prop={'size': 8})
            # self.show()
            return

        # self.circles(xvals,yvals, s=rvals, color=self.color_by_celltype[cell_type_index], alpha=self.alpha_value)
        self.circles(xval,yval, s=rval, c='k')

        self.ax0.set_aspect(1.0)

        self.ax0.set_xlim(self.plot_xmin, self.plot_xmax)
        self.ax0.set_ylim(self.plot_ymin, self.plot_ymax)

        # self.update_plots()
        self.canvas.update()
        self.canvas.draw()

        # print("event.inaxes", event.inaxes)
        # print("x", event.x)
        # print("y", event.y)

    # def on_release(self, event):
    #     print("release:")
    #     print("event.xdata", event.xdata)
    #     print("event.ydata", event.ydata)
    #     print("event.inaxes", event.inaxes)
    #     print("x", event.x)
    #     print("y", event.y)

    # def on_move(self, event):
    #     print("move")
    #     print("event.xdata", event.xdata)
    #     print("event.ydata", event.ydata)
    #     print("event.inaxes", event.inaxes)
    #     print("x", event.x)
    #     print("y", event.y)


if __name__ == "__main__":
    import sys

    app = QtWidgets.QApplication(sys.argv)
    w = MplWidget()
    w.show()
    sys.exit(app.exec_())
