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

        self.instructions = f"`c` clear; `u` undo; `n` new type; `s` save"
        self.cell_type_name = "epi"
        self.cell_volume = 2494.0
        self.cell_radius = (self.cell_volume * 0.75 / np.pi) ** (1./3)
        print(f"default self.cell_radius= {self.cell_radius}")
        # self.cell_radius = 8.412710547954228   # from PhysiCell_phenotype.cpp
        self.color_by_celltype = ['gray','red','green','yellow','cyan','magenta','blue','brown','black','orange','seagreen','gold']

        self.canvas.axes = self.canvas.figure.add_subplot(111)
        self.ax0 = self.canvas.axes
        self.ax0.set_title(self.instructions, fontsize=10)

        self.num_pts = 0
        self.bezpts = np.zeros((4,2))
        self.bez_plot = []
        print("self.bezpts=",self.bezpts)
        self.num_eval = 10

        self.plot_xmin = -500
        self.plot_xmax = 500
        self.plot_ymin = -500
        self.plot_ymax = 500

        self.xv = None  # xvalues of points along Bezier curve
        self.yv = None  # yvalues of points along Bezier curve

        self.ax0.set_aspect(1.0)
        self.ax0.set_xlim(self.plot_xmin, self.plot_xmax)
        self.ax0.set_ylim(self.plot_ymin, self.plot_ymax)

        self.canvas.mpl_connect("button_press_event", self.button_press)
        # self.canvas.mpl_connect("button_release_event", self.on_release)
        # self.canvas.mpl_connect("motion_notify_event", self.on_move)

    def keyPressEvent(self, event):
        # if event.key() == Qt.Key_Space:
            # self.test_method()
        if event.key() == Qt.Key_C:
            self.ax0.cla()
            self.ax0.set_aspect(1.0)
            self.ax0.set_xlim(self.plot_xmin, self.plot_xmax)
            self.ax0.set_ylim(self.plot_ymin, self.plot_ymax)
            self.ax0.set_title(self.instructions, fontsize=10)
            self.canvas.update()
            self.canvas.draw()
        elif event.key() == Qt.Key_U:
            # print("----- Undo is TBD...")
            print("# plots= ",len(self.bez_plot))
            print(type(self.bez_plot))
            print(self.bez_plot)
            # self.bez_plot = self.bez_plot.pop()
            self.bez_plot.pop()
            # self.bez_plot.del[-1]
            print("# plots= ",len(self.bez_plot))
            self.update_plots()
        elif event.key() == Qt.Key_N:
            # print("----- New cell type ...")
            self.cell_type_name = input("Enter cell type name:    (disregard this--> ")
            print("-- new cell type name is ", self.cell_type_name)
            volume = input("Enter cell type volume:     (disregard this--> ")
            self.cell_volume = float(volume)
            self.cell_radius = (self.cell_volume * 0.75 / np.pi) ** (1./3)
            print(f"-- its volume is {self.cell_volume} (and r={self.cell_radius}")
        elif event.key() == Qt.Key_S:
            fname = "curvy.csv"
            print(f"----- Writing to {fname}...")
            with open(fname, 'w') as f:
                f.write('x,y,z,type,volume,cycle entry,custom:GFP,custom:sample\n')
                for idx in range(len(self.bez_plot)):
                    xvals = self.bez_plot[idx][:,0]
                    yvals = self.bez_plot[idx][:,1]
                    for jdx in range(len(xvals)):
                        f.write(f"{xvals[jdx]},{yvals[jdx]},0,{self.cell_type_name}\n")

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
    # def update_plots(self,xvals,yvals):
    def update_plots(self):
        print("----------------  update_plots -----------------")
        self.ax0.cla()

        # print(f"------- update_plots(): type(xvals)={type(xvals)}")
        # print(f"        type(xvals)={type(xvals)}")
        # print(f"        xvals.shape={xvals.shape}")
        # print(f"        self.bezpts={self.bezpts}")

        print(f"        len(self.bez_plot={len(self.bez_plot)}")
        rval = 12.0
        for idx in range(len(self.bez_plot)):
            xvals = self.bez_plot[idx][:,0]
            yvals = self.bez_plot[idx][:,1]
            print("xvals=", xvals)
            print("yvals=", yvals)
            self.ax0.plot(xvals, yvals)   # plot 4-pt control polygon

            # self.xv, self.yv= self.bezier_curve(self.bezpts,40)  # eval Bezier
            # self.xv, self.yv= self.bezier_curve(self.bez_plot[idx],40)  # eval Bezier
            xvals, yvals= self.bezier_curve(self.bez_plot[idx],40)  # eval Bezier

            # self.circles(self.xv,self.yv, s=rval, c='r', edgecolor='red')
            self.circles(xvals,yvals, s=rval, c='r', edgecolor='red')

            # self.xv2, self.yv2= self.cell_spacing(self.xv,self.yv)  # eval Bezier

            self.ax0.plot(xvals, yvals)   # plot 4-pt control polygon

        # self.ax0.set_xlabel('time (mins)')
        # self.ax0.set_ylabel('# of cells')
        # self.ax0.set_title("'c' to clear all.", fontsize=10)
        self.ax0.set_title(self.instructions, fontsize=10)

        self.ax0.set_aspect(1.0)
        self.ax0.set_xlim(self.plot_xmin, self.plot_xmax)
        self.ax0.set_ylim(self.plot_ymin, self.plot_ymax)
        # self.ax0.legend(loc='center right', prop={'size': 8})
        self.canvas.update()
        self.canvas.draw()
        return

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
    def cell_spacing(self, xv,yv):  # eval Bezier
        return xv, yv

    #---------------------------------------------------------------------------
    def button_press(self, event):
        print("----------------  button_press -----------------")
        xval = event.xdata
        yval = event.ydata
        if xval is None:
            return
        print("xval,yval=", xval,yval)
        self.bezpts[self.num_pts] = [xval,yval]
        print("self.bezpts=",self.bezpts)

        # xvals = self.bezpts[self.num_pts][0]
        # yvals = self.bezpts[self.num_pts][1]
        xvals = self.bezpts[:,0]
        yvals = self.bezpts[:,1]
        print("xvals=", xvals)
        print("yvals=", yvals)

        self.num_pts += 1

        # if self.num_pts > 1:
        #     self.ax0.plot(xvals, yvals)

        rval = 8.0
        if self.num_pts == 4:
            rval = 12.0
            self.circles(xval,yval, s=rval, c='r', edgecolor='red')
            print("------- plot Bezier!")
            # xv, yv= bezier_curve(bezpts, nTimes=num_eval)  # eval Bezier
            # xv, yv= self.bezier_curve(self.bezpts, nTimes=self.num_eval)  # eval Bezier
            # self.xv, self.yv= self.bezier_curve(self.bezpts,40)  # eval Bezier

            self.bez_plot.append(self.bezpts.copy())
            print(f"        type(self.bez_plot)={type(self.bez_plot)}")
            print(f"        len(self.bez_plot)={len(self.bez_plot)}")
            print(f"        self.bez_plot={self.bez_plot}")

            # print("xv=",xv)
            # print("yv=",yv)
            # print("len(xv)=",len(self.xv))
            # self.csv_array = np.empty([1,4])  # should probably *just* np.delete, but meh
            # self.csv_array = np.delete(self.csv_array,0,0)
            # for idx in range(len(xv)):
            #     print("idx=",idx)
            self.num_pts = 0

            # self.population_plot[self.discrete_scalar].ax0.plot(xv, yv, label=ctname, linewidth=lw, color=ctcolor)
            # self.ax0.plot(xv, yv, 'o')

            # self.circles(self.xv,self.yv, s=rval, c='r', edgecolor='red')
            # self.xv2, self.yv2= self.cell_spacing(self.xv,self.yv)  # eval Bezier

            # self.ax0.plot(xvals, yvals)   # plot 4-pt control polygon
            # # self.ax0.set_xlabel('time (mins)')
            # # self.ax0.set_ylabel('# of cells')
            # self.ax0.set_title("Press 'c' to clear all.", fontsize=10)
            # # self.ax0.legend(loc='center right', prop={'size': 8})
            # self.canvas.update()
            # self.canvas.draw()
            # # self.population_plot[self.discrete_scalar].ax0.legend(loc='center right', prop={'size': 8})
            # # self.show()

            # self.update_plots(xvals,yvals)
            self.update_plots()
            return

        # self.circles(xvals,yvals, s=rvals, color=self.color_by_celltype[cell_type_index], alpha=self.alpha_value)
        self.circles(xval,yval, s=rval, c='k')

        self.ax0.set_aspect(1.0)

        self.ax0.set_xlim(self.plot_xmin, self.plot_xmax)
        self.ax0.set_ylim(self.plot_ymin, self.plot_ymax)

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
