# Filename:  Plot2DHeat.py
# Purpose:   
#

# imports
import      sys
import      math as mth
import      matplotlib.pyplot as plt
import      matplotlib.animation as animation
import      matplotlib
import      timeit
import      numpy as np
from        matplotlib import cm
from        mpl_toolkits.mplot3d import Axes3D
from        IPython.display import HTML

class HeatOnRectangle():
    def __init__(self, xmax= 100, ymax= 100, tmax=10, dx = 1, dy= 1, dt=1, k = 4.0, s = 0.25, icTemp = 500, bcXAxisMinY = 0.0, bcXAxisMaxY = 0.0, bcYAxisMinX = 0.0, bcYAxisMaxX = 0.0, frames=50, bDebug = False):
    
        #Global debug
        self.bDebug = bDebug
        
        #Init physical grid for calculation
        self.dt = dt
        self.dx = dx
        self.dy = dy
        self.Tmax = tmax
        self.Xmax = xmax
        self.Ymax = ymax

        self.Tseg = mth.floor(tmax/dt)
        self.Xseg = mth.ceil(xmax/dx)
        self.Yseg = mth.ceil(ymax/dy)
        
        # Set physical parameters
        self.k = k                  # thermal diffusivity
        self.s = s                  # vonNeumann 
        
        # Boundary Conditions
        self.bcXAxisMinY = bcXAxisMinY
        self.bcXAxisMaxY = bcXAxisMaxY
        self.bcYAxisMinX = bcYAxisMinX 
        self.bcYAxisMaxX = bcYAxisMaxX         
        
        # Initial Condition
        self.icTemp = icTemp
        
        # for video
        self.frames = frames

        # initialize matrices
        self.uAnimation_t = np.zeros((self.Xseg, self.Yseg), float)
        self.uAnimation_t1 = np.zeros((self.Xseg, self.Yseg), float)
        self.x = np.arange(0, self.Xseg, 1)
        self.y = np.arange(0, self.Yseg, 1)
        
        self._initialize()
    # end of __init__
    
    def __repr__(self):
        '''
            This method is used to print the results of the function u.
        '''
        print ("u = {:.1f}".format(u))    
    # end of __repr__

    def _initialize(self):
        if self.bDebug == True:
            print ("_initialize")
            
        #init graphical output
        self.fig = plt.figure()                                 # Create figure object
        self.fig.set_size_inches(8, 10)
        
        self.ax = self.fig.add_subplot(projection='3d')
        #self.ax = Axes3D(self.fig, auto_add_to_figure=True)    # Create axes object
        self.fig.add_axes(self.ax)
        self.ax.set_xlabel("x (cm)")
        self.ax.set_ylabel("y (cm)")
        self.ax.set_zlabel("Temp (C) u(x,y,t)")      
    
        return 
    # end of _initialize    

    def _applyInitialConditionsAnim(self):
        '''
            applyInitialConditionsAnim:
            updates the temperature array with the temperature set at certain points (center) 
            at the first time segment.
            
            return - none
        '''
    
        if (self.bDebug == True):
            print("_applyInitialConditionsAnim")
        
        #find some grid element near the center
        a = np.int32(np.floor(len(self.x)/2))
        b = np.int32(np.floor(len(self.y)/2))
                
        self.uAnimation_t[a-1:a+1, b-1:b+1] = self.icTemp
        
        if (self.bDebug == True):
            print (f"uAnimation_t = ")
            print (self.uAnimation_t)
            print (f"uAnimation_t1 = ")
            print (self.uAnimation_t1)
        
        return 
    # end of _applyInitialConditionsAnim
    
    def _applyBoundaryConditionsAnim(self):
        '''
            applyBoundaryConditionsAnim:
            updates the temperature array with the temperature set at certain points along the 
            edge(s) of the surface at some time segment t.

            return:  none
        '''
        if debug == True:
            print(f"In _applyBoundaryConditionsAnim")

        #along the x-axis
        self.uAnimation_t[:, 0] = self.bcXAxisMinY        
        if (self.bDebug == True):
            print (f'applyBoundaryConditionsAnim:  bcXAxisMinY')
            print (self.uAnimation_t)
        
        #along the x-axis on the upper-side of the rectangle
        self.uAnimation_t[:, -1] = self.bcXAxisMaxY
        if (self.bDebug == True):
            print (f'applyBoundaryConditionsAnim:  bcXAxisMaxY')
            print (self.uAnimation_t)
        
        #along the y-axis
        self.uAnimation_t[0, 1:-2] = self.bcYAxisMinX
        if (self.bDebug == True):
            print (f'applyBoundaryConditionsAnim:  bcYAxisMinX')
            print (self.uAnimation_t)
        
        #along the y-axis on the right-side of the rectangle
        self.uAnimation_t[-1, -1:-2] = self.bcYAxisMaxX
        if (self.bDebug == True):
            print (f'applyBoundaryConditionsAnim:  bcYAxisMaxX')
            print (self.uAnimation_t)

        self.uAnimation_t1 = self.uAnimation_t.copy()
        
        if (self.bDebug == True):
            print (f"uAnimation_t = ")
            print (self.uAnimation_t)
            print (f"uAnimation_t1 = ")
            print (self.uAnimation_t1)
        
        return
    # end of _applyBoundaryConditionsAnim
    
    def _initPlot(self):
        '''
            calcTempNextTimeSegmentAnim2:
            updates the temperature array for each internal point of the serface
            at some time segment t.  The equilibrium solution needs to be included.

            Note - for now the equilibrium solution is supported in the x direction only.
            return:  none
        '''
        if (self.bDebug == True):
            print(f'initPlot')
            print (f'x = {self.x}')
            print (f'y = {self.y}')

        #self._applyInitialConditionsAnim()
        self._applyBoundaryConditionsAnim()        
        
        if (self.bDebug == True):
            print (f"uAnimation_t = ")
            print (self.uAnimation_t)
            print (f"uAnimation_t1 = ")
            print (self.uAnimation_t1)
        
        return
    # end if _initPlot
    
    def _calcTempNextTimeSegmentAnim2(self, frameN, X, Y, plot):
        '''
            calcTempNextTimeSegmentAnim2:
            updates the temperature array for each internal point of the serface
            at some time segment t.  The equilibrium solution needs to be included.

            Note - for now the equilibrium solution is supported in the x direction only.
            return:  none
        '''
        if (self.bDebug == True):
            print (f'In _calcTempNextTimeSegmentAnim2')
            print(f'...before mutation')
            print (f't = {frameN}')
            print (f'x = {self.x}')
            print (f'y = {self.y}')
            print(self.uAnimation_t[:, :])
            print(self.uAnimation_t1[:, :])

        self.uAnimation_t1[1:-1 , 1:-1] = \
                self.uAnimation_t[1:-1 , 1:-1] + \
                (self.s) * (self.uAnimation_t[2: , 1:-1] + \
                self.uAnimation_t[:-2  ,1:-1] + \
                self.uAnimation_t[1:-1 ,2: ] + \
                self.uAnimation_t[1:-1 , :-2] - \
                4.0 * self.uAnimation_t[1:-1 , 1:-1])
                
        self.uAnimation_t = self.uAnimation_t1.copy()
        
        if (self.bDebug == True):
            print(f'...after mutation at time {frameN * self.dt * 1000} ms')
            print (self.uAnimation_t[:, :])
            print (self.uAnimation_t1[:, :])

        self._updatePlotObject(frameN, X, Y, plot)

        return 
    # end of _calcTempNextTimeSegmentAnim2
    
    def _calcTempNextTimeSegmentAnim(self, frameN, X, Y, plot):
        '''
            calcTempNextTimeSegmentAnim:
            updates the temperature array for each internal point of the serface
            at some time segment t.  The equilibrium solution needs to be included.

            Note - for now the equilibrium solution is supported in the x direction only.
            return:  none
        '''
        if (self.bDebug == True):
            print (f'In _calcTempNextTimeSegmentAnim')
            print (f'...before mutation')
            print (f't = {frameN}')
            print (f'x = {self.x}')
            print (f'y = {self.y}')
            print (self.uAnimation_t[:, :])
            print (self.uAnimation_t1[:, :])

        # Slow method using nested for-loops
        for j in range(1, len(self.x)-1):
            for l in  range(1, len(self.y)-1):
                self.uAnimation_t1[j,l] = (self.s) * (self.uAnimation_t[j+1,l] + \
                                            self.uAnimation_t[j-1,l] + \
                                            self.uAnimation_t[j,l+1] + \
                                            self.uAnimation_t[j,l-1] - \
                                            4.0 * self.uAnimation_t[j,l]) + \
                                            self.uAnimation_t[j,l] 

        self.uAnimation_t = self.uAnimation_t1.copy()

        if (self.bDebug == True):
            print(f'...after mutation')
            print(self.uAnimation_t[:, :])
            print(self.uAnimation_t1[:, :])
        
        self._updatePlotObject(frameN, X, Y, plot)

        return
    # end of _calcTempNextTimeSegmentAnim
    
    def _updatePlotObject(self, frameN, X, Y, plot):
        '''
            _updatePlotObject:
            
            return:  none
        '''
        if debug == True:
            print(f"In _updatePlotObject")
    
        plot[0].remove()
        self.ax.set_title(f'Heat distribution on the plane\nTime segment {frameN}', fontsize = 15)
        plot[0] = self.ax.plot_surface(X, Y, self.uAnimation_t[:,:], cmap=cm.coolwarm, linewidth=0, antialiased=False)

        return 
    # end of _updatePlotObject
    
    def _animateSurfacePlot(self):
        '''
            animateSurfacePlot:
            
            return:  none
        '''
        if debug == True:
            print(f"In _animateSurfacePlot")
        
        # setup the plot parameters
        X, Y, plot = self._createPlotSurface()
        
        # call the animator.  blit=True means only re-draw the parts that have changed.
        anim = animation.FuncAnimation(self.fig, self._calcTempNextTimeSegmentAnim2, self.Tseg, fargs=(X, Y, plot), init_func = self._initPlot, blit=False, repeat = False)   
        #anim = animation.FuncAnimation(self.fig, self._calcTempNextTimeSegmentAnim, self.Tseg, fargs=(X, Y, plot), init_func = self._initPlot, blit=False)               
        
        # save the animation as an mp4.  This requires ffmpeg or mencoder to be
        # installed.  The extra_args ensure that the x264 codec is used, so that
        # the video can be embedded in html5.  You may need to adjust this for
        # your system: for more information, see
        # http://matplotlib.sourceforge.net/api/animation_api.html
        
        #anim.save('basic_animation.mp4', writer='ffmpeg', fps=30)

        #html = HTML(anim.to_html5_video())
        #display(html)
        #self._addColorBar(plot)
        plt.show()
        #plt.close()
        
        return
    # end of _animateSurfacePlot
    
    def _createPlotSurface(self):
        '''
            _createPlotSurface:
            Creates the resoruces to plot the results on a surface plot.

            return:  meshgrid, and the plot list
        '''
        if debug == True:
            print(f"In _createPlotSurface")
    
        X, Y = np.meshgrid(self.x, self.y)

        # create a surface plot object.
        plot = [self.ax.plot_surface(X, Y, self.uAnimation_t[:,:], cmap=cm.coolwarm, linewidth=0, antialiased=False)]

        #plot = self._addColorBar(plot)
        
        return X, Y, plot,
    # end of _createPlotSurface

    def _addColorBar(self, plot):
        if debug == True:
            print(f"In _addColorBar")
    
        # Add a color bar which maps values to colors.
        cb = self.fig.colorbar(plot[0], shrink=0.5, orientation='horizontal', label = 'Temperature', extend = 'both' )
        cb.set_label("Temp (C)", labelpad =-1)
        
        return plot
    # end of _addColorBar
    
# end of class HeatOnRectangle

def main(k = 4.0, T = 10, Length = 10, Height = 10, deltaX = 0.1, deltaY = 0.1, bcXMinY = 0.0, bcXMaxY = 0.0, bcYMinX = 0.0, bcYMaxX = 0.0, icTemp = 0.0, s = 0.25, nFrames = 100, debug = False, animate = False):

    if debug == True:
        print(f"In Main")

    Nx = mth.ceil(Length/dx)           # Number of x-axis segments
    Ny = mth.ceil(Height/dy)           # Number of y-axis segments
    dt = (s *(dx)**2)/(k)              # max size of time segments (s)
    nFrames = mth.floor(T/dt)          # max number of frames = number of time segments
    
    matplotlib.matplotlib_fname()

    #1. Compute all data to be plotted.
    #2. Explicitly define figure object.
    r = HeatOnRectangle(Length, Height, T, deltaX, deltaY, dt, k, s, icTemp, bcXAxisMinY, bcXAxisMaxY, bcYAxisMinX, bcYAxisMaxX, nFrames, debug)

    r._animateSurfacePlot()
    
    print ("Done!") 

    sys.exit("done")
    
    pass
    
# end of main()

# *****
# Python entry point
# *****
if __name__ == "__main__":
    '''
        Define the parameters for the simulation
    '''
    
    # physical input parameters
    k = 4.0                             # thermal diffusivity of steel (4 mm2/s)
    T = 20                             # Number of time frames and max time for simulation (s)
    L = 10                              # Length of x-axis side (mm)
    H = 10                              # Length of y-axis side (mm)
    dx = 1                              # size of x-axis segments (mm)
    dy = 1                              # size of y-axis segments (mm)
    s = 0.25                            # von Neumann max value allowed for dt/(k dx) to allow stability
                                        # this will remain a user-defined parameter so that the user can    
                                        # observe the differences in the output.
    frames = 100                        # frames       
    initCondTemp = 0.0                  # Initial condition around the center of the plane
    
    # Boundary Conditions
    bcXAxisMinY = 5000.0                # BC along the x-axis when y= 0
    bcXAxisMaxY = 4000.0                # BC along the x-axis when y= L
    bcYAxisMinX = 0.0                   # BC along the y-axis when x = 0
    bcYAxisMaxX = 0.0                   # BC along the y-axis when x = H
    
    debug = False
    animate = False

    main(k, T, L, H, dx, dy, bcXAxisMinY, bcXAxisMaxY, bcYAxisMinX, bcYAxisMaxX, initCondTemp, s, frames, debug, animate)