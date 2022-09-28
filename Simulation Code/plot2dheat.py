# Filename:     Plot2DHeat.py
# Purpose:      Plot the temperature distribution for a planar surface over time.
#               This simulation allows the user to compare the total time to compute
#               the temperature distribution between using native numpy arrays and
#               for-loops.
#
# Author:       Joe Dumont
# Date:         4-May-2022
#

# imports
import      sys
import      matplotlib
import      math                    as mth
import      matplotlib.pyplot       as plt
import      matplotlib.animation    as animation
import      numpy                   as np
import      time                    as stopwatch
from        matplotlib              import cm
from        mpl_toolkits.mplot3d    import Axes3D
from        IPython.display         import HTML

class HeatOnRectangle():
    def __init__(self, xmax= 100, ymax= 100, tmax=10, dx = 1, dy= 1, k = 4.0, s = 0.25, icTemp = 500, bcXAxisMinY = 0.0, bcXAxisMaxY = 0.0, bAnimate = True, bForLoopMethod = False, bDebug = False):

        #Global debug
        self.bDebug = bDebug

        #Init physical grid for calculation
        self.dt = (s *(dx)**2)/(k)         # max size of time segments (s)
        self.dx = dx
        self.dy = dy
        self.Tmax = tmax
        self.Xmax = xmax
        self.Ymax = ymax

        self.Tseg = mth.floor(tmax/self.dt)
        self.Xseg = mth.ceil(xmax/dx)
        self.Yseg = mth.ceil(ymax/dy)

        # Set calculation method
        self.bForLoopMethod = bForLoopMethod

        # Set physical parameters
        self.k = k                         # thermal diffusivity
        self.s = s                         # vonNeumann 

        # Boundary Conditions
        self.bcXAxisMinY = bcXAxisMinY
        self.bcXAxisMaxY = bcXAxisMaxY

        # Initial Condition
        self.icTemp = icTemp

        # for video
        self.frames = mth.floor(tmax/self.dt)   # max number of frames = number of time segments
        self.bAnimate = bAnimate

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
        '''
            _initialize:  initializes variables
            
            ARGS
            Input:      none
            Return:     none
        '''
        if self.bDebug == True:
            print ("_initialize")

        #init graphical output
        self.fig = plt.figure()                                 # Create figure object
        self.fig.set_size_inches(8, 10)

        self.ax = self.fig.add_subplot(projection='3d')
        self.fig.add_axes(self.ax)
        self.ax.set_xlabel("x (cm)")
        self.ax.set_ylabel("y (cm)")
        self.ax.set_zlabel("Temp (C) u(x,y,t)")      

        return 
    # end of _initialize    

    def _applyInitialConditions(self):
        '''
            _applyInitialConditions:
            updates the temperature array with the temperature set at certain points (center) 
            at the first time segment.
            
            ARGS
            Input:      none
            Return:     none
        '''
    
        if (self.bDebug == True):
            print("_applyInitialConditions")

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
    # end of _applyInitialConditions

    def _applyBoundaryConditions(self):
        '''
            _applyBoundaryConditions:  Applied the values along the 4 edges to the temperature matrix.

            ARGS
            Input:      none    
            Return:     none
        '''
        if self.bDebug == True:
            print(f"In _applyBoundaryConditions")

        # along the x-axis, y = 0
        self.uAnimation_t[:, 0] = self.bcXAxisMinY        
        if (self.bDebug == True):
            print (f'applyBoundaryConditionsAnim:  bcXAxisMinY')
            print (self.uAnimation_t)

        # along the x-axis on the upper-side of the rectangle, y = H
        self.uAnimation_t[:, -1] = self.bcXAxisMaxY
        if (self.bDebug == True):
            print (f'applyBoundaryConditionsAnim:  bcXAxisMaxY')
            print (self.uAnimation_t)

        c = (self.bcXAxisMaxY - self.bcXAxisMinY)/self.Xmax

        # along the y-axis; x = 0 
        # Calcualted based on the BCs of the physical problem - 
        # see the Jupyter notebook example for non-homogeneous BCs
        #self.uAnimation_t[0, 1:-2] = self.bcYAxisMinX
        for j in range(1, len(self.y) - 1): 
            self.uAnimation_t[0, j] = c * (j * self.dy) + self.bcXAxisMinY
            self.uAnimation_t[-1, j] = c * (j * self.dy) + self.bcXAxisMinY

        if (self.bDebug == True):
            print (f'applyBoundaryConditionsAnim:  bcYAxisMinX')
            print (self.uAnimation_t)

        #along the y-axis on the right-side of the rectangle; x = L
        #self.uAnimation_t[-1, -1:-2] = self.bcYAxisMaxX
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
    # end of _applyBoundaryConditions
    
    def _initPlot(self):
        '''
            _initPlot:  Initializes the values of the matrix
            
            ARGS
            Input:      none
            Return:     none
        '''
        if (self.bDebug == True):
            print(f'initPlot')
            print (f'x = {self.x}')
            print (f'y = {self.y}')

        #self._applyInitialConditions()
        self._applyBoundaryConditions()        
        
        if (self.bDebug == True):
            print (f"uAnimation_t = ")
            print (self.uAnimation_t)
            print (f"uAnimation_t1 = ")
            print (self.uAnimation_t1)
        
        return
    # end if _initPlot

    def _calcTempNextTimeSegmentForLoop(self, frameN):
        '''
            _calcTempNextTimeSegmentForLoop: Updates the temperature array for each internal point of the surface
                    at some time segment t using for-loops.

            ARGS
            Input:      frameN:int - the number of the time segment the plot is currently displaying
            Return:     none
        '''
        if (self.bDebug == True):
            print (f'In _calcTempNextTimeSegmentForLoop')
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
        
        return
    # end of _calcTempNextTimeSegmentForLoop
    
    def _calcTempNextTimeSegmentOpt(self, frameN):
        '''
            _calcTempNextTimeSegmentOpt: Updates the temperature array for each internal point of the surface
                    at some time segment t using numpy's C-optimized vector functions.  The values at the boundary 
                    are added at the start of the simulation when the boundary conditions are processed.  This approach
                    keeps the calculation focused on the parts of the temperature matrix that change in time.

            ARGS
            Input:      frameN:int - the number of the time segment the plot is currently displaying
            Return:     none
        '''
        if (self.bDebug == True):
            print (f'In _calcTempNextTimeSegmentAndAnimateOpt')
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
    
        return
    # end of _calcTempNextTimeSegmentOpt
    
    def _calcTempNextTimeSegmentAndAnimateOpt(self, frameN, X, Y, plot):
        '''
            _calcTempNextTimeSegmentAndAnimateOpt:  Wrapper function to get the matrix results and then update 
                        the plot object.

            ARGS
            Input:      frameN:int - the number of the time segment the plot is currently displaying
                        X:vector - the X values
                        Y:vector - the Y values
                        plot:plot_surface object - the plot to be mutated
            Return:     none
        '''
            
        if self.bDebug == True:
            print (f'In _calcTempNextTimeSegmentAndAnimateOpt')
       
        self._calcTempNextTimeSegmentOpt(frameN)
        self._updatePlotObject(frameN, X, Y, plot)

        return 
    # end of _calcTempNextTimeSegmentAndAnimateOpt
    
    def _calcTempNextTimeSegmentAndAnimateForLoop(self, frameN, X, Y, plot):
        '''
            _calcTempNextTimeSegmentAndAnimateForLoop:  Wrapper function to get the updated matrix for the next time 
                        segment and then update the plot object.
            ARGS
            Input:      frameN:int - the number of the time segment the plot is currently displaying
                        X:vector - the X values
                        Y:vector - the Y values
                        plot:plot_surface object - the plot to be mutated
            Return:     none
        '''
        if self.bDebug == True:
            print (f'In _calcTempNextTimeSegmentAndAnimateForLoop')
            
        self._calcTempNextTimeSegmentForLoop(frameN)
        self._updatePlotObject(frameN, X, Y, plot)

        return
    # end of _calcTempNextTimeSegmentAnim
    
    def _updatePlotObject(self, frameN, X, Y, plot):
        '''
            _updatePlotObject:  Updates the title and the the plot object.
            
            ARGS
            Input:      frameN:int - the number of the time segment the plot is currently displaying
                        X:vector - the X values
                        Y:vector - the Y values
                        plot:plot_surface object - the plot to be mutated
            Return:     none
        '''
        if self.bDebug == True:
            print(f"In _updatePlotObject")
    
        plot[0].remove()

        self.ax.set_title(f'Heat distribution on {self.Xmax} cm x {self.Ymax} cm plane\nFrame {frameN}\n Time = {frameN * self.dt:.2f} sec', fontsize = 15)
        plot[0] = self.ax.plot_surface(X, Y, self.uAnimation_t[:,:], cmap=cm.coolwarm, linewidth=0, antialiased=False)

        return 
    # end of _updatePlotObject
    
    def _animateSurfacePlot(self):
        '''
            _animateSurfacePlot:  Creates an animation of the plot over time.
            
            ARGS
            Input:      none
            Return:     none
        '''
        if self.bDebug == True:
            print(f"In _animateSurfacePlot")

        # setup the plot parameters
        X, Y, plot = self._createPlotSurface()
        
        # start the timer after the resources are set
        startT = stopwatch.perf_counter()
        
        # call the animator.  blit=True means only re-draw the parts that have changed.
        # this will not change the performance calculation as FuncAnimation is a single call. 
        if self.bForLoopMethod == True:
            anim = animation.FuncAnimation(self.fig, self._calcTempNextTimeSegmentAndAnimateForLoop, self.Tseg + 1, fargs=(X, Y, plot), init_func = self._initPlot, blit=False)               
        else:
            anim = animation.FuncAnimation(self.fig, self._calcTempNextTimeSegmentAndAnimateOpt, self.Tseg + 1, fargs=(X, Y, plot), init_func = self._initPlot, blit=False, repeat = False)   

        # Stop timer and get the time diff - do not include the time it takes to display the animation.
        endT = stopwatch.perf_counter()
        deltaT = endT - startT

        print(f'Simulation completed...after mutation at time {self.frames * self.dt * 1000} ms')
        print(f'Time to run simulation with FuncAnimation:  {deltaT:.3E} ms')
        
        plt.show()
        

        # save the animation as an mp4.  This requires ffmpeg or mencoder to be
        # installed.  The extra_args ensure that the x264 codec is used, so that
        # the video can be embedded in html5.  You may need to adjust this for
        # your system: for more information, see
        # http://matplotlib.sourceforge.net/api/animation_api.html
        anim.save('basic_animation.mp4', fps=30)

        #html = HTML(anim.to_html5_video())
        #display(html)
        #plt.close()
        
        return
    # end of _animateSurfacePlot
    
    def _generateTemperatureResults(self):
        '''
            _generateTemperatureResults:  Run the simulation without animation.
            
            ARGS
            Input:      none
            Return:     none
        '''
        if self.bDebug == True:
            print(f"In _generateTemperatureResults")

        startT = stopwatch.perf_counter()
        
        self._initPlot()
        
        # to avoid adding a conditional step to the time calculation, the for loop
        # is du-plicated
        if self.bForLoopMethod == True:
            for t in range(1, self.frames + 1):
                self._calcTempNextTimeSegmentForLoop(t)
        else:
            for t in range(1, self.frames + 1):   
                self._calcTempNextTimeSegmentOpt(t)        

        endT = stopwatch.perf_counter()
        deltaT = endT - startT
        
        print(f'Simulation completed...after mutation at time {self.frames * self.dt * 1000} ms')
        print(f'Time to run simulation {deltaT:.3E} ms')
        print(self.uAnimation_t[:, :])
        
        return
    # end of _generateTemperatureResults
    
    def _createPlotSurface(self):
        '''
            _createPlotSurface:  Creates the resoruces to plot the results on a surface plot.

            ARGS
            Input:      none
            Return:     meshgrid, and the plot list
        '''
        if self.bDebug == True:
            print(f"In _createPlotSurface")
    
        X, Y = np.meshgrid(self.x, self.y)

        # create a surface plot object.
        plot = [self.ax.plot_surface(X, Y, self.uAnimation_t[:,:], cmap=cm.coolwarm, linewidth=0, antialiased=False)]

        #plot = self._addColorBar(plot)
        
        return X, Y, plot,
    # end of _createPlotSurface

    def _addColorBar(self, plot):
        '''
            _addColorBar:  Adds a colorbar to the plot.  Given that the data changes over time, the simulation must be
                        run to completion to save the minimum and maximum values for u(x,y,t).
                        
            ARGS
            Input:      plot:plot_surface - plot_surface object
            Return:     plot - updated plot surface object.
        '''
        if self.bDebug == True:
            print(f"In _addColorBar")
    
        # Add a color bar which maps values to colors.
        cb = self.fig.colorbar(plot[0], shrink=0.5, orientation='horizontal', label = 'Temperature', extend = 'both' )
        cb.set_label("Temp (C)", labelpad =-1)
        
        return plot
    # end of _addColorBar
    
# end of class HeatOnRectangle

def runSimulation(L, H, T, dx, dy, k, s, icTemp, bcXAxisMinY, bcXAxisMaxY, bAnimate, bForLoopMethod, bDebug):
    '''
        runSimulation (Public function):  Creates the HeatOnRectangle object and initiates the simulation.
        
        ARGS
        Input:      L:float - length (L) of the plane; x will range from [0 to L] 
                    H:float - height (H) of the plane; y will range from [0 to H]
                    T:float - maximum time (s); t will range from [0 to T] 
                    dx:float - size of the differential segment along the x-axis
                    dy:float - size of the differential segment along the y-axis
                    k:float - diffusion constant
                    s:float - numerical convergence constant (derived from von Neumann procedure)
                    icTemp:float - initial condition temperature
                    bcXAxisMinY:float - boundary condition temperature along the x-axis at y = 0
                    bcXAxisMaxY:float - boundary condition temperature along the x-axis at y = H
                    bAnimate:bool - True to use FuncAnimation; False otherwise
                    bForLoopMethod: bool -  True will use the for loop function to calculate the next time segment
                                            False will use the numpy optimized method
                    bDebug:bool - include debug statements in stdio 
        Return:     none
    '''
    if bDebug == True:
        print (f'In _runSimulation')

    r = HeatOnRectangle(L, H, T, dx, dy, k, s, icTemp, bcXAxisMinY, bcXAxisMaxY, bAnimate, bForLoopMethod, bDebug)
        
    if bAnimate == True:
        r._animateSurfacePlot()
    else:
        r._generateTemperatureResults()
        
    return
# end of _runSimulation
    
# -----------------------------------
def main(L = 10, H = 10, T = 10, dx = 0.1, dy = 0.1, bcXMinY = 0.0, bcXMaxY = 0.0, icTemp = 0.0, k = 4.0, s = 0.25, bAnimate = False, bForLoopMethod = False, bDebug = False):

    if debug == True:
        print(f"In Main")
    
    runSimulation(L, H, T, dx, dy, k, s, icTemp, bcXAxisMinY, bcXAxisMaxY, bAnimate, bForLoopMethod, bDebug)
    
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
    T = 50                              # Max time for simulation (s) - not runtime.
    L = 30                              # Length of x-axis side (mm)
    H = 30                              # Length of y-axis side (mm)
    dx = 1                              # size of x-axis segments (mm)
    dy = 1                              # size of y-axis segments (mm)
    k = 4.0                             # thermal diffusivity of steel (4 mm2/s)
    s = 0.25                            # von Neumann max value allowed for dt/(k dx) to allow stability
                                        # this will remain a user-defined parameter so that the user can    
                                        # observe the differences in the output.

    # Initial Conditions (predetermined spots on plane)
    initCondTemp = 0.0                  
    
    # Boundary Conditions
    bcXAxisMinY = 5000.0                # BC along the x-axis when y=0
    bcXAxisMaxY = 2000.0                # BC along the x-axis when y=H
    
    debug = False
    animate = True
    useForLoopMethod = False

    main(L, H, T, dx, dy, bcXAxisMinY, bcXAxisMaxY, initCondTemp, k, s, animate, useForLoopMethod, debug)