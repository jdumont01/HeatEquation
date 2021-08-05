# Filename:  Plot2DHeat.py
# Purpose:   
#

# imports
import      sys
import      matplotlib.pyplot as plt
import      plotly as pl
import      numpy as np
from        mpl_toolkits.mplot3d import Axes3D
from        matplotlib.animation import FuncAnimation
import      matplotlib.animation as animation
from        matplotlib import cm
from        matplotlib.ticker import LinearLocator, FormatStrFormatter
import      matplotlib
from        IPython.display import HTML

class HeatOnRectangle():
    def __init__(self, xmax= 100, ymax= 100, tmax=10, dx = 0.1, dy= 0.1, dt=0.01, frames=500, bDebug = False):
        #Global debug
        self.bDebug = bDebug
        
        #Init physical grid for calculation
        self.dt = dt
        self.dx = dx
        self.dy = dy
        self.Tmax = tmax
        self.Xmax = xmax
        self.Ymax = ymax
        self.frames = frames
        
        self.t = np.linspace(0, tmax, 1)
        self.x = np.linspace(0, xmax, xmax+1)
        self.y = np.linspace(0, ymax, ymax+1)   
        self.u = np.zeros([len(self.x), len(self.y)], float)
        self.uAnimation = np.zeros([len(self.x), len(self.y), self.frames + 1], float)
        
        #init graphical output
        self.fig = plt.figure()                                 # Create figure object
        #self.fig.set_size_inches(6, 6)
        
        self.ax = Axes3D(self.fig, auto_add_to_figure=False)    # Create axes object
        self.fig.add_axes(self.ax)
        self.surf, = plt.plot([], [], [])                       # Create empty 3D line object
        plt.title("Heat distribution on a plane.")
        self.ax.set_xlabel("x (cm)")
        self.ax.set_ylabel("y (cm)")
        self.ax.set_zlabel("Temp (C) u(x,y,t)")      
        self.X_ax, self.Y_ax = np.meshgrid(self.x,self.y)
        
    # end of __init__
    
    def __str__(self):
        
        print ("u = %d".format(u))    
    # end of __str__

    def initSurface(self, x, y):
        '''
            x and y are lists (vectors) that define the length and width of the rectangle
            
            return an array initialized with zero's
        '''
        self.u = np.zeros([len(x), len(y)], float)

        return 

    def applyInitialConditionsAnim(self):
    
        #print("applyInitialConditionsAnim")
        #find some grid element near the center
        a = np.int32(np.floor(len(self.x)/2))
        b = np.int32(np.floor(len(self.y)/2))
                
        self.uAnimation[a-10:a+10, b-5:b+20, 0] = 50000.0

        if (self.bDebug == True):
            print (self.uAnimation[a, b, 0])
            print (self.uAnimation)
        
        return 

    def applyBoundaryConditionsAnim(self, t):
        #along the x-axis
        self.uAnimation[:, 0, t] = 0.0
        
        #along the x-axis on the upper-side of the rectangle
        self.uAnimation[:, len(self.y)-1, t] = 0.0
        
        #along the y-axis
        self.uAnimation[0,:,t] = 10000.0
        
        #along the y-axis on the right-side of the rectangle
        self.uAnimation[len(self.x)-1, :, t] = 10000.0

        if (self.bDebug == True):
            print (self.uAnimation)
        
        return
        
    def calcTempNextTimeSegmentAnim(self, t):
        for j in range(1, len(self.x)-1):
            if (self.bDebug == True):
                print ("j = ", j)
            for l in  range(1, len(self.y)-1):
                if (self.bDebug == True):
                    print ("l = ", l)
                self.uAnimation[j,l,t] = (1.0 / 4.0) * (self.uAnimation[j+1,l,t-1] + 
                                            self.uAnimation[j-1,l,t-1] + 
                                            self.uAnimation[j,l+1,t-1] + 
                                            self.uAnimation[j,l-1,t-1]) 

            if (self.bDebug == True):
                print (self.uAnimation)

        return

    def calcDataAnimation(self, tMax):

        #Generate the data throughout the time interval
        for t in range(1, tMax):        
            self.applyBoundaryConditionsAnim(t)
            self.calcTempNextTimeSegmentAnim(t)
        
        return 
    
    def animatePlot(self, nFrame, plot):
        print("nFrame =", nFrame)
        self.applyBoundaryConditionsAnim(nFrame)
        self.calcTempNextTimeSegmentAnim(nFrame)
        
        #plot.remove()
        self.ax.clear()
        plot = self.ax.plot_surface(X, Y, self.uAnimation[:,:,nFrame], cmap=cm.coolwarm, linewidth=0, antialiased=False)
        plt.title("Heat distribution on a plane.")
        
        return plot,
    
    def calcDataForAnimation(self, nFrame, X, Y, plot):
    
        plot[0].remove()
        #self.ax.clear()
        plt.title("Heat distribution on a plane.")
        plot[0] = self.ax.plot_surface(X, Y, self.uAnimation[:,:,nFrame], cmap=cm.coolwarm, linewidth=0, antialiased=False)
        plt.title("Heat distribution on a plane.")
        
        return 

    def animateSurfacePlot(self, T, X, Y, plot):
        #print("In animateSurfacePlot.  T =", T)
        # call the animator.  blit=True means only re-draw the parts that have changed.
        anim = animation.FuncAnimation(self.fig, self.calcDataForAnimation, T, fargs=(X, Y, plot), interval=20, blit=False)
        #anim = animation.FuncAnimation(r.fig, r.animatePlot, T, fargs=(X,Y,plot), interval=20, blit=False)

        # save the animation as an mp4.  This requires ffmpeg or mencoder to be
        # installed.  The extra_args ensure that the x264 codec is used, so that
        # the video can be embedded in html5.  You may need to adjust this for
        # your system: for more information, see
        # http://matplotlib.sourceforge.net/api/animation_api.html
        anim.save('basic_animation.mp4', fps=30)

        HTML(anim.to_html5_video())
        #display(html)
        #plt.close()
        print("Leaving animateSurfacePlot")
        
        return
    
    def plotSurface(self):

        #print("In PlotSurface")
        X, Y = np.meshgrid(self.x,self.y)

        # Plot the surface.
        plot = [self.ax.plot_surface(self.X_ax, self.Y_ax, self.uAnimation[:,:,0], cmap=cm.coolwarm, linewidth=0, antialiased=False)]

        # Add a color bar which maps values to colors.
        cb = self.fig.colorbar(plot[0], shrink=0.5, aspect=5)
        cb.set_label("Temp (C)")
        
        return X, Y, plot,

# end of class HeatOnRectangle

 
# *****
# Python entry point
# *****
if __name__ == "__main__":
    '''
    Process to create an animiated graphic using FuncAnimation (from http://www.acme.byu.edu/wp-content/uploads/2018/09/Animation.pdf)
    1. Compute all data to be plotted.
    2. Explicitly define figure object.
    3. Define line objects to be altered dynamically.
    4. Create function to update line objects.
    5. Create FuncAnimation object.
    6. Display using plt.show().

    Approach from the following sources:
    https://stackoverflow.com/questions/45712099/updating-z-data-on-a-surface-plot-in-matplotlib-animation
    https://pythonmatplotlibtips.blogspot.com/2018/11/animation-3d-surface-plot-artistanimation-matplotlib.html
    '''
    s = -0.20
    T = 100
    N = 100
    L = 100
    Nx = 1
    Ny = 1
    Nt = 1
    dt = T/Nt
    dx = N/Nx
    dy = L/Ny    
    debug = False
    matplotlib.matplotlib_fname()

    #1. Compute all data to be plotted.
    #2. Explicitly define figure object.
    r = HeatOnRectangle(N, L, T, dx, dy, dt, T, debug)
    r.applyInitialConditionsAnim()
    r.calcDataAnimation(T)
    #print (r.uAnimation)

    X, Y, plot = r.plotSurface()
    r.animateSurfacePlot(T, X, Y, plot )
    
    #plt.show()

    print ("Done!")    
