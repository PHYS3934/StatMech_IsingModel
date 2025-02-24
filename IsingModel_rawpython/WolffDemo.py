# Import neceesary libraries:
import numpy as np
from numpy.random import rand
from scipy.stats import pearsonr
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import matplotlib
import tkinter as tk
from tkinter import ttk
matplotlib.use('TkAgg')  
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

# Import neceesary functions:
from myNeighbors import myNeighbors
from WolffIteration import WolffIteration

class WolffDemo:
    def __init__(self, N=100):
        # Parameters
        self.N = N
        self.J = 1
        self.kT = 2/np.log(1 + np.sqrt(2)) * 1.0  # T_c
        self.p = 1 - np.exp(-2*self.J/self.kT)

        # Initialize grids with -1 and 1
        self.grid0 = np.sign(0.5 - rand(N, N))
        self.grid = self.grid0.copy()
        
        # Get adjacency information
        self.adj = myNeighbors(np.arange(N*N), N)
        
        # Initialize state tracking
        self.correlations = [pearsonr(self.grid.flatten(), self.grid0.flatten())[0]]
        self.state = 'initial'
        self.iteration = 0
        
        # Create main window
        self.root = tk.Tk()
        self.root.title("Wolff Algorithm Demo")

        # Bind the close window event
        self.root.protocol("WM_DELETE_WINDOW", self.on_closing)
        
        # Create custom style for red buttons
        self.style = ttk.Style()
        self.style.configure('Red.TButton', 
                           background='red',
                           foreground='white',
                           font=('Helvetica', 15),
                           padding=(20, 10))
        
        # Make buttons larger when hovered
        self.style.map('Red.TButton',
                      background=[('active', 'darkred')])
        
        # Initialize autoplay state
        self.playing = False
        self.min_speed = 10  # minimum delay (fastest)
        self.max_speed = 100  # maximum delay (slowest)
        self.play_speed = 500  # initial delay in milliseconds
        
        # Create main container frame
        self.main_frame = ttk.Frame(self.root)
        self.main_frame.pack(fill=tk.BOTH, expand=True)

        # Create control frame at the top
        self.control_frame = ttk.Frame(self.main_frame)
        self.control_frame.pack(side=tk.TOP, fill=tk.X, padx=10, pady=5)
        
        # Create buttons in the control frame
        self.next_button = ttk.Button(self.control_frame, text='Next Step', command=self.step, style='Red.TButton')
        self.next_button.pack(side=tk.LEFT, padx=10, pady=10)
        
        self.play_button = ttk.Button(self.control_frame, text='Play', command=self.toggle_play, style='Red.TButton')
        self.play_button.pack(side=tk.LEFT, padx=10, pady=10)
        
        self.reset_button = ttk.Button(self.control_frame, text='Reset', command=self.reset, style='Red.TButton')
        self.reset_button.pack(side=tk.LEFT, padx=10, pady=10)
        
        # Create speed control
        ttk.Label(self.control_frame, text="Speed:").pack(side=tk.LEFT, padx=5)
        
        self.speed_scale = ttk.Scale(self.control_frame, 
                                   from_=0, to=100, 
                                   orient=tk.HORIZONTAL,
                                   length=100,
                                   command=self.update_speed)
        self.speed_scale.set(50)  # Set initial value
        self.speed_scale.pack(side=tk.LEFT, padx=5)
        
        # Create a single state label that gets updated with both pieces of info
        self.state_label = ttk.Label(self.control_frame, text='Initial State    Correlation: 1.000')
        self.state_label.pack(side=tk.LEFT, padx=15)
        
        self.status_label = ttk.Label(self.control_frame, text='Initial State')
        self.status_label.pack(side=tk.LEFT, padx=15, pady=5)
        
        self.corr_label = ttk.Label(self.control_frame, text=f'Correlation: {self.correlations[-1]:.3f}')
        self.corr_label.pack(side=tk.LEFT, padx=5, pady=5)

        # Create figure
        self.fig = plt.figure(figsize=(8, 10))
        self.gs = self.fig.add_gridspec(5, 1, height_ratios=[4, 4, 4, 4, 2])  # Changed last value from 1 to 2
        
        # Create canvas and pack it below the controls
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.main_frame)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True, padx=10, pady=5)
        
        # Initial display
        self.grid_display = np.where(self.grid == -1, 0, 1)
        self.update_display()
        
        print(f'Iteratively add neighbor to cluster with probability {self.p:.3f}')

    def update_display(self):
        """Update the plot display"""
        self.fig.clear()
        
        # Main plot
        ax1 = self.fig.add_subplot(self.gs[0:4])
        
        # Create custom colormap for the grid
        colors = ['white', 'black', '#E6550D', '#FEE6CE']
        custom_cmap = plt.cm.colors.ListedColormap(colors)
        
        # Display grid
        im = ax1.imshow(self.grid_display, cmap=custom_cmap, aspect='equal', vmin=0, vmax=4)
        cbar = plt.colorbar(im, ax=ax1, ticks=[0.5, 1.5, 2.5, 3.5], fraction=0.046, pad=0.04)
        cbar.set_ticklabels(['down', 'up', 'seed', 'cluster'])
        ax1.set_xticks([])
        ax1.set_yticks([])
        
        # Set title based on state
        if self.state == 'seed':
            ax1.set_title('SEED')
        elif self.state == 'cluster':
            ax1.set_title('CLUSTER')
        else:
            ax1.set_title(f'Grid State - Iteration {self.iteration}')
        
        # Correlation plot
        ax2 = self.fig.add_subplot(self.gs[4])
        ax2.plot(self.correlations, '.-k', linewidth=2)  
        ax2.set_xlabel('Iteration')
        ax2.set_ylabel('Corr0')
        ax2.set_yticks([0, 0.5, 1])
        ax2.set_ylim(0, 1)  
        ax2.grid(False) 
        ax2.tick_params(top=True, right=True)
        ax2.set_title(f'{self.correlations[-1]:.3f}')
                
        plt.subplots_adjust(hspace=0.2)  # Reduce space between plots
        self.canvas.draw()

    def step(self):
        """Advance the simulation one step"""
        if self.state == 'initial':
            # Get cluster update
            self.C, self.seed = WolffIteration(self.N, self.p, self.grid, self.adj)
            self.theNewSpin = -self.grid[np.unravel_index(self.seed-1, (self.N,self.N))]
            
            # Show seed
            self.grid_display = np.where(self.grid == -1, 0, 1)
            self.grid_display[np.unravel_index(self.seed-1, (self.N,self.N))] = 2
            self.status_label.configure(text='Seed Selected')
            self.state = 'seed'
            
        elif self.state == 'seed':
            # Show cluster
            self.grid_display = np.where(self.grid == -1, 0, 1)
            for c in self.C:
                self.grid_display[np.unravel_index(c-1, (self.N,self.N))] = 3
            self.grid_display[np.unravel_index(self.seed-1, (self.N,self.N))] = 2
            self.status_label.configure(text='Cluster Found')
            self.state = 'cluster'
            
        elif self.state == 'cluster':
            # Update spins
            for c in self.C:
                self.grid[np.unravel_index(c-1, (self.N,self.N))] = self.theNewSpin
            self.grid_display = np.where(self.grid == -1, 0, 1)
            
            # Update correlation
            self.correlations.append(
                pearsonr(self.grid.flatten(), self.grid0.flatten())[0]
            )
            self.corr_label.configure(text=f'Correlation: {self.correlations[-1]:.3f}')
            
            self.status_label.configure(text='Updated State')
            self.state = 'initial'
            self.iteration += 1
            
        self.update_display()

    def update_speed(self, value):
        """Update the animation speed based on slider value"""
        # Convert slider value (0-100) to speed (max_speed to min_speed)
        # We invert the scale so sliding right makes it faster
        normalized = float(value) / 100
        self.play_speed = int(self.max_speed - normalized * (self.max_speed - self.min_speed))
    
    def toggle_play(self):
        """Toggle autoplay state"""
        self.playing = not self.playing
        self.play_button.configure(text='Pause' if self.playing else 'Play')
        if self.playing:
            self.auto_step()
    
    def auto_step(self):
        """Automatically perform steps while playing"""
        if self.playing:
            self.step()
            self.root.after(self.play_speed, self.auto_step)
    
    def reset(self):
        """Reset the simulation to initial state"""
        self.playing = False
        self.play_button.configure(text='Play')
        self.grid = self.grid0.copy()
        self.correlations = [pearsonr(self.grid.flatten(), self.grid0.flatten())[0]]
        self.state = 'initial'
        self.iteration = 0
        self.grid_display = np.where(self.grid == -1, 0, 1)
        self.status_label.configure(text='Initial State')
        self.corr_label.configure(text=f'Correlation: {self.correlations[-1]:.3f}')
        self.update_display()

    def on_closing(self):
        """Handle window close event"""
        plt.close('all')  # Close all matplotlib windows
        self.root.destroy()  # Destroy the tkinter window
        self.root.quit()  # Stop the mainloop

    def run(self):
        """Start the main event loop"""
        self.root.mainloop()

if __name__ == "__main__":
    # Create and run the demo
    demo = WolffDemo()
    demo.run()
