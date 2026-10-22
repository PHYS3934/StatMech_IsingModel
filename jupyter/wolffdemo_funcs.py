# Interactive demo of the Wolff cluster algorithm for the WolffDemo notebook:
# step through choosing a seed spin, growing its cluster, and flipping it,
# while tracking the correlation with the initial grid.

# Import necessary libraries
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from scipy.stats import pearsonr
import ipywidgets as widgets
from IPython.display import display, clear_output

# Import the Wolff cluster rule and the neighbor lookup (shared with IsingSim.ipynb)
from isingsim_funcs import WolffIteration, myNeighbors


class WolffDemo:
    def __init__(self, N=100):
        # Parameters
        self.N = N
        self.J = 1
        self.kT = 2/np.log(1 + np.sqrt(2)) * 1.0  # T_c
        # Probability of adding a parallel neighbor to the cluster (activating a bond):
        self.p = 1 - np.exp(-2*self.J/self.kT)

        # Random initial grid of -1 and 1 (kept as grid0 to compare against)
        self.grid0 = np.sign(0.5 - np.random.random((N, N)))
        self.grid = self.grid0.copy()

        # Linear indices of the four neighbors of each spin
        self.adj = myNeighbors(np.arange(N*N), N)

        # The demo cycles through three states:
        # 'initial' (choose a seed) -> 'seed' (grow its cluster) -> 'cluster' (flip it)
        self.state = 'initial'
        self.iteration = 0
        # Correlation with the initial grid shows how quickly the grid decorrelates:
        self.correlations = [1.0]

        # Create widgets
        self.next_button = widgets.Button(description='Next Step')
        self.next_button.on_click(self.step)

        self.reset_button = widgets.Button(description='Reset')
        self.reset_button.on_click(self.reset)

        self.status_label = widgets.Label(value='Initial State')
        self.corr_label = widgets.Label(value=f'Correlation with initial grid: {self.correlations[-1]:.3f}')

        # Create layout
        self.output = widgets.Output()
        self.container = widgets.VBox([
            widgets.VBox([self.status_label, self.corr_label]),
            self.output,
            widgets.HBox([self.next_button, self.reset_button])
        ])

        # Initial display (0 = down, 1 = up, 2 = seed, 3 = cluster)
        self.grid_display = np.where(self.grid == -1, 0, 1)
        self.update_display()
        display(self.container)
        print(f'Iteratively add neighbor to cluster with probability {self.p:.3f}')

    def update_display(self):
        """Update the plot display"""
        with self.output:
            clear_output(wait=True)

            fig = plt.figure(figsize=(8, 10))
            gs = plt.GridSpec(5, 1, height_ratios=[4, 4, 4, 4, 1], figure=fig)

            # Main plot: the grid, with the seed and cluster highlighted
            ax1 = fig.add_subplot(gs[0:4])
            custom_cmap = ListedColormap(['white', 'black', '#E6550D', '#FEE6CE'])
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

            # Correlation with the initial grid after each cluster flip
            ax2 = fig.add_subplot(gs[4])
            ax2.plot(self.correlations, '.-k')
            ax2.set_xlabel('Iteration')
            ax2.set_ylabel('Corr. with\ninitial grid')
            ax2.set_ylim(-0.1, 1)
            ax2.set_title(f'{self.correlations[-1]:.3f}')

            plt.subplots_adjust(hspace=-0.25)  # Reduce space between plots
            display(fig)
            plt.close(fig)

    def step(self, b=None):
        """Advance the demo one step"""
        if self.state == 'initial':
            # Grow a cluster (the display shows its seed first, then the whole cluster)
            self.cluster, self.seed = WolffIteration(self.N, self.p, self.grid, self.adj)

            # Show seed
            self.grid_display = np.where(self.grid == -1, 0, 1)
            self.grid_display.flat[self.seed] = 2
            self.status_label.value = 'Seed Selected'
            self.state = 'seed'

        elif self.state == 'seed':
            # Show cluster
            self.grid_display = np.where(self.grid == -1, 0, 1)
            self.grid_display.flat[self.cluster] = 3
            self.grid_display.flat[self.seed] = 2
            self.status_label.value = 'Cluster Identified'
            self.state = 'cluster'

        elif self.state == 'cluster':
            # Flip every spin in the cluster
            self.grid.flat[self.cluster] = -self.grid.flat[self.cluster]
            self.grid_display = np.where(self.grid == -1, 0, 1)

            # Update correlation with the initial grid
            self.correlations.append(pearsonr(self.grid.flatten(), self.grid0.flatten())[0])
            self.corr_label.value = f'Correlation with initial grid: {self.correlations[-1]:.3f}'

            self.status_label.value = 'Updated State'
            self.state = 'initial'
            self.iteration += 1

        self.update_display()

    def reset(self, b=None):
        """Reset the demo to the initial grid"""
        self.grid = self.grid0.copy()
        self.correlations = [1.0]
        self.state = 'initial'
        self.iteration = 0
        self.grid_display = np.where(self.grid == -1, 0, 1)
        self.status_label.value = 'Initial State'
        self.corr_label.value = f'Correlation with initial grid: {self.correlations[-1]:.3f}'
        self.update_display()
