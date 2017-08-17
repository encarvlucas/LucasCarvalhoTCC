import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

fig, ax = plt.subplots()
# fig.set_tight_layout(True)

# Query the figure's on-screen size and DPI. Note that when saving the figure to
# a file, we need to provide a DPI for that separately.
print('fig size: {0} DPI, size in inches {1}'.format(
    fig.get_dpi(), fig.get_size_inches()))

# Plot a scatter that persists (isn't redrawn) and the initial line.
x = np.arange(0, 20, 1)
y = []
for i in range(10):
    y.append(np.sin(x/np.pi + i))
y = np.array(y)
print y
# ax.scatter(x, x + np.random.normal(0, 3.0, len(x)))
# line, = ax.plot(x, x - 5, 'r-', linewidth=2)

def update(i):
    # label = 'timestep {0}'.format(i)
    # print(label)
    # Update the line and the axes (with a new xlabel). Return a tuple of
    # "artists" that have to be redrawn for this frame.
    # plt.plot(x,np.sin(x/np.pi))
    plt.plot(x,i)
    # ax.set_xlabel(label)
    return ax

# if __name__ == '__main__':
    # FuncAnimation will call the 'update' function for each frame; here
    # animating over 10 frames, with an interval of 200ms between frames.
anim = FuncAnimation(fig, update, frames=y, interval=100)
# if len(sys.argv) > 1 and sys.argv[1] == 'save':
#     anim.save('line.gif', dpi=80, writer='imagemagick')
# else:
    # plt.show() will just loop the animation forever.
plt.show()