import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from scipy.stats import norm

def dummy_histogram_y(f, s, x, width=0.02):
    size1 = s**2
    size2 = 2*s*(1-s)
    size3 = (1-s)**2

    peak1 = (f + (1-f)*s)**2
    peak2 = (1-f)*s*(f + (1-f)*s)
    peak3 = (1-f)**2 * s**2

    y1 = size1 * norm.pdf(x, peak1, width)
    y2 = size2 * norm.pdf(x, peak2, width)
    y3 = size3 * norm.pdf(x, peak3, width)

    return y1 + y2 + y3

fig, ax = plt.subplots(figsize=(6,6))
line, = ax.plot([], [], color='k', lw=2)
ax.set_xlim(0, 1)
ax.set_ylim(0, 10)
ax.set_xlabel(r'$d$')
ax.set_ylabel('Probability')

#### F OVER TIME
# PARAMS
s = 0.8
x = np.linspace(0, 1, 1000)
f_values = np.linspace(0, 1, 50)

def animate(i):
    f = f_values[i]
    y_total = dummy_histogram_y(f, s, x)
    line.set_data(x, y_total)
    ax.set_title(f'f={f:.2f}, s={s}')
    return line,

anim = FuncAnimation(fig, animate, frames=len(f_values), interval=200)
plt.show()

# #### S OVER TIME
# # PARAMS
# f = 0.5
# x = np.linspace(0, 1, 1000)
# s_values = np.linspace(0, 1, 50)
# 
# def animate(i):
#     s = s_values[i]
#     y_total = dummy_histogram_y(f, s, x)
#     line.set_data(x, y_total)
#     ax.set_title(f'f={f}, s={s:.2f}')
#     return line,
#
# anim = FuncAnimation(fig, animate, frames=len(s_values), interval=200)
# plt.show()
