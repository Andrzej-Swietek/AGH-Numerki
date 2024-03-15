import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


output = pd.read_csv('wyniki3.txt', sep=' ', engine='python', header=None, names=["x","y"])



# Plot data
plt.plot(output["x"], output['y'], label='Solution 3')
plt.xlabel('t')
plt.ylabel('x(t)')
plt.title('beta=0.4, F_0 = 0.1')
plt.legend()
plt.grid(True)
plt.savefig('solution_plot.png')
plt.show()

