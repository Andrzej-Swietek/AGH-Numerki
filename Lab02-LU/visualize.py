import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


Xa = 1/2
Xb = 2.0

output = pd.read_csv('data/output.txt', sep=' ', engine='python', header=None, names=['index', 'value'])
exact_values = pd.read_csv('data/exact_solution.txt', sep=' ', engine='python', header=None, names=['index', 'exact_value'])

mse = np.mean((output['value'] - exact_values['exact_value'])**2)
print("Mean Square Error:", mse)


# Plot data
plt.plot(output['index'], output['value'], label='Numerical Solution')
plt.plot(exact_values['index'], exact_values['exact_value'], label='Exact Solution')
plt.xlabel('Index')
plt.ylabel('Value')
plt.title('Comparison of Numerical and Exact Solutions')
plt.legend()
plt.grid(True)
plt.savefig('solution_plot.png')
plt.show()
