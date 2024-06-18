import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

result_df = df.read_csv(f"data/data.csv", sep=" ")
print(result_df.head())
plt.xlabel('i-ty przedział')
plt.ylabel('')
plt.grid(True)
plt.legend()
plt.gca().tick_params(width=1.5, length=5)
plt.title('Monte Carlo')

plt.hist()
plt.show()