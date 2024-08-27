import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib as mlp

mlp.rcParams["pdf.fonttype"]=42
plt.rcParams.update({'font.family':'arial'})

# Create the DataFrame with your data
data = {
    'Pan-cancer': [217, 798, 47352],
    'GTEx': [62, 961, 64283]
}
index_labels = ['PSG-PSG', 'PSG-Other', 'Other-Other']
df = pd.DataFrame(data, index=index_labels)

# Apply log transformation to the data (adding 1 to avoid log(0))
df_log = np.log10(df + 1)

# Set up the figure and 3D axis
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

# Define the position of bars
x_labels = df.columns.tolist()
y_labels = df.index.tolist()

x = np.arange(len(x_labels))
y = np.arange(len(y_labels))
x, y = np.meshgrid(x, y)

# Flatten the meshgrid arrays
x = x.flatten().astype(float)  # Convert x to float for proper subtraction
y = y.flatten()

# Reduce spacing between Pan-cancer and GTEx by shifting GTEx bars closer to Pan-cancer
x[x == 1] -= 0.6  # Shift GTEx bars to the left by 0.3 units

# Flatten the data array for the bar heights
z = np.zeros_like(x)  # z represents the base height of bars (which is 0)
dx = dy = 0.15  # Reduced bar width and depth for thinner bars
dz = df_log.values.flatten()  # height of bars

# Define colors based on the interaction types
colors = {
    'PSG-PSG': '#b2e2e2',  # Blue
    'PSG-Other': '#66c2a4',  # Orange
    'Other-Other': '#2ca25f'  # Green
}

# Assign colors to bars based on their interaction type
bar_colors = [colors[y_labels[i]] for i in y]

# Plot the bars
ax.bar3d(x, y, z, dx, dy, dz, color=bar_colors, shade=True)

# Add the actual values on top of each bar
for i in range(len(x)):
    ax.text(x[i], y[i], dz[i], f'{int(df.values.flatten()[i])}', color='black', ha='center', va='bottom')

# Set the ticks and labels
ax.set_xticks(np.arange(len(x_labels)))
ax.set_xticklabels(x_labels, rotation=45)
ax.set_yticks(np.arange(len(y_labels)))
ax.set_yticklabels(y_labels)

ax.set_xlabel('Datasets')
ax.set_ylabel('Interaction Types')
ax.set_zlabel('Log(Count)')

# Add a title
ax.set_title('3D Bar Chart of Log-transformed Coexpression Pairs')

# Add a legend using dummy plots (this avoids the empty array error)
for interaction_type, color in colors.items():
    ax.plot([], [], [], color=color, label=interaction_type, linestyle='None', marker='o')

ax.legend(title='Interaction Types')

# Save the plot as a file or display it
plt.tight_layout()
plt.savefig(r"D:\essential_gene_evo\new_figure_1\figure3b.pdf", dpi=300)