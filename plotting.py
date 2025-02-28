import matplotlib.pyplot as plt
import numpy as np

# Function to read the txt file and extract vector data (valueX and valueY)
def read_vectors(file_path):
    vectors = []
    with open(file_path, 'r') as f:
        for line in f:
            parts = line.split()
            x = int(parts[0])
            y = int(parts[1])
            vector_x = float(parts[2])
            vector_y = float(parts[3])
            vectors.append((x, y, vector_x, vector_y))
    return vectors

# Function to read the txt file and extract strain data (strainXX, strainYY, strainXY)
def read_strain(file_path):
    strains = []
    with open(file_path, 'r') as f:
        for line in f:
            parts = line.split()
            x = int(parts[0])
            y = int(parts[1])
            strain_xx = float(parts[2])
            strain_yy = float(parts[3])
            strain_xy = float(parts[4])
            strains.append((x, y, strain_xx, strain_yy, strain_xy))
    return strains

# Function to plot the magnitude of vector components (valueX and valueY)
def plot_magnitudes(vectors):
    # Extract x, y, vector_x, and vector_y components for easy handling
    x_coords = [vec[0] for vec in vectors]
    y_coords = [vec[1] for vec in vectors]
    vector_x = [vec[2] for vec in vectors]
    vector_y = [vec[3] for vec in vectors]

    # Create the first plot for valueX magnitudes
    plt.subplot(2, 3, 1)
    plt.scatter(x_coords, y_coords, c=vector_x, cmap='viridis', marker='o')
    plt.colorbar(label='Motion X')
    plt.title('Motion X')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.gca().set_aspect('equal', adjustable='box')  # Ensures equal scaling of both axes

    # Create the second plot for valueY magnitudes
    plt.subplot(2, 3, 2)
    plt.scatter(x_coords, y_coords, c=vector_y, cmap='viridis', marker='o')
    plt.colorbar(label='Motion Y')
    plt.title('Motion Y')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.gca().set_aspect('equal', adjustable='box')  # Ensures equal scaling of both axes

# Function to plot strain components (strainXX, strainYY, strainXY)
def plot_strain(strains):
    # Extract x, y, strainXX, strainYY, strainXY components for easy handling
    x_coords = [strain[0] for strain in strains]
    y_coords = [strain[1] for strain in strains]
    strain_xx = [strain[2] for strain in strains]
    strain_yy = [strain[3] for strain in strains]
    strain_xy = [strain[4] for strain in strains]

    # Create the third plot for strainXX
    plt.subplot(2, 3, 3)
    plt.scatter(x_coords, y_coords, c=strain_xx, cmap='plasma', marker='s')
    plt.colorbar(label='strainXX')
    plt.title('strainXX')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.gca().set_aspect('equal', adjustable='box')  # Ensures equal scaling of both axes

    # Create the fourth plot for strainYY
    plt.subplot(2, 3, 4)
    plt.scatter(x_coords, y_coords, c=strain_yy, cmap='plasma', marker='s')
    plt.colorbar(label='strainYY')
    plt.title('strainYY')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.gca().set_aspect('equal', adjustable='box')  # Ensures equal scaling of both axes

    # Create the fifth plot for strainXY
    plt.subplot(2, 3, 5)
    plt.scatter(x_coords, y_coords, c=strain_xy, cmap='plasma', marker='s')
    plt.colorbar(label='strainXY')
    plt.title('strainXY')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.gca().set_aspect('equal', adjustable='box')  # Ensures equal scaling of both axes

# Example usage
vector_file_path = 'motion.txt'  # Path to your vector data file (valueX, valueY)
strain_file_path = 'strain.txt'   # Path to your strain data file (strainXX, strainYY, strainXY)

# Read data from both files
vectors = read_vectors(vector_file_path)
strains = read_strain(strain_file_path)

# Create a single figure for all plots
plt.figure(figsize=(15, 10))

# Plot the vector magnitudes and strain components as subplots
plot_magnitudes(vectors)
plot_strain(strains)

# Show the plot
plt.tight_layout()  # Automatically adjusts subplots to fit
plt.show()
