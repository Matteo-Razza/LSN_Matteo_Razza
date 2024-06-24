import matplotlib.pyplot as plt
import numpy as np
import imageio

# Print initialization message
print("Generating images...")

# Load generation and best length data
gen, best = np.loadtxt("circ_fastest.dat", usecols=(0, 1), delimiter=' ', unpack=True)

# Load initial city coordinates
X, Y = np.loadtxt("circ_map.dat", usecols=(1, 2), delimiter=' ', unpack=True)

# Define the number of images to generate
tot = 500

# Prepare to store filenames
filenames = []

# Loop through each generation and generate the corresponding image
for j in range(tot):
    X, Y = np.loadtxt(f"Results/Circumference/path_{j+1}", usecols=(1, 2), delimiter=' ', unpack=True)
    
    fig = plt.figure(figsize=(8, 8))
    plt.plot(X, Y, label=f"Generation {j+1}, Length = {best[j]}")
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title('Best path among 34 cities randomly placed on a circumference')
    plt.ylim([-1.1, 1.1])
    plt.xlim([-1.1, 1.1])
    plt.grid(True)
    plt.legend(loc='upper right')
    
    # Update progress
    percent = (j + 1) * 100 / tot
    print(f"\rExecution progression: {percent:.2f}%", end="")
    
    filename = f"Images/Circumference/{j+1}.png"
    plt.savefig(filename, format="png", dpi=150)
    filenames.append(filename)
    plt.close(fig)

# Create a GIF from the saved images
with imageio.get_writer('circumference.gif', mode='I') as writer:
    print("\nGenerating GIF...")
    for i, filename in enumerate(filenames):
        image = imageio.imread(filename)
        writer.append_data(image)
        
        # Update progress
        percent = (i + 1) * 100 / len(filenames)
        print(f"\rWriting GIF frames: {percent:.2f}%", end="")

print("\nGif successfully created!")