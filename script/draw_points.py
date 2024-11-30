import matplotlib.pyplot as plt
import sys


def read_points(file_path):
    points = []
    with open(file_path, "r") as file:
        # Read the first line to get the number of points and dimensions
        _, dimensions = map(int, file.readline().split())

        if dimensions != 2:
            print("This script only supports 2D points.")
            sys.exit(1)

        # Read the remaining lines to get the points
        for line in file:
            point = tuple(map(int, line.split()))
            points.append(point)
    return points


def plot_points(points):
    # Separate the points into x and y coordinates
    x_coords, y_coords = zip(*points)

    # Plot the points
    plt.figure(figsize=(10, 6))
    plt.scatter(x_coords, y_coords, s=1)  # Use smaller marker size for large datasets
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.title("Scatter Plot of Points")
    plt.grid(True)
    plt.show()
    plt.savefig("n1000_10000.png")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python draw_points.py <file_path>")
        sys.exit(1)

    file_path = sys.argv[1]
    points = read_points(file_path)
    plot_points(points)
