import matplotlib.pyplot as plt
import csv

# Initialize lists to store x and f(x)
x_vals = []
f_vals = []

try:
    # Read the CSV file
    with open('./sol.csv', mode='r') as file:
        csv_reader = csv.reader(file)
        next(csv_reader)  # Skip the header row
        for row in csv_reader:
            x_vals.append(float(row[0]))  # x value
            f_vals.append(float(row[1]))  # f(x) value
except FileNotFoundError:
    print(f"Error: The solution file (\"sol.csv\" in main directory) does not exist.")
    exit(1)  # Exit the program with an error code

# Plot the data
plt.plot(x_vals, f_vals, linestyle='-', color='b', label="f(x)")

# Add labels and title
plt.xlabel('x')
plt.ylabel('f(x)')
plt.title('Plot of f(x) vs. x')

# Show grid
plt.grid(True)

# Show the legend
plt.legend()

# Display the plot
plt.show()
