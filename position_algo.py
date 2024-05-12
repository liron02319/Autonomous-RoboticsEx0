import csv
import numpy as np
import os

# Constants
C = 299792458  # Speed of light (m/s)
OMEGA_E = 7.2921151467e-5  # Earth's rotation rate (rad/s)

# Read the CSV file
def read_csv(filename):
    data = []
    with open(filename, 'r', encoding='utf-8-sig') as file:
        reader = csv.reader(file)
        next(reader)
        for row in reader:
            data.append(row)
    return data

# Compute satellite positions
def compute_sat_pos(data, gps_time):
    sat_positions = []
    for row in data:
        if row[0] == gps_time:
            sat_x = float(row[2])
            sat_y = float(row[3])
            sat_z = float(row[4])
            sat_positions.append((sat_x, sat_y, sat_z))
    return sat_positions

# Iterative least-squares algorithm
def compute_position(data, gps_time):
    sat_positions = compute_sat_pos(data, gps_time)
    num_sats = len(sat_positions)

    # Check if there are enough satellites for positioning
    if num_sats < 4:
        print("Not enough satellites for positioning.")
        return None

    # Set up the least-squares problem
    X = np.zeros((num_sats, 4))
    Y = np.zeros(num_sats)

    for i, (sat_x, sat_y, sat_z) in enumerate(sat_positions):
        X[i, 0] = sat_x
        X[i, 1] = sat_y
        X[i, 2] = sat_z
        X[i, 3] = 1
        Y[i] = np.sqrt(sat_x**2 + sat_y**2 + sat_z**2)

    # Solve the least-squares problem using NumPy
    beta = np.linalg.lstsq(X, Y, rcond=None)[0]
    receiver_pos = beta[:3]

    return receiver_pos

# Compute positions for all GPS times
def compute_positions(data):
    gps_times = set([row[0] for row in data])
    positions = {}

    for gps_time in gps_times:
        receiver_pos = compute_position(data, gps_time)
        if receiver_pos is not None:
            positions[gps_time] = (receiver_pos)

    return positions

# Get the current working directory
cwd = os.getcwd()

# Construct the file path
file_path = os.path.join(cwd, 'plot_data.csv')

# Read the CSV file
data = read_csv(file_path)
positions = compute_positions(data)

for gps_time, position in positions.items():
    receiver_pos = position
    print(f"GPS Time: {gps_time}")
    print(f"Receiver Position (X, Y, Z): {receiver_pos}")
    print()
