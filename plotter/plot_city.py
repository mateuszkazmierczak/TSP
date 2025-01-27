import sys
import argparse
import plotly.graph_objs as go
from plotly.offline import plot
import os

def read_cities(file_path):
    """
    Reads city data from a file.

    :param file_path: Path to the cities file.
    :return: List of tuples [(id, x, y), ...]
    """
    cities = []
    with open(file_path, 'r') as f:
        lines = f.readlines()
        if not lines:
            print("City file is empty.")
            return cities
        try:
            total_cities = int(lines[0].strip())
        except ValueError:
            print("First line of city file should be the number of cities.")
            sys.exit(1)
        for line in lines[1:]:
            parts = line.strip().split()
            if len(parts) < 3:
                continue  # Skip incomplete lines
            try:
                city_id = int(parts[0])
                x = float(parts[1])
                y = float(parts[2])
                cities.append((city_id, x, y))
            except ValueError:
                print(f"Invalid line in city file: {line}")
                continue
    if len(cities) != total_cities:
        print(f"Warning: Expected {total_cities} cities, but found {len(cities)}.")
    return cities

def read_roads(file_path):
    """
    Reads road data from a file where the list of city IDs forms
    a path. Consecutive city IDs will be connected, and the last
    city will also connect to the first city (forming a cycle).

    :param file_path: Path to the roads file.
    :return: List of tuples [(city_id1, city_id2), ...]
    """
    roads = []
    with open(file_path, 'r') as f:
        tokens = f.read().split()

    # If there's fewer than 2 city IDs, we cannot form a road.
    if len(tokens) < 2:
        print("Not enough city IDs to form any roads.")
        return roads

    # Convert each token to int
    try:
        city_ids = list(map(int, tokens))
    except ValueError:
        print("All tokens in the road file must be integer city IDs.")
        sys.exit(1)

    # Connect every consecutive pair
    for i in range(len(city_ids) - 1):
        roads.append((city_ids[i], city_ids[i+1]))

    # Also connect the last city with the first (making a cycle)
    roads.append((city_ids[-1], city_ids[0]))

    return roads

def draw_cities_and_roads(cities, roads=None, output_filename="plot.html"):
    """
    Draws cities and roads using Plotly.

    :param cities: List of tuples [(id, x, y), ...]
    :param roads: List of tuples [(city_id1, city_id2), ...]
    :param output_filename: Filename for the output HTML plot.
    """
    if not cities:
        print("No cities to plot.")
        sys.exit(1)
    
    # Create a dictionary for quick lookup: city_id -> (x, y)
    city_dict = {city[0]: (city[1], city[2]) for city in cities}
    
    # Extract coordinates and IDs for city scatter
    ids = [city[0] for city in cities]
    xs = [city[1] for city in cities]
    ys = [city[2] for city in cities]
    
    # Create city scatter trace
    city_scatter = go.Scatter(
        x=xs,
        y=ys,
        mode='markers+text',
        text=[str(city_id) for city_id in ids],
        textposition='top right',
        marker=dict(color='blue', size=8),
        textfont=dict(color='red', size=10),
        name='Cities'
    )
    
    # Initialize figure
    fig = go.Figure()
    fig.add_trace(city_scatter)
    
    # Add roads if provided
    if roads:
        for idx, (city_id1, city_id2) in enumerate(roads, start=1):
            # Validate city IDs
            if city_id1 not in city_dict or city_id2 not in city_dict:
                print(f"Skipping road {idx}: Invalid city IDs {city_id1} or {city_id2}.")
                continue

            x1, y1 = city_dict[city_id1]
            x2, y2 = city_dict[city_id2]
            fig.add_trace(
                go.Scatter(
                    x=[x1, x2],
                    y=[y1, y2],
                    mode='lines',
                    line=dict(color='green', width=2),
                    showlegend=False,
                    name=f"Road {city_id1}-{city_id2}"
                )
            )
    
    # Update layout
    fig.update_layout(
        title="Visualization of Cities and Roads",
        xaxis_title="X Coordinate",
        yaxis_title="Y Coordinate",
        xaxis=dict(showgrid=True),
        yaxis=dict(showgrid=True),
        width=1000,
        height=800
    )
    
    # Save to HTML
    plot(fig, filename=output_filename, auto_open=False)
    print(f"Plot successfully saved to '{output_filename}'.")

def main():
    # Setup argument parser
    parser = argparse.ArgumentParser(
        description="Plot cities and optional roads using Plotly. Consecutive city IDs in roads form connections."
    )
    parser.add_argument('cities_file', help='Path to the cities data file.')
    parser.add_argument('-r', '--roads_file', help='Path to the roads data file.', default=None)
    parser.add_argument('-o', '--output', help='Output HTML filename.', default='plot.html')
    
    args = parser.parse_args()
    
    # Check if cities file exists
    if not os.path.isfile(args.cities_file):
        print(f"Cities file '{args.cities_file}' does not exist.")
        sys.exit(1)
    
    # Read cities
    cities = read_cities(args.cities_file)
    
    # Read roads if provided
    roads = None
    if args.roads_file:
        if not os.path.isfile(args.roads_file):
            print(f"Roads file '{args.roads_file}' does not exist.")
            sys.exit(1)
        roads = read_roads(args.roads_file)
    
    # Draw plot
    draw_cities_and_roads(cities, roads, args.output)

if __name__ == "__main__":
    main()
