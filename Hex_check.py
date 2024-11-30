import json
import math
import sys
import argparse

def load_json_data(file_path=None, json_string=None):
    try:
        if file_path:
            with open(file_path, 'r') as f:
                return json.load(f)
        elif json_string:
            return json.loads(json_string)
        else:
            raise ValueError("Either file_path or json_string must be provided.")
    except (FileNotFoundError, json.JSONDecodeError) as e:
        print(f"Error: {str(e)}")
        sys.exit(1)

def distance(p1, p2):
    return math.sqrt(sum((p1[i] - p2[i])**2 for i in range(3)))

def angle_between(v1, v2):
    dot_prod = sum(a * b for a, b in zip(v1, v2))
    mag1 = math.sqrt(sum(a**2 for a in v1))
    mag2 = math.sqrt(sum(b**2 for b in v2))
    cos_angle = max(min(dot_prod / (mag1 * mag2), 1.0), -1.0)
    return math.degrees(math.acos(cos_angle))

def calculate_outer_diameter(face_vertices, n_sides):
    # Calculate vertex-to-vertex diameter (longest diagonal)
    vertex_to_vertex_diameter = max(
        distance(face_vertices[i], face_vertices[j]) 
        for i in range(len(face_vertices)) 
        for j in range(i + 1, len(face_vertices))
    )
    # Convert vertex-to-vertex diameter to flat-to-flat diameter
    flat_to_flat_diameter = vertex_to_vertex_diameter * math.cos(math.pi / n_sides)
    return flat_to_flat_diameter

def detect_polygonal_prism(json_data, tolerance=1e-6):
    
    solids = json_data.get("solids", [])
    if len(solids) != 1:
        return (False, {})

    solid = solids[0]
    faces = [face for face in json_data.get("faces", []) if face["solid_index"] == solid["index"]]

    # Sort faces by Z-coordinate of their center of mass
    faces_sorted_by_z = sorted(faces, key=lambda f: f["center_of_mass"][2])
    bottom_face, top_face = faces_sorted_by_z[0], faces_sorted_by_z[-1]

    vertices = {vertex["index"]: vertex["point"] for vertex in json_data.get("vertices", [])}

    def get_face_vertices(face):
        return [vertices[v_idx] for v_idx in face["bounding_vertices"]]

    # Get vertices for top and bottom faces
    top_vertices = get_face_vertices(top_face)
    bottom_vertices = get_face_vertices(bottom_face)

    # Detect number of sides
    n_sides = len(top_vertices)
    if n_sides < 3:
        return (False, {})

    # Validate regularity of the polygon
    side_lengths = [distance(top_vertices[i], top_vertices[(i + 1) % n_sides]) for i in range(n_sides)]
    if max(side_lengths) - min(side_lengths) > tolerance:
        return (False, {})

    # Calculate outer diameter (flat-to-flat)
    outer_diameter = calculate_outer_diameter(top_vertices, n_sides)
    thickness = abs(top_face["center_of_mass"][2] - bottom_face["center_of_mass"][2])

    parameters = {
        "type": f"{n_sides}-Sided Prism",
        "n_sides": n_sides,
        "outer_diameter": outer_diameter,
        "thickness": thickness
    }

    return (True, parameters)


def main():
    parser = argparse.ArgumentParser(description="Detect geometric shapes from JSON data.")
    parser.add_argument('file', type=str, help='Path to the JSON file containing the solid data.')

    json_data = load_json_data(file_path=r"C:\Users\ASUS\Downloads\geometry_data.json")

    is_hex_prism, hex_params = detect_polygonal_prism(json_data)
    if is_hex_prism:
        print("Detected: Hexagonal Prism")
        print(f"Outer Diameter(Across Flats): {hex_params['outer_diameter']:.6f} mm")
        print(f"Thickness: {hex_params['thickness']:.6f} mm")
    else:
        print("Not a Hexagonal Prism.")

if __name__ == "__main__":
    main()
