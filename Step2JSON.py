import sys

FREECAD_LIB_PATH_BIN = r"C:\Program Files\FreeCAD 1.0\bin"
FREECAD_LIB_PATH_LIB = r"C:\Program Files\FreeCAD 1.0\lib"

sys.path.append(FREECAD_LIB_PATH_BIN)
sys.path.append(FREECAD_LIB_PATH_LIB)

import FreeCAD
import Part
import os
import math
import json

def load_step_file(input_step_file):
    """
    Loads a STEP file.
    Args:
        input_step_file (str): Path to the STEP file to load.
    Returns:
        Part.Shape: The loaded shape.
    """
    shape = Part.Shape()
    shape.read(input_step_file)
    print(f"STEP file '{input_step_file}' loaded successfully.")
    return shape

def vector_to_list(vector):
    """
    Converts a FreeCAD Vector to a list.
    Args:
        vector (FreeCAD.Vector): The vector to convert.
    Returns:
        list: A list representation of the vector.
    """
    return [vector.x, vector.y, vector.z]

def extract_geometry_data(shape):
    """
    Extracts geometric and topological data from the shape.
    Args:
        shape (Part.Shape): The shape to extract data from.
    Returns:
        dict: A dictionary containing geometric and topological data.
    """

    geometry_data = {
        "solids": [],
        "shells": [],
        "faces": [],
        "edges": [],
        "vertices": [],
    }

    # Global mappings and counters
    vertex_indices = {}
    edge_indices = {}
    face_indices = {}
    vertex_counter = 1
    edge_counter = 1
    face_counter = 1

    def get_vertex_key(vertex):
        """
        Generates a unique key for a vertex based on its coordinates.
        Rounds coordinates to avoid floating-point precision issues.
        """
        x = round(vertex.Point.x, 6)
        y = round(vertex.Point.y, 6)
        z = round(vertex.Point.z, 6)
        return (x, y, z)

    def get_edge_key(edge):
        """
        Generates a unique key for an edge based on its start and end vertices.
        Ensures consistent ordering regardless of direction.
        """
        start_vertex_key = get_vertex_key(edge.Vertexes[0])
        end_vertex_key = get_vertex_key(edge.Vertexes[-1])

        # Ensure consistent ordering
        vertex_keys = sorted([start_vertex_key, end_vertex_key])
        return (vertex_keys[0], vertex_keys[1])

    # Extract Solids and their Faces
    # Process solids
    for solid_index, solid in enumerate(shape.Solids, start=1):
        solid_faces = []
        for face in solid.Faces:
            face_key = face.hashCode(1000000)  # Use face's unique hash code as key
            if face_key not in face_indices:
                face_index = face_counter
                face_counter += 1
                face_indices[face_key] = face_index

                surface = face.Surface

                face_edges = []
                face_vertices = []

                # Process edges
                for edge in face.Edges:
                    edge_key = get_edge_key(edge)
                    if edge_key not in edge_indices:
                        edge_index = edge_counter
                        edge_counter += 1
                        edge_indices[edge_key] = edge_index

                        curve = edge.Curve

                        # Process vertices
                        start_vertex = edge.Vertexes[0]
                        end_vertex = edge.Vertexes[-1]

                        start_vertex_key = get_vertex_key(start_vertex)
                        if start_vertex_key not in vertex_indices:
                            vertex_index = vertex_counter
                            vertex_counter += 1
                            vertex_indices[start_vertex_key] = vertex_index
                            geometry_data["vertices"].append({
                                "index": vertex_index,
                                "point": vector_to_list(start_vertex.Point),
                                "connected_edges": [],
                                "connected_faces": [],
                            })
                        else:
                            vertex_index = vertex_indices[start_vertex_key]
                        start_vertex_index = vertex_index

                        end_vertex_key = get_vertex_key(end_vertex)
                        if end_vertex_key not in vertex_indices:
                            vertex_index = vertex_counter
                            vertex_counter += 1
                            vertex_indices[end_vertex_key] = vertex_index
                            geometry_data["vertices"].append({
                                "index": vertex_index,
                                "point": vector_to_list(end_vertex.Point),
                                "connected_edges": [],
                                "connected_faces": [],
                            })
                        else:
                            vertex_index = vertex_indices[end_vertex_key]
                        end_vertex_index = vertex_index

                        # Edge data
                        edge_data = {
                            "index": edge_index,
                            "length": edge.Length,
                            "curve_type": type(curve).__name__,
                            "start_vertex": start_vertex_index,
                            "end_vertex": end_vertex_index,
                            "parameters": get_curve_parameters(curve, start_vertex.Point, end_vertex.Point),
                            "connected_faces": [face_index],
                        }
                        geometry_data["edges"].append(edge_data)

                        # Update vertex connectivity
                        geometry_data["vertices"][start_vertex_index - 1]["connected_edges"].append(edge_index)
                        geometry_data["vertices"][start_vertex_index - 1]["connected_faces"].append(face_index)
                        geometry_data["vertices"][end_vertex_index - 1]["connected_edges"].append(edge_index)
                        geometry_data["vertices"][end_vertex_index - 1]["connected_faces"].append(face_index)
                    else:
                        edge_index = edge_indices[edge_key]
                        # Update connected_faces
                        geometry_data["edges"][edge_index - 1]["connected_faces"].append(face_index)

                    face_edges.append(edge_index)

                # Face vertices (collect unique vertices for the face)
                for vertex in face.Vertexes:
                    vertex_key = get_vertex_key(vertex)
                    vertex_index = vertex_indices[vertex_key]
                    face_vertices.append(vertex_index)
                    # Ensure vertex is connected to this face
                    if face_index not in geometry_data["vertices"][vertex_index - 1]["connected_faces"]:
                        geometry_data["vertices"][vertex_index - 1]["connected_faces"].append(face_index)

                # Face data
                face_data = {
                    "index": face_index,
                    "solid_index": solid_index,
                    "area": face.Area,
                    "center_of_mass": vector_to_list(face.CenterOfMass),
                    "surface_type": type(surface).__name__,
                    "parameters": get_surface_parameters(surface),
                    "bounding_edges": face_edges,
                    "bounding_vertices": face_vertices,
                }
                geometry_data["faces"].append(face_data)
                solid_faces.append(face_index)
            else:
                # Face already processed
                face_index = face_indices[face_key]
                solid_faces.append(face_index)

        # Solid data
        solid_data = {
            "index": solid_index,
            "volume": solid.Volume,
            "center_of_mass": vector_to_list(solid.CenterOfMass),
            "faces": solid_faces,
        }
        geometry_data["solids"].append(solid_data)

    return geometry_data

def get_surface_parameters(surface):
    """
    Retrieves parameters of a surface based on its type.
    Args:
        surface (Part.GeomSurface): The surface to analyze.
    Returns:
        dict: Parameters of the surface.
    """
    params = {}
    if isinstance(surface, Part.Plane):
        params = {
            "type": "Plane",
            "position": vector_to_list(surface.Position),
            "normal": vector_to_list(surface.Axis),
        }
    elif isinstance(surface, Part.Cylinder):
        params = {
            "type": "Cylinder",
            "radius": surface.Radius,
            "axis": vector_to_list(surface.Axis),
            "center": vector_to_list(surface.Center),
        }
    elif isinstance(surface, Part.Cone):
        params = {
            "type": "Cone",
            "radius": surface.Radius,
            "semi_angle": math.degrees(surface.SemiAngle),
            "axis": vector_to_list(surface.Axis),
            "center": vector_to_list(surface.Center),
        }
    elif isinstance(surface, Part.Sphere):
        params = {
            "type": "Sphere",
            "radius": surface.Radius,
            "center": vector_to_list(surface.Center),
        }
    elif isinstance(surface, Part.Toroid):
        params = {
            "type": "Torus",
            "major_radius": surface.MajorRadius,
            "minor_radius": surface.MinorRadius,
            "axis": vector_to_list(surface.Axis),
            "center": vector_to_list(surface.Center),
        }
    elif isinstance(surface, Part.BSplineSurface):
        params = {
            "type": "BSplineSurface",
            "degree_u": surface.UDegree,
            "degree_v": surface.VDegree,
            "is_rational": surface.isRational(),
            "num_poles_u": len(surface.getPoles()),
            "num_knots_u": len(surface.getUKnots()),
            "num_poles_v": len(surface.getPoles()[0]),
            "num_knots_v": len(surface.getVKnots()),
        }
    else:
        params = {
            "type": "Unknown",
        }
    return params

def get_curve_parameters(curve, start_point, end_point):
    """
    Retrieves parameters of a curve based on its type.
    Args:
        curve (Part.GeomCurve): The curve to analyze.
        start_point (FreeCAD.Vector): The start point of the edge.
        end_point (FreeCAD.Vector): The end point of the edge.
    Returns:
        dict: Parameters of the curve.
    """
    params = {}
    if isinstance(curve, Part.Line):
        params = {
            "type": "Line",
            "start_point": vector_to_list(start_point),
            "end_point": vector_to_list(end_point),
            "direction": vector_to_list(curve.Direction),
        }
    elif isinstance(curve, Part.Circle):
        # Calculate start and end angles
        start_angle = math.degrees(curve.parameter(start_point))
        end_angle = math.degrees(curve.parameter(end_point))
        params = {
            "type": "Circle",
            "radius": curve.Radius,
            "center": vector_to_list(curve.Center),
            "axis": vector_to_list(curve.Axis),
            "start_angle": start_angle,
            "end_angle": end_angle,
        }
    elif isinstance(curve, Part.Ellipse):
        params = {
            "type": "Ellipse",
            "major_radius": curve.MajorRadius,
            "minor_radius": curve.MinorRadius,
            "center": vector_to_list(curve.Center),
            "axis": vector_to_list(curve.Axis),
            "start_point": vector_to_list(start_point),
            "end_point": vector_to_list(end_point),
        }
    elif isinstance(curve, Part.BSplineCurve):
        params = {
            "type": "BSplineCurve",
            "degree": curve.Degree,
            "is_rational": curve.isRational(),
            "num_poles": len(curve.getPoles()),
            "start_point": vector_to_list(start_point),
            "end_point": vector_to_list(end_point),
        }
    else:
        params = {
            "type": "Other",
            "curve_type": type(curve).__name__,
            "start_point": vector_to_list(start_point),
            "end_point": vector_to_list(end_point),
        }
    return params

def save_geometry_data_to_json(geometry_data, output_file):
    """
    Saves the extracted geometric data to a JSON file.
    Args:
        geometry_data (dict): The geometric data to save.
        output_file (str): The path to the output JSON file.
    """
    with open(output_file, 'w') as f:
        json.dump(geometry_data, f, indent=4)
    print(f"\nGeometry data saved to '{output_file}'.")

if __name__ == "__main__":
    # Specify the input and output files directly here
    input_step_file = r"C:\Users\ASUS\Downloads\flange.step"
    output_file = r"C:\Users\ASUS\Downloads\geometry_data.json"      

    if not os.path.isfile(input_step_file):
        print(f"Input STEP file '{input_step_file}' does not exist.")
        sys.exit(1)

    # Load STEP file
    shape = load_step_file(input_step_file)

    # Extract geometric and topological data
    geometry_data = extract_geometry_data(shape)

    # Save the extracted geometric data to a JSON file
    save_geometry_data_to_json(geometry_data, output_file)