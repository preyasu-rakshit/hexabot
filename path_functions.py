from math import sin, cos, tan, pi
import matplotlib.pyplot as plt

def multiply_matrix_vector(matrix, vector):
    if len(matrix) != 3 or len(matrix[0]) != 3:
        raise ValueError("Matrix should be 3x3")
    if len(vector) != 3:
        raise ValueError("Vector should have 3 elements")

    result = [0, 0, 0]

    for i in range(3):
        for j in range(3):
            result[i] += matrix[i][j] * vector[j]

    return result


def rotate(theta, vector):
    return multiply_matrix_vector([[cos(theta), -sin(theta), 0], [sin(theta), cos(theta), 0], [0, 0, 1]], vector)


def interpolate_3d_points(point1, point2, num_intermediate_points):
    if len(point1) != 3 or len(point2) != 3:
        raise ValueError("Both input points should be 3D points (have 3 coordinates)")

    intermediate_points = []
    for i in range(num_intermediate_points + 2):
        t = i / (num_intermediate_points + 1)  # Calculate interpolation parameter
        x = (1 - t) * point1[0] + t * point2[0]
        y = (1 - t) * point1[1] + t * point2[1]
        z = (1 - t) * point1[2] + t * point2[2]
        intermediate_points.append([x, y, z])

    return intermediate_points


def generate_3d_bezier_curve(start_point, end_point, control_point1, control_point2, num_points):
    if (
        len(start_point) != 3
        or len(end_point) != 3
        or len(control_point1) != 3
        or len(control_point2) != 3
    ):
        raise ValueError("All input points should be 3D points (have 3 coordinates)")

    if start_point[0] != end_point[0]:
        raise ValueError("X-coordinates of start and end points must be the same.")

    bezier_curve = []

    for t in range(num_points + 1):
        t_normalized = t / num_points
        x = start_point[0]
        y = (
            (1 - t_normalized) ** 3 * start_point[1]
            + 3 * (1 - t_normalized) ** 2 * t_normalized * control_point1[1]
            + 3 * (1 - t_normalized) * t_normalized ** 2 * control_point2[1]
            + t_normalized ** 3 * end_point[1]
        )
        z = (
            (1 - t_normalized) ** 3 * start_point[2]
            + 3 * (1 - t_normalized) ** 2 * t_normalized * control_point1[2]
            + 3 * (1 - t_normalized) * t_normalized ** 2 * control_point2[2]
            + t_normalized ** 3 * end_point[2]
        )
        bezier_curve.append([x, y, z])

    return bezier_curve


def rotate_trajectory(angle, trajectory):
    out = []
    for i in trajectory:
        out.append(rotate(angle, i))
    return out


def generate_gait_points(direction=0, magnitude=1, res=50):
    control_point1 = [0, -0.8*magnitude, 5*magnitude]
    control_point2 = [0, 0.8*magnitude, 5*magnitude]
    line = interpolate_3d_points([0, 1*magnitude, 0], [0, -1*magnitude, 0], res)
    curve = generate_3d_bezier_curve([0, -1*magnitude, 0], [0, 1*magnitude, 0], control_point1, control_point2, int(res/1.1))
    path = line + curve
    final_path = rotate_trajectory(direction, path)
    return final_path

def plot_3d_points(points):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    x = [point[0] for point in points]
    y = [point[1] for point in points]
    z = [point[2] for point in points]
    ax.plot(x, y, z, marker='o', linestyle='-')
    plt.xlabel('X')
    plt.ylabel('Y')
    ax.set_zlabel('Z')
    plt.title('3D Bezier Curve')
    plt.show()
    

if __name__ == "__main__":
    plot_3d_points(generate_gait_points())