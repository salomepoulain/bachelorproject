def calculate_distance(point1, point2):
    x1, y1, z1 = point1
    x2, y2, z2 = point2
    dx = x2 - x1
    dy = y2 - y1
    dz = z2 - z1
    distance = (dx**2 + dy**2 + dz**2)**0.5
    return distance

point1 = (9.22,5.03,8.99)
point2 = (9.7,4.21,9.3)

print(calculate_distance(point1, point2))

point3 = (28.763, 20.257, 2.364)
point4 = (29.218, 19.471, 2.17)

print(calculate_distance(point3, point4))


# 1.8416568627190022
# 0.999449848666755