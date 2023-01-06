import math
import numpy as np
import logging

logger = logging.getLogger(__name__)


def angle3Dchain(liste):
    logger.debug(f"liste {liste}")
    X0 = np.asarray(liste[-2])
    X1 = np.asarray(liste[-1])
    [x, y, z] = X1-X0

    logger.debug(f"X0:{X0}; X1:{X1}; [x,y,z]: {[x,y,z]}")

    length = (x*x + y*y + z*z)**0.5
    phi = math.asin(z/length)

    if y >= 0 and (x*x+y*y != 0):
        theta = math.acos(x/(x*x+y*y)**0.5)
    elif y >= 0 and (x*x+y*y == 0):
        theta = 0.0
    else:
        theta = 2*math.pi-math.acos(x/(x*x+y*y)**0.5)
    return theta, phi


def get_volume_of_sphere(radius: float) -> float:
    """ Return volume of a sphere given radius """
    return 4/3 * math.pi * radius**3


def distance_multiple_sets(x1: list[float], y1: list[float], z1: list[float], x2: list[float], y2: list[float], z2: list[float]):
    """ return Euclidean distances """
    p1 = np.array([x1, y1, z1])
    p2 = np.array([x2, y2, z2])
    squared_dist = np.sum((p1-p2)**2, axis=0)
    return np.sqrt(squared_dist)


def distance_between_points(p1: tuple[float], p2: tuple[float]):
    """ return Euclidean distance """
    p1 = np.array(p1)
    p2 = np.array(p2)
    return np.linalg.norm(p1 - p2)
