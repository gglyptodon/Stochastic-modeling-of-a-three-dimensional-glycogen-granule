import json
from pathlib import Path
from unittest import TestCase
from glycogen_module.core import GlycogenStructure, Model
from glycogen_module.utils import distance_between_points, distance_multiple_sets
import pytest


class TestDistMethods(TestCase):
    SCRIPTDIR = Path(__file__).parent.resolve()
    PARAMS_PATH = Path(f"{SCRIPTDIR}/testdata/parameters_1.json")
    PREMADE_G = Path(f"{SCRIPTDIR}/testdata/test_a_b_zero.json")
    PREMADE_G_7_CHAINS = Path(f"{SCRIPTDIR}/testdata/test_7chains.json")
    PREMADE_TAB1 = Path(
        f"{SCRIPTDIR}/testdata/g_N2000_gamma_0.2_1.0_exported_tab1_seed123.json")
    PREMADE_FIG5 = Path(
        f"{SCRIPTDIR}/testdata/g_N1000_gamma_0.2_1.0_exported_fig5_seed123.json")

    def test_dist_single_same(self):
        pos1 = (1, 2, 3)
        pos2 = (1, 2, 3)
        d1_2 = distance_between_points(pos1, pos2)
        print(d1_2)
        self.assertEqual(d1_2, 0)

    def test_dist_single(self):
        pos1 = (2, 3, 4)
        pos2 = (1, 2, 3)
        d1_2 = distance_between_points(pos1, pos2)
        self.assertAlmostEqual(d1_2, 1.7320508)

    def test_dist_single_origin(self):
        pos1 = (2, 3, 4)
        origin = (0, 0, 0)
        d1_2 = distance_between_points(pos1, origin)
        print(d1_2)
        self.assertAlmostEqual(d1_2, 5.3851648)

    def test_dist_multi(self):
        # Some premade glycogen structure
        g = GlycogenStructure.from_json_file(
            TestDistMethods.PREMADE_FIG5, no_init=True)
        print(g)
        x, y, z = g.get_glucose_positions_by_dimension()
        dists_orig = distance_multiple_sets(x, y, z, [0], [0], [0])

        for i in range(len(x)):
            xp, yp, zp = x[i], y[i], z[i]
            self.assertAlmostEqual(
                dists_orig[i], distance_between_points((xp, yp, zp), (0, 0, 0)))

    def test_get_number_of_units_within_radius_n(self, radius=10):
        # Some premade glycogen structure
        g = GlycogenStructure.from_json_file(
            TestDistMethods.PREMADE_FIG5, no_init=True)
        print(g)
        cx, cy, cz = g.get_mean_by_dimension()
        x, y, z = g.get_glucose_positions_by_dimension()
        distances = distance_multiple_sets(x, y, z, [cx], [cy], [cz])
        print(distances)
        within_radius = len([d for d in distances if d < radius])
        self.assertEqual(within_radius, 14)

    def test_get_number_of_units_within_radius(self, radius=10):
        # Some premade glycogen structure
        g = GlycogenStructure.from_json_file(
            TestDistMethods.PREMADE_FIG5, no_init=True)
        print(g)
        cx, cy, cz = g.get_mean_by_dimension()
        x, y, z = g.get_glucose_positions_by_dimension()
        distances = distance_multiple_sets(x, y, z, [cx], [cy], [cz])
        print(distances)
        within_radius = len([d for d in distances if d < radius])
        self.assertEqual(within_radius, 14)


    #def test_get_number_of_units_within_radius_cmp(self, radius=10):
    #    # Some premade glycogen structure
    #    g = GlycogenStructure.from_json_file(
    #        TestDistMethods.PREMADE_FIG5, no_init=True)
    #    print(g)
    #    old = g.get_number_of_units_within_radius()
    #    new = g.get_number_of_units_within_radius_n()
#
#        print(old, new)
 #       self.assertEqual(old, new)