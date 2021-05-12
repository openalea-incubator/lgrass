"""
Class interface for stand generation.
Adapted from adelwheat. C Fournier
"""

from math import sqrt
from numpy import cos, sin, pi
from scipy.interpolate import interp1d
from operator import itemgetter
from random import random, sample


def randomise_position(position, radius):
    az = random() * 2 * pi
    r = random() * radius
    dx = r * cos(az)
    dy = r * sin(az)
    x, y, z = position
    return x + dx, y + dy, z


def regular(nb_plants, nb_rank, dx, dy, nx=None):
    if nx is None:
        nx = int(nb_plants / nb_rank)
    ny = nb_rank
    domain = ((0, 0), (nx * dx, ny * dy))
    return [(i * dx + dx / 2., j * dy + dy / 2., 0.) for j in range(ny) for i in range(nx)], domain


def regular_plot(inter_plant, inter_row, nrow, plant_per_row, noise=0, convunit=100, center_scene=True):
    dx = inter_plant * convunit
    dy = inter_row * convunit
    positions, domain = regular(nrow * plant_per_row, int(nrow), dx, dy, int(plant_per_row))
    domain_area = abs(domain[1][0] - domain[0][0]) / convunit * abs(domain[1][1] - domain[0][1]) / convunit

    # sorting by ranks
    positions = sorted(positions, key=itemgetter(1, 0))
    # add noise
    if noise > 0:
        positions = map(lambda x: randomise_position(x, noise * convunit), positions)
    if center_scene:
        xc = float(domain[1][0] + domain[0][0]) / 2
        yc = float(domain[1][1] + domain[0][1]) / 2
        positions = [(x - xc, y - yc, z) for x, y, z in positions]
        domain = ((domain[0][0] - xc, domain[0][1] - yc), (domain[1][0] - xc, domain[1][1] - yc))

    return positions, domain, domain_area


class AgronomicStand(object):
    def __init__(self, sowing_density=10, plant_density=10, inter_row=0.8,
                 noise=0, density_curve_data=None):
        self.sowing_density = sowing_density
        self.inter_row = inter_row
        self.plant_density = plant_density
        self.inter_plant = 1. / inter_row / sowing_density
        self.noise = noise
        self.density_curve_data = density_curve_data
        df = density_curve_data
        if df is None:
            self.density_curve = None
        else:
            # hs_curve = interp1d(df['HS'], df['density'])
            TT_curve = interp1d(df['TT'], df['density'])
            # self.density_curve = {'hs_curve':hs_curve,'TT_curve':TT_curve}
            self.density_curve = TT_curve

    def smart_stand(self, nplants=1, at=None, convunit=100):
        """ return an (almost) square stand that match inter-row, current density and nplants in the stand,
             but (dynamicaly) adjusting inter-plant to solve the problem
        """

        density = self.plant_density
        if at is not None:
            if self.density_curve is not None:
                density = self.density_curve(at)

        # find a square design for sowing
        nsown = nplants * 1. * self.sowing_density / density
        side = sqrt(1. / self.sowing_density * nsown)
        nrow = int(max(1, round(side / self.inter_row)))
        plant_per_row = int(max(1, round(float(nsown) / nrow)))
        while nplants > (nrow * plant_per_row):
            plant_per_row += 1
        domain_area = nrow * self.inter_row * plant_per_row * self.inter_plant
        # adjust inter_plant spacing so that n_emerged / domain_area match plant density
        n_emerged = int(round(domain_area * density))
        # assert(n_emerged >= nplants)
        target_domain_area = 1. * n_emerged / density
        inter_plant = target_domain_area / (
                plant_per_row * nrow * self.inter_row)

        positions, domain, domain_area = regular_plot(inter_plant,
                                                      self.inter_row, nrow,
                                                      plant_per_row,
                                                      noise=self.noise,
                                                      convunit=convunit)

        positions = sample(positions, int(n_emerged))
        return n_emerged, domain, positions, domain_area, nrow, plant_per_row
