import math
import random

import numpy as np
from scipy import stats
from scipy.stats import rv_discrete

from utils import isprime


# encrypt_debug_info = ""


class LWE:

    def __init__(self, n, p):
        if not (isprime(p) and (n ** 2) < p < (2 * (n ** 2)) and p >= 2):
            raise Exception("p is not a prime ")

        self.n = n
        self.p = p
        self.m = int(1.1 * n * np.log(self.p))
        self.distribution = self.__generate_chi_distribution()

        self.secret = np.random.randint(0, self.p, size=self.n)
        self.equation_vectors = [
            np.random.randint(0, self.p, size=self.n) for _ in range(self.m)
        ]
        self.noise = self.distribution.rvs(size=self.m)

        self.equation_results = [
            (self.equation_vectors[i].dot(self.secret) + self.noise[i]) % self.p for i in range(self.m)
        ]

    def __generate_chi_distribution(self):
        standard_deviation = self.p / (2 * (self.n ** 2))
        norm = stats.norm(scale=standard_deviation, loc=0)
        distribution_range = range(-math.ceil(self.p / 2), self.p // 2)
        sample_points = []

        discrete_values = ([x for x in distribution_range], [0.] * self.p)

        for i in distribution_range:
            if i == -math.ceil(self.p / 2):
                sample_points.append(i - 0.5)
            sample_points.append(i + 0.5)

        cumulative_distribution_func_values = norm.cdf(sample_points)
        for i in range(1, self.p + 1):
            discrete_values[1][i - 1] = (cumulative_distribution_func_values[i] -
                                         cumulative_distribution_func_values[i - 1])

        return rv_discrete(a=-math.ceil(self.p / 2), b=self.p // 2, values=discrete_values)

    def encrypt(self, bit):
        # global encrypt_debug_info
        if bit != 0 and bit != 1:
            raise Exception("Bit must be 0 or 1")

        subset_size = random.randint(1, self.m)
        subset = list(set(random.sample(range(self.m), subset_size)))

        sum_equation = [0] * self.n
        sum_no_mod = [0] * self.n
        for i in subset:
            sum_equation = (sum_equation + self.equation_vectors[i]) % self.p
            sum_no_mod = sum_no_mod + self.equation_vectors[i]
        sum_result = ((0 if bit == 0 else self.p // 2) + sum([self.equation_results[i] for i in subset])) % self.p

        # Debug info for testing if needed
        #         encrypt_debug_info = f"""Bit: {bit}, subset: {subset}, size: {subset_size}
        # Secret: {self.secret}
        # Noise: {self.noise}
        # Equation vectors: {[self.equation_vectors[i] for i in subset]}
        # Equation vectors sum: {sum_no_mod}
        # Equation vectors sum mod p: {sum_equation}
        # Results: {[self.equation_results[i] for i in subset]}
        # Results sum: {(0 if bit == 0 else self.p // 2) + sum([self.equation_results[i] for i in subset])}
        # Results sum mod p: {((0 if bit == 0 else self.p // 2) + sum([self.equation_results[i] for i in subset])) % self.p}
        # """

        return sum_equation, sum_result

    def decrypt(self, sum_equation, sum_result):
        decoded = (sum_result - (sum_equation.dot(self.secret) % self.p)) % self.p
        mid = self.p // 2
        threshold = int(np.round(self.p / 4))
        dist_from_mid = decoded if decoded <= mid else self.p - decoded

        return 0 if dist_from_mid < threshold else 1
