import math
import random

import galois
import numpy as np

from utils import isprime


class NTRU:

    def __init__(self, n: int, p: int, q: int, d: int):
        if n <= 0 or p <= 0 or q <= 0 or d <= 0:
            raise Exception("All values must be positive")

        gcd_nq = math.gcd(n, q)
        gcd_pq = math.gcd(p, q)
        if gcd_nq != gcd_pq != 1 or q <= (6 * d + 1) * p:
            raise Exception("q must hold gcd(n,q)=gcd(p,q)=1 and q > (6d+1)p")

        if not isprime(p):
            raise Exception("p must be prime")

        self.n = n
        self.p = p
        self.q = q
        self.d = d
        self.field_p = galois.GF(p)
        self.field_q = galois.GF(q)
        self.base_poly_values = [1] + [0] * (self.n - 1) + [-1]
        self.base_poly_in_p = galois.Poly(self.base_poly_values, field=self.field_p)
        self.base_poly_in_q = galois.Poly(self.base_poly_values, field=self.field_q)

        self.f_values = self.__generate_ternary_poly(d + 1, d, True)
        self.g_value = self.__generate_ternary_poly(d, d)

        self.f_p_inv_poly = galois.egcd(galois.Poly(self.f_values, field=self.field_p),
                                        galois.Poly(self.base_poly_values, field=self.field_p))[1] % self.base_poly_in_p
        self.f_q_inv_poly = galois.egcd(galois.Poly(self.f_values, field=self.field_q),
                                        galois.Poly(self.base_poly_values, field=self.field_q))[1] % self.base_poly_in_q

        self.h = self.f_q_inv_poly * galois.Poly(self.g_value, field=self.field_q)

    def __generate_ternary_poly(self, d1: int, d2: int, must_have_inv: bool = False):
        """
        Generates a polynomial which holds:
        - it has `d1` number of coefficients = 1
        - it has `d2` number of coefficients = -1
        - has an inverse in R_p and R_q according to `has_inverse`
        - its degree = n
        :param d1: number of coefficients = 1
        :param d2: number of coefficients = 2
        :param must_have_inv: does the polynomial have to have an inverse in R_p and R_q
        :return: a polynomial that matches this requirement
        """
        values = [-2] * self.n
        count = 0
        while count <= d1:
            pos = random.randint(0, self.n - 1)
            if values[pos] == -2:
                values[pos] = 1
                count += 1
        count = 0
        while count <= d2:
            pos = random.randint(0, self.n - 1)
            if values[pos] == -2:
                values[pos] = -1
                count += 1

        for x in range(self.n):
            if values[x] == -2:
                values[x] = 0

        if not must_have_inv:
            return values

        poly_p = galois.Poly(values, field=self.field_p)
        base_poly_p = galois.Poly(self.base_poly_values, field=self.field_p)
        if galois.gcd(poly_p, base_poly_p) != 1:  # Need to generate new polynomials
            return self.__generate_ternary_poly(d1, d2, must_have_inv)

        poly_q = galois.Poly(values, field=self.field_q)
        base_poly_q = galois.Poly(self.base_poly_values, field=self.field_q)
        if galois.gcd(poly_q, base_poly_q) != 1:  # Need to generate new polynomials
            return self.__generate_ternary_poly(d1, d2, must_have_inv)

        return values