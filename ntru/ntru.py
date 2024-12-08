import math
import random

from galois import Poly, GF, egcd, gcd

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
        self.field_p = GF(p)
        self.field_q = GF(q)
        self.base_poly_values = [1] + [0] * (self.n - 1) + [-1]
        self.base_poly_in_p = Poly(self.base_poly_values, field=self.field_p)
        self.base_poly_in_q = Poly(self.base_poly_values, field=self.field_q)

        self.f_values = self._generate_ternary_poly_coeffs(d + 1, d, True)
        self.g_value = self._generate_ternary_poly_coeffs(d, d)

        [self.center_lift_mapping_p, self.center_lift_mapping_q] = self.__generate_center_lift_mapping()

        self.f_p_inv_poly = egcd(Poly(self.f_values, field=self.field_p),
                                 Poly(self.base_poly_values, field=self.field_p))[1] % self.base_poly_in_p
        self.f_q_inv_poly = egcd(Poly(self.f_values, field=self.field_q),
                                 Poly(self.base_poly_values, field=self.field_q))[1] % self.base_poly_in_q

        self.h = self.f_q_inv_poly * Poly(self.g_value, field=self.field_q) % self.base_poly_in_q

    def _generate_ternary_poly_coeffs(self, d1: int, d2: int, must_have_inv: bool = False):
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

        poly_p = Poly(values, field=self.field_p)
        base_poly_p = Poly(self.base_poly_values, field=self.field_p)
        if gcd(poly_p, base_poly_p) != 1:  # Need to generate new polynomials
            return self._generate_ternary_poly_coeffs(d1, d2, must_have_inv)

        poly_q = Poly(values, field=self.field_q)
        base_poly_q = Poly(self.base_poly_values, field=self.field_q)
        if gcd(poly_q, base_poly_q) != 1:  # Need to generate new polynomials
            return self._generate_ternary_poly_coeffs(d1, d2, must_have_inv)

        return values

    def __generate_center_lift_mapping(self) -> list[dict[int, int]]:
        return [{(x + self.p) % self.p: x for x in range(-1 * math.ceil(self.p / 2), math.ceil(self.p / 2))},
                {(x + self.q) % self.q: x for x in range(-1 * math.ceil(self.q / 2), math.ceil(self.q / 2))}]

    def _center_lift_poly_coeffs(self, coeffs: list[int], fromP: bool = False) -> list[int]:
        try:
            next(iter([x for x in coeffs if x < 0 or x >= self.p]))
            raise Exception("To center raise we move from Rp to R, poly is not in Rp")
        except Exception as e:
            if e.args == "To center raise we move from Rp to R, poly is not in Rp":
                raise e
            pass
        return [self.center_lift_mapping_p[x] if fromP else self.center_lift_mapping_q[x] for x in coeffs]

    def encrypt(self, message: Poly) -> Poly:
        """
        Encrypts as follows:
        1. Center lift m to R from Rp
        2. Generate a random ephemeral value from Ternary poly class of d,d
        3. compute r * h = x
        4. multiply x by p (p*x = y)
        5. add y to m (m+y=e)
        6. return e mod q
        :param message: the message to encrypt
        :return: a polynomial in Rq that is the encrypted message
        """
        if message.field != self.field_p:
            raise Exception("Message must be from Rp")
        cl_m_poly_in_q = Poly(self._center_lift_poly_coeffs(message.coeffs.tolist(), True), field=self.field_q)
        r_poly_in_q = Poly(self._generate_ternary_poly_coeffs(self.d, self.d), field=self.field_q)
        p_scalar = Poly([self.p], field=self.field_q)  # same same as multiplying by scalar

        return p_scalar * r_poly_in_q * self.h % self.base_poly_in_q + cl_m_poly_in_q

    def decrypt(self, cipher: Poly) -> Poly:
        """
        Decrypts as follows:
        1. Multiplies f by the ciphertext in Rq, a = f*e
        2. Centerlift a to Rp (First centerlift from Rq to R and then mod p) a'
        3. Multiply a' by the inverse of F in p -> m'
        4. m' is the edecrypted data
        :param cipher:
        :return:
        """
        if cipher.field != self.field_q:
            raise Exception("Ciphertext fields is not Rq, can't proceed")

        f_poly_in_q = Poly(self.f_values, field=self.field_q)
        a = f_poly_in_q * cipher % self.base_poly_in_q

        a_poly_in_p_after_cl = Poly([x % self.p for x in self._center_lift_poly_coeffs(a.coeffs.tolist(), False)],
                                    field=self.field_p)

        return self.f_p_inv_poly * a_poly_in_p_after_cl % self.base_poly_in_p
