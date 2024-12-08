import unittest
from unittest.mock import patch, MagicMock

import galois

from ntru.ntru import NTRU


class NTRUTest(unittest.TestCase):

    @patch('ntru.ntru.NTRU._generate_ternary_poly_coeffs')
    def test_known_correct(self, mock_gen_poly_coeffs: MagicMock):
        mock_gen_poly_coeffs.side_effect = [
            [1, 0, -1, 1, 1, 0, -1], # f(x) = x^6-x^4+x^3+x^2-1
            [1, 0, 1, 0, -1, -1, 0], # g(x) = x^+x^4-x^2-x
            [1, -1, 0, 0, 0, 1, -1]  # r(x) = x^6 -x^6+x-1
        ]

        instance = NTRU(7, 3, 41, 2)
        msg_values = [0, -1, 0, 1, 1, -1, 1]
        msg = galois.Poly(msg_values, field=galois.GF(3))
        cipher = instance.encrypt(msg)
        clear = instance.decrypt(cipher)

        self.assertListEqual(msg.coeffs.tolist(), clear.coeffs.tolist())

    def test_correct_result(self):
        instance = NTRU(23, 11, 313, 4)
        msg_values = [3,6,1,-3,6,3,-2,4,2,0,1,5,-4,1,-1,0,0,0,0,2,0,1,1]
        msg = galois.Poly(msg_values, field=galois.GF(11))
        cipher = instance.encrypt(msg)
        clear = instance.decrypt(cipher)
        self.assertListEqual(msg.coeffs.tolist(), clear.coeffs.tolist())


if __name__ == '__main__':
    unittest.main()
