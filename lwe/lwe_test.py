import random

import numpy as np

from lwe.lwe import LWE

pairs = [(10, 127), (14, 337), (19, 613), (40, 3001), (21, 751), (51, 4889), (69, 7829)]

for _ in range(10_000):
    n, p = random.choice(pairs)
    bit = np.random.randint(0, 2, size=1)[0]
    instance = LWE(n, p)
    eqs, res = instance.encrypt(bit)
    decrypted = instance.decrypt(eqs, res)
    assert decrypted == bit
