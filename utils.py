import random


def mod_pow(b: int, n: int, m: int):
    return pow(b, n, mod=m)


def isprime(p: int, rounds: int = 10) -> bool:
    """
    Miller Rabin test for primality - 10 rounds
    Source: https://en.wikipedia.org/wiki/Miller–Rabin_primality_test#Miller–Rabin_test
    :param p: the value to test
    :return: true if it is, false otherwise
    """
    # Lower side edge cases
    if (p % 2 == 0 and p != 2) or p == 1: return False
    if p == 2 or p == 3: return True

    p_copy = p - 1
    pow_counter = 0  # s
    while p_copy % 2 == 0:
        p_copy /= 2
        pow_counter += 1
    mul_remainder = p_copy  # d

    for i in range(rounds):
        a = random.randint(2, p - 2)
        x = mod_pow(a, mul_remainder, p)
        for j in range(pow_counter):
            y = mod_pow(x,2,p)
            if y == 1 and x != 1 and x != p-1:
                return False
            x = y
        if y != 1:
            return False

    return True