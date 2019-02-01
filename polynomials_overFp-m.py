#
# polygf.py : GF(p) で多項式を計算する (p は素数)
#
#             Copyright (C) 2018 Makoto Hiroi
#
import numpy as np
import functools
import itertools
import math

"""
拡張ユークリッドの互除法
引数を整数とし、a * x + b * y = gcd(a, b)
を満たす(x, y)の組を求める。戻り値は
gcd(a, b), x, y
"""
def ext_euclid(a, b):
    xs = a, 1, 0
    ys = b, 0, 1
    while ys[0] != 0:
        q, z = divmod(xs[0], ys[0])
        xs, ys = ys, (z, xs[1] - q * ys[1], xs[2] - q * ys[2])
    return xs
print(ext_euclid(24, 6))
print(ext_euclid(28, 63))
print(ext_euclid(24, 27))
print(ext_euclid(-24, 6))

# 逆数表の生成
def make_inv(n):
    inv = [0] * n
    for x in range(1, n):
        g, y, _ = ext_euclid(x, n)
        if g != 1:
            raise Exception("not prime number!")
        inv[x] = (y + n) % n
    return inv
print(make_inv(7))
print(make_inv(3))
# print(make_inv(24))
# "not prime number!"

"""
GF(p)で多項式を計算する（p：素数）
"""
class GF:
    def __init__(self, n):
        self.num = n
        self.inv = make_inv(n)

    # 正規化
    def normalize(self, xs):
        xs %= self.num
        xs = np.trim_zeros(xs, 'f')
        if not len(xs):
            xs = np.array([0])
        return xs

    # 足し算
    def polyadd(self, xs, ys):
        xn = len(xs)
        yn = len(ys)
        if xn >= yn:
            zs = xs.copy()
            zs[-yn:] += ys
        else:
            zs = ys.copy()
            zs[-xn:] += xs
        return self.normalize(zs)
    
    # 引き算
    def polysub(self, xs, ys):
        return self.polyadd(xs, self.num - ys)
    
    # 掛け算
    def polymul(self, xs, ys):
        xn = len(xs)
        yn = len(ys)
        zn = xn + yn - 1
        zs = np.full(zn, 0, dtype=np.int32)
        for i in range(yn):
            zs[i:i+xn] += ys[i] * xs
        return self.normalize(zs)
    
    # 割り算
    def polydiv(self, xs, ys):
        xn = len(xs)
        yn = len(ys)
        zs = xs.copy()
        qs = []
        for _ in range(xn - yn + 1):
            temp = (zs[0] * self.inv[ys[0]]) % self.num
            zs[:yn] += self.num - temp * ys
            qs.append(temp)
            zs = zs[1:] % self.num
        if qs == []: qs = [0]
        return np.array(qs), self.normalize(zs)


gf3 = GF(3)

fx = np.array([1, 2, 1]) # x^2 + 2 x +1
gx = np.array([1, 0, 2]) # x^2 + 2
print(gf3.polyadd(fx, gx))
print(gf3.polysub(fx, gx))
print(gf3.polymul(fx, gx))
p, q = gf3.polydiv(fx, gx)#商と余り
print(p, q)

fx = np.array([1, 0,  2]) # x^2 - 1 = (x + 1)(x - 1) 
gx = np.array([1, 1]) # x + 1
print(gf3.polyadd(fx, gx))
print(gf3.polysub(fx, gx))
print(gf3.polymul(fx, gx))
p, q = gf3.polydiv(fx, gx)#商と余り
print(p, q)