"""
参考：「お気楽NumPyプログラミング超入門」
http://www.geocities.jp/m_hiroi/light/numpy06.html
"""

import numpy as np

"""
リスト：多項式の四則演算
"""
# 正規化
def normalize(xs):
    xs = np.trim_zeros(xs, 'f')
    #リストが空だった場合、0を返すようにする
    if not len(xs):
        xs = np.array([0.])
    return xs

# 足し算
def polyadd(xs, ys):
    xn = len(xs)
    yn = len(ys)
    # zs[-yn:]の意味は「python リスト コロン」
    # でググるとすぐ出てくる
    if xn >= yn:
        zs = xs.copy()
        zs[-yn:] += ys
    else:
        zs = ys.copy()
        zs[-xn:] += xs
    return normalize(zs)

# 引き算
def polysub(xs, ys):
    return polyadd(xs, -ys)

# 掛け算
def polymul(xs, ys):
    xn = len(xs)
    yn = len(ys)
    zn = xn + yn - 1
    zs = np.full(zn, 0.)
    for i in range(yn):
        zs[i:i + xn] += ys[i] * xs
    return normalize(zs)

# 割り算
def polydiv(xs, ys):
    xn = len(xs)
    yn = len(ys)
    zs = xs.copy()
    qs = []
    for _ in range(xn - yn + 1):
        temp = zs[0] / ys[0]
        zs[:yn] -= temp * ys
        qs.append(temp)
        zs = zs[1:]
    if qs == []: qs = [0.]
    return np.array(qs), normalize(zs)

"""
実行例
"""
fx = np.array([1., 1.]) # x + 1
gx = np.array([1., 1.]) # x + 1
print(polyadd(fx, gx))
print(polysub(fx, gx))
print(polymul(fx, gx))
print(polydiv(fx, gx))

fx = np.array([2., 1.]) # x^2 + 1
gx = np.array([2., 1.]) # x^2 + 1
print(polyadd(fx, gx))
print(polysub(fx, gx))
print(polymul(fx, gx))
print(polydiv(fx, gx))

fx = np.array([1., 5., 6.]) # x^2 + 5x + 6 = (x + 2)(x + 3)
gx = np.array([1., 2.]) # x + 2
print(polyadd(fx, gx))
print(polysub(fx, gx))
print(polymul(fx, gx))
print(polydiv(fx, gx))

