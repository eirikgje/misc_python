from __future__ import division
import numpy as np
import numpy.random as rn

def roll(num, typ, add = 0):
    tot = 0
    for i in range(num):
        curr = np.ceil(rn.random() * typ)
        print curr
        tot += curr
    tot = tot + add
    return tot

def avghd(num, typ):
    return np.floor(num * sum(np.arange(1, typ + 1)) / typ)

def fudgeroll(add=0):
    tot = 0
    for i in range(4):
        res = roll(1, 3)
        if res == 3:
            tot += 1
        elif res == 1:
            tot -= 1
    return tot + add
