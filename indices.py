import random
import math
import pickle
import numpy as np
from scipy.stats import norm, uniform
# random seed generator for RandUCB
from numpy.random import Generator, PCG64


seed = 100

# Gittins index
G = pickle.load(open("./gittins.pkl", "rb"))

def GI(sk0, skt, fk0, fkt, T, t, K, **kw):
    return G[sk0 + skt][fk0 + fkt]

# Current belief
def CB(sk0, skt, fk0, fkt, T, t, K, **kw):
    return (sk0 + skt) / (sk0 + fk0 + skt + fkt)

# Fixed (equal) randomization
def FR(sk0, skt, fk0, fkt, T, t, K, **kw):
    return random.random()

# Raw Thompson Sampling, just sampling from beta distribution
def RTS(sk0, skt, fk0, fkt, T, t, K, **kw):
    return random.betavariate(sk0 + skt, fk0 + fkt)

# Upper condfidence bound index
def UCB(sk0, skt, fk0, fkt, T, t, K, **kw):
    return (sk0 + skt) / (sk0 + fk0 + skt + fkt) + math.sqrt(2 * math.log(max(1,t))/(sk0 + fk0 + skt + fkt))

# Randomized belief index
def RBI(sk0, skt, fk0, fkt, T, t, K, **kw):
    return (sk0 + skt) / (sk0 + fk0 + skt + fkt) + random.expovariate(1/K) * K / (sk0 + fk0 + skt + fkt)

# Randomized Gittins index
def RGI(sk0, skt, fk0, fkt, T, t, K, **kw):
    return G[sk0 + skt][fk0 + fkt] + random.expovariate(1/K) * K / (sk0 + fk0 + skt + fkt)

# Thompson sampling
pre_calc_indices_TS = []
def TS(sk0, skt, fk0, fkt, T, t, K, **kw):
    global pre_calc_indices_TS
    S = kw["S"]
    F = kw["F"]
    n = kw["n"]
    arms_prior = kw["arms_prior"]
    k = kw["k"]

    simulation_time = 200

    if k == 0:
        probs = []
        s = []
        f = []
        for kk in range(K):
            s.append(arms_prior[kk][0] + S[n][t][kk])
            f.append(arms_prior[kk][1] + F[n][t][kk])

        max_time = [0] * K
        for _ in range(simulation_time):
            sample = []
            for kk in range(K):
                sample.append(random.betavariate(s[kk], f[kk]))
            max_time[argmax(sample)] += 1
        pre_calc_indices_TS = sample_from_probs(normalize(list(map(lambda p: math.pow(p, t/(2*T)), normalize(max_time)))))

    return pre_calc_indices_TS[k]


# return a one-hot vector shows the choice of arms
def sample_from_probs(probs):
    assert abs(sum(probs) - 1.0) < 0.000001
    K = len(probs)
    r = random.random()
    indices = [0.0 for _ in range(K)]
    acc = 0.0
    for kk in range(K):
        acc += probs[kk]
        if acc >= r:
            indices[kk] = 1.0
            break

    return indices
    
def argmax(a):
    return max(list(range(len(a))), key=lambda x: a[x])
    
    
def argmax_list(a):
    max_a = max(a)
    l = len(a)
    return list(filter(lambda x: a[x] == max_a, list(range(l))))
    
def normalize(probs):
    s = sum(probs)
    for p in probs:
        assert p >= 0
    if s == 0:
        s = 1;
    return list(map(lambda x: x/s, probs))


# Randomized Upper condfidence bound index
def choice(options,probs):
    x = np.random.rand()
    cum = 0
    for i,p in enumerate(probs):
        cum += p
        if x < cum:
            break
    return options[i]

def RandUCB(sk0, skt, fk0, fkt, T, t, K, probs, **kw):
    mu_hat = (sk0 + skt) / (sk0 + fk0 + skt + fkt)
    # m = np.random.choice(20, p=probs)
    m = choice(range(20), probs)
    zt = m/19 
    c = math.sqrt(1/(sk0 + fk0 + skt + fkt)) 
    return  mu_hat + zt*c


