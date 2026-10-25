"""Exact disease expectations per order; only order integration is sampled.

The full permanent-SI kernel and restricted source-target kernel are batched
without importing the existing implementations, allowing independent checks.
One initially clean worker; uniform initial source including unvisited beds;
independent perfect before/after hygiene; no environment or turnover.
"""
from collections import Counter
import math
from numbers import Integral, Real
import time
import numpy as np


def encode_sequence(sequence):
    if isinstance(sequence,(str,bytes)):
        raise ValueError('use an iterable of string labels')
    x=tuple(sequence)
    if not x or any(not isinstance(b,str) or not b.strip() for b in x):
        raise ValueError('nonempty labels required')
    labels=tuple(dict.fromkeys(x));index={b:i for i,b in enumerate(labels)}
    return np.array([index[b] for b in x],dtype=np.int16),labels


def unique_multiset_orders(sequence,limit=10000):
    encoded,labels=encode_sequence(sequence)
    counts=Counter(map(int,encoded));n=len(encoded)
    number=math.factorial(n)//math.prod(math.factorial(c) for c in counts.values())
    if number>limit:
        raise ValueError('unique-order enumeration exceeds declared bound')
    orders=[];path=[]
    def visit():
        if len(path)==n:
            orders.append(path.copy());return
        for label in counts:
            if counts[label]:
                counts[label]-=1;path.append(label);visit();path.pop();counts[label]+=1
    visit()
    return np.array(orders,dtype=np.int16),labels


def sample_multiset_orders(sequence,number,seed):
    if not isinstance(number,Integral) or number<1 or not isinstance(seed,Integral) or seed<0:
        raise ValueError('positive count and nonnegative integer seed required')
    encoded,labels=encode_sequence(sequence)
    rng=np.random.default_rng(int(seed))
    # Each independently shuffled list projects uniformly to multiset orders:
    # each distinct order has the same product of label-count factorials.
    orders=rng.permuted(np.broadcast_to(encoded,(number,len(encoded))),axis=1)
    return orders,labels


def expected_order_outcomes(orders,h,alpha,beta,*,source_pool_size=15,secondary=True,chunk_size=256,deadline=None):
    x=np.asarray(orders)
    if x.ndim!=2 or not x.shape[0] or not x.shape[1] or not np.issubdtype(x.dtype,np.integer) or np.any(x<0):
        raise ValueError('orders must be a nonempty integer matrix')
    m=int(x.max())+1
    if m>8 or set(np.unique(x))!=set(range(m)):
        raise ValueError('use contiguous labels and at most8 visited beds')
    reference=np.bincount(x[0],minlength=m)
    if any(not np.array_equal(np.bincount(row,minlength=m),reference) for row in x):
        raise ValueError('all orders must have the same multiset')
    for value in (h,alpha,beta):
        if not isinstance(value,Real) or not math.isfinite(value) or not 0<=value<=1:
            raise ValueError('finite probabilities required')
    if not isinstance(source_pool_size,Integral) or isinstance(source_pool_size,bool) or source_pool_size<m:
        raise ValueError('source universe must include visited beds')
    if not isinstance(chunk_size,Integral) or chunk_size<1 or not isinstance(secondary,bool):
        raise ValueError('positive chunk size and boolean mechanism required')
    result=np.zeros(len(x));r=1-h
    for start in range(0,len(x),chunk_size):
        if deadline is not None and time.perf_counter()>=deadline:
            raise TimeoutError('phase2 wall budget exhausted')
        paths=x[start:start+chunk_size];count=len(paths)
        total=np.zeros(count)
        if secondary:
            masks=np.arange(1<<m,dtype=np.int32)[None,:]
            clean=np.zeros((count,1<<m));dirty=np.zeros_like(clean)
            clean[:,1<<np.arange(m)]=1/source_pool_size
            clean[:,0]=(source_pool_size-m)/source_pool_size
            for label in paths.T:
                if deadline is not None and time.perf_counter()>=deadline:
                    raise TimeoutError('phase2 wall budget exhausted')
                clean+=h*dirty;dirty*=r
                bit=(1<<label)[:,None];colonized=(masks&bit)!=0
                acquisition=alpha*clean*colonized
                transfer=beta*dirty*~colonized
                clean-=acquisition;dirty+=acquisition-transfer
                dirty+=np.take_along_axis(transfer,masks^bit,axis=1)
                total+=transfer.sum(axis=1)
                clean+=h*dirty;dirty*=r
            if not np.allclose((clean+dirty).sum(axis=1),1,rtol=1e-11,atol=1e-12):
                raise ArithmeticError('full-SI mass not conserved')
        else:
            pairs=[(s,t) for s in range(m) for t in range(m) if s!=t]
            if not pairs:
                continue
            sources=np.array([s for s,t in pairs])[None,:]
            targets=np.array([t for s,t in pairs])[None,:]
            clean=np.ones((count,len(pairs)));dirty=np.zeros_like(clean)
            for label in paths.T:
                if deadline is not None and time.perf_counter()>=deadline:
                    raise TimeoutError('phase2 wall budget exhausted')
                clean+=h*dirty;dirty*=r
                acquisition=alpha*clean*(label[:,None]==sources)
                transfer=beta*dirty*(label[:,None]==targets)
                clean-=acquisition;dirty+=acquisition-transfer
                total+=transfer.sum(axis=1)/source_pool_size
                clean+=h*dirty;dirty*=r
        result[start:start+count]=total
    return result
