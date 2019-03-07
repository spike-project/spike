#!/usr/bin/env python 
# encoding: utf-8

"""
counter.py

This utility declare two useful decorators

counting and timeit

usage :
@counting
def func()
    ...

counting : 
    The number of time the decorated function is called is counted,  and the count is printed at the end of the prgm execution
timeit :
    As counting, but additionally, the total (elapsed) time spent in the function is totalized and printed at the end of the prgm


Created by Marc-André on 2012-11-05.

Copyright (c) 2012 IGBMC. All rights reserved.
"""

################################################################################
# UTILITIES

from __future__ import print_function

class Counter(object):
    """
    defines an incremental counter
    c = Counter()     # actually defines a function c()
    v = c()           # each time c() is called, the internal value is incremented
    value is found at c.count
    """
    def __init__(self, ini=0):
        self.count = ini
    def __call__(self):
        self.count += 1
        return self.count
def counting(func_to_decorate):
    """define a decorator that count the nmber of function/method call, and display at the end of the prgm"""
    import functools
    import atexit
    count = Counter()
    @atexit.register
    def report():
        print("function %s called %d times"%(func_to_decorate.__name__, count.count))
    @functools.wraps(func_to_decorate)
    def wrapper(*args, **kw):
        result = func_to_decorate(*args, **kw)
        count()
        return result
    return wrapper
def timeit(func_to_decorate):
    """
    define a decorator that count the total time passed in the function/method call, and display at the end of the prgm
    NOT PRECISE for very short times
    """
    import functools
    import atexit
    import time
    count = Counter()
    count.tim = 0.0
    @atexit.register
    def report():
        if count.count>0:
            print("function %s  time : %f sec    called %d times    time/call : %f sec"%(func_to_decorate.__name__, count.tim, count.count, count.tim/count.count))
        else:
            print("function %s called %d times"%(func_to_decorate.__name__, count.count))
    @functools.wraps(func_to_decorate)
    def wrapper(*args, **kw):
        t0 = time.time()
        result = func_to_decorate(*args, **kw)
        count.tim += time.time() - t0
        count()
        return result
    return wrapper

#### TESTS and examples #############
def main():
    """
    should print :
    function test2  time : 3.000320 sec    called 3   time/call : 1.000107 sec
    function test called 13 times
    """
    # define functions locally in main() so as not to clutter the callers name space
    @counting
    def test(x):
        "example for counting"
        return x

    @timeit
    def test2(t):
        "example for timeit"
        import time
        time.sleep(t)
        return test(t)

    for i in range(10):
        test(i)
    for i in range(3):
        test2(i)

if __name__ == '__main__':
    main()

