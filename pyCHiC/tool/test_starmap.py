

import time
import multiprocessing
from functools import partial

def f(a, b, c):
    time.sleep(5)
    print("{} {} {}".format(a["z"], b, c))

def main():
    iterable = [1, 2, 3, 4, 5]
    pool = multiprocessing.Pool(3)
    a = {"z":"caca"}
    b = "there"
    func = partial(f, a, b)
    time.sleep(5)
    pool.map(func, iterable)
    pool.close()
    pool.join()

if __name__ == "__main__":
    main()