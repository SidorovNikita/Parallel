import pandas as pd
import os
import numpy as np
import time
import multiprocessing as mp

def task(args):
    pid = args[0]
    num_threads = args[1]
    test = args[2]
    capacity = args[3]

    N = test.shape[0]
    values = test[:, 0]
    prices = test[:, 1]

    s = 2**N // num_threads
    start = s * pid
    count = s if (pid != (num_threads - 1)) else s + 2**N % num_threads

    max_price = 0
    best_set = 0

    for i in range(start, start + count):
        idx = np.array(list(map(int, list(format(i, "0" + str(N) + "b")))))
        price = prices[idx == 1].sum()
        if values[idx == 1].sum() <= capacity and price > max_price:
            max_price = price
            best_set = i

    return best_set, max_price


if __name__ == '__main__':
    files = os.listdir("Tests")
    for filename in files:
        test = pd.read_csv("Tests/" + filename, header=None, delimiter=";").values
        start_time = time.time()
        num_threads = 3

        N = test.shape[0]
        values = test[:, 0]
        capacity = values.sum() / 2

        with mp.Pool(num_threads) as p:
            res = p.map(task, [(pid, num_threads, test, capacity) for pid in range(num_threads)])
            res = sorted(res, key=lambda x: x[1], reverse=True)

        print(filename, format(res[0][0], "0" + str(N) + "b"), time.time() - start_time)
        answer = format(res[0][0], "0" + str(N) + "b").replace('0', '0;').replace('1','1;')
        file = open("Results/bpresult_" + filename, "w")
        file.write(answer[:-1])