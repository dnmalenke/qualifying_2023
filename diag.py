import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib


def butterfly(N):
    num_stages = int(math.log2(N))
    for stage in range(num_stages):
        num_blocks = int(N/2**(stage+1))
        for block in range(num_blocks):
            for n in range(2**stage):
                index1 = n + 2**stage*block*2
                index2 = n + 2**stage*block*2 + 2**stage
                omega = np.exp(-2*np.pi*1j*n/(2**(stage+1)))
                plt.arrow(stage, index1, 1, 0, head_width=0.2, head_length=0.1, fc='k', ec='k')
                plt.arrow(stage, index2, 1, 0, head_width=0.2, head_length=0.1, fc='k', ec='k')
                plt.arrow(stage, index1, 1, 2**stage, head_width=0.2, head_length=0.1, fc='k', ec='k', linestyle='--')
                plt.arrow(stage, index2, 1, -2**stage, head_width=0.2, head_length=0.1, fc='k', ec='k', linestyle='--')
                # plt.arrow(stage, index1, (2**stage)*omega.imag, (2**stage)*omega.real, head_width=0.2, head_length=0.1, fc='r', ec='r')
                # plt.arrow(stage, index2, (2**stage)*omega.imag, (2**stage)*omega.real, head_width=0.2, head_length=0.1, fc='b', ec='b')

    plt.xlim([0, num_stages])
    plt.ylim([0, N-1])
    plt.gca().invert_yaxis()
    plt.gca().invert_xaxis()
    plt.show()

if __name__ == "__main__":
    matplotlib.use('TkAgg')
    butterfly(64)