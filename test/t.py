import cmath
import math

def bit_reverse_permutation(x):
    N = len(x)
    bits = int(math.log2(N))

    reversed_indices = [0] * N
    for i in range(N):
        reversed_indices[i] = int('{:0{width}b}'.format(i, width=bits)[::-1], 2)

    # Apply bit-reversal permutation
    result = [0] * N
    for i in range(N):
        result[i] = x[reversed_indices[i]]

    return result

# Example usage:
x = [0, 1, 2, 3, 4, 5, 6, 7]  # Input sequence

# Design the pipe operations
temp_x = x

# Rearrange input sequence using bit-reversal permutation
temp_x = bit_reverse_permutation(temp_x)

# Perform radix-2 butterflies
N = len(temp_x)
levels = int(math.log2(N))

for level in range(1, levels + 1):
    step = 2 ** level
    w = cmath.exp(-2j * cmath.pi / step)

    for k in range(0, N, step):
        for j in range(k, k + step // 2):
            u = temp_x[j]
            v = temp_x[j + step // 2] * cmath.exp(-1j * cmath.pi * (j - k) / (step // 2))

            temp_x[j] = u + v
            temp_x[j + step // 2] = u - v

X = temp_x

print(X)
