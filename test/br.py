# Generate a list of 1024 10-bit reversed numbers
output = []
for i in range(1024):
    # Reverse the bits of i
    reversed_num = int(bin(i)[2:].zfill(10)[::-1], 2)
    output.append(reversed_num)

# Print the generated list
for i, num in enumerate(output):
    print(f"output[{i}] = {num}")


