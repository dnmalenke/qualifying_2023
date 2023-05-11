with open("1.txt", "r") as f1, open("2.txt", "r") as f2:
    nums = set(line.split()[1] for line in f1)
    result = [line.split()[0] for line in f2 if line.split()[1] not in nums]

print(result)
