import bct_py as bct

# Initialize a matrix with the data in example.dat
m = []
f = open("example.dat")
for i in range(30):
    m.append([])
    for j in range(30):
        m[i].append(int(f.readline()))

# Display the matrix
bct.printf(m, "%g")

# Calculate degree distribution
deg, id, od = bct.degrees_dir(m)

# Display the results
bct.printf(id, "%g")
bct.printf(od, "%g")
bct.printf(deg, "%g")

# No memory management necessary
