import bct_gsl as bct

# Initialize a matrix with the data in example.dat
m_py = []
f = open("example.dat")
for i in range(30):
    m_py.append([])
    for j in range(30):
        m_py[i].append(int(f.readline()))
m = bct.to_gslm(m_py)

# Display the matrix
bct.printf(m, "%g")

# Calculate degree distribution
deg, id, od = bct.degrees_dir(m)

# Display the results
bct.printf(id, "%g")
bct.printf(od, "%g")
bct.printf(deg, "%g")

# Free all memory
bct.gsl_free(m)
bct.gsl_free(id)
bct.gsl_free(od)
bct.gsl_free(deg)
