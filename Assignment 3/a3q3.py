import numpy as np 
from scipy import integrate
from matplotlib import pyplot as plt


dat = np.loadtxt("dish_zenith.txt").T # Importing raw data


# Storing it in separate array


# Plotting the initial raw data
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(dat[0], dat[1], dat[2])
ax.set_xlabel('X-axis')
ax.set_ylabel('Y-axis')
ax.set_zlabel('Z-axis')
plt.savefig("figs/q3_rawdata.jpg")
#plt.show()

# Our x, y, z data
xs = dat[0]
ys = dat[1]
zs = dat[2]


cross = np.zeros([len(xs), len(ys)])




# Defining an array with our 4 different functions so we can iterate through them when making matrix A
funcs = np.array(
    [
        lambda x, y : x**2 + y**2,
        lambda x, y: x,
        lambda x, y: y,
        lambda x, y: 1,
    ]
)



A = np.zeros([len(xs), len(funcs)])

for func in funcs:

    func_data = []

    for i in range(len(xs)):
        func_data.append([func(xs[i], j) for j in ys])

    A.append(func_data)

A = np.array(A)



print(np.shape(A.T))

lhs = A.T@A
rhs = A.T@z
fitp = np.linalg.inv(lhs)@rhs
pred=A@fitp

print(fitp)

assert(1==0)
plt.clf()
plt.plot(t,x,'*')
plt.plot(t,pred)
plt.show()
print('squared residual is ',np.sum((x-pred)**2))

