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
plt.show()

# Our x, y, z data
xs = dat[0]
ys = dat[1]
zs = dat[2]

xy_vec = np.array([xs, ys]).T



# Defining an array with our 4 different functions so we can iterate through them when making matrix A
funcs = np.array(
    [
        lambda x : x[0]**2 + x[1]**2,
        lambda x : x[0],
        lambda x : x[1],
        lambda x : 1,
    ]
)

#print(np.array([funcs[1](j) for j in xy_vec]))

#assert(0==1)


A = np.zeros([len(xs), len(funcs)])


for i in range(len(funcs)):
    A[:, i] = np.array([funcs[i](j) for j in xy_vec])




#print(np.shape(A.T))

lhs = A.T@A
rhs = A.T@zs
fitp = np.linalg.inv(lhs)@rhs
#pred = A@fitp


def paraboloid(x, y, fit):
    return fit[0]*(x**2 + y**2) + fit[1]*x + fit[2]*y + fit[3]




# Plotting the surface of the paraboloid with the original data
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
x = np.linspace(min(xs), max(xs), 100)
y = np.linspace(min(ys), max(ys), 100)
X, Y = np.meshgrid(x, y)
zs = np.array(paraboloid(np.ravel(X), np.ravel(Y), fitp))
Z = zs.reshape(X.shape)

ax.plot_surface(X, Y, Z)

ax.scatter(
    dat[0], 
    dat[1], 
    dat[2],
    c = 'r'
    )

ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
ax.legend()
plt.savefig('figs/a3q3_paraboloid_data.jpg')
plt.show()

