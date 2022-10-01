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


# Creating A matrix
A = np.zeros([len(xs), len(funcs)])

for i in range(len(funcs)):
    A[:, i] = np.array([funcs[i](j) for j in xy_vec])



# Linear algebra operations
lhs = A.T@A
rhs = A.T@zs
fitp = np.linalg.inv(lhs)@rhs


# Defining paraboloid function
def paraboloid(x, y, fit):
    return fit[0]*(x**2 + y**2) + fit[1]*x + fit[2]*y + fit[3]




# Plotting the surface of the paraboloid with the original data

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
x = np.linspace(min(xs), max(xs), 100)
y = np.linspace(min(ys), max(ys), 100)
X, Y = np.meshgrid(x, y)
zshape = np.array(paraboloid(np.ravel(X), np.ravel(Y), fitp))
Z = zshape.reshape(X.shape)

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
plt.savefig('figs/a3q3_paraboloid_data.jpg')
plt.show()

print("Our best fit parameters are:", fitp)


# Coordinate transformation function
def coord_trans(x):
    A, B, C, D = x
    return np.array([
        A,
        -B/(2*A),
        -C/(2*A),
        -(B**2 + C**2 - 4*A*D)/(4*A)
    ])

print("Our best fit parameters (in the original coordinate system) are:", coord_trans(fitp))



# Now we find the error in the data

original_fit = coord_trans(fitp)

noise  = np.std(np.array([paraboloid(j[0], j[1], fitp) for j in xy_vec]) - zs)

trial_err = 0.0000024
fitp[0] = fitp[0]+trial_err

noise2  = np.std(np.array([paraboloid(j[0], j[1], fitp) for j in xy_vec]) - zs)

print("The noise is:", noise)
print("New noise is:", noise2)

print("The estimated error on a is", trial_err)

print("focal length is:", 1/original_fit[0]/4, '\pm', 4/original_fit[0]**2*trial_err)