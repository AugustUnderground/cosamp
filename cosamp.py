import numpy as np
import scipy
import matplotlib.pyplot as plt

# def cosamp(Phi, u, s, tol=1e-10, max_iter=1000):
def cosamp(Phi, u, s, tol=1e-10, max_iter=10):
    max_iter -= 1 # Correct the while loop
    num_precision = 1e-12
    a = np.zeros(Phi.shape[1])
    v = u
    iter = 0
    halt = False
    while not halt:
        iter += 1
        y = abs(np.dot(np.transpose(Phi), v))
        Omega = [i for (i, val) in enumerate(y) if val > np.sort(y)[::-1][2*s] and val > num_precision] # equivalent to below
        #Omega = np.argwhere(y >= np.sort(y)[::-1][2*s] and y > num_precision)
        T = np.union1d(Omega, a.nonzero()[0])
        #T = np.union1d(Omega, T)
        b = np.dot( np.linalg.pinv(Phi[:,T]), u )
        iGood = (abs(b) > np.sort(abs(b))[::-1][s]) & (abs(b) > num_precision)
        T = T[iGood]
        a[T] = b[iGood]
        v = u - np.dot(Phi[:,T], b[iGood])
        halt = np.linalg.norm(v)/np.linalg.norm(u) < tol or \
               iter > max_iter
        print(f"Iteration {iter}: a = {a.nonzero()[0]}\n")
    return a

n = 4096 # number of measurements
t = np.linspace(0.0, 1.0, num=n)

x = np.sin(91*2*np.pi*t) + np.sin(412*2*np.pi*t) # original signal (to be reconstructed)

# randomly sample signal
p = 128 # random sampling (Note that this is one eighth of the Shannon–Nyquist rate!)
aquis = np.round((n-1) * np.random.rand(p)).astype(int)
y = x[aquis] # our compressed measurement from the random sampling

# Here {y} = [C]{x} = [C][Phi]{s}, where Phi is the inverse discrete cosine transform
Phi = scipy.fft.dct(np.eye(n), axis=0, norm='ortho')
CPhi = Phi[aquis,:]

# l1 minimization (through linear programming)
s = cosamp(CPhi, y, 10) # obtain the sparse vector through CoSaMP algorithm
xrec = scipy.fft.idct(s, axis=0, norm='ortho') # Reconstructed signal

s[s.nonzero()[0]]


ss = np.array([6.936376032725912, 27.474683924431055, -27.578700285137383, -11.337040416677038, -7.6673838856754575, -4.563676616381389, -4.8639213347650365, -6.087959111673908, 6.127047371932055, 10.446992690237028, 33.76464111472499, -14.967453149814236, -16.20380157256696, -4.713715149349071, 5.0523166106586945 ])



scipy.fft.idct(ss, axis=0, norm='ortho')

figw, figh = 7.0, 5.0 # figure width and height
plt.figure(figsize=(figw, figh))
plt.plot(t, s)
plt.title('Sparse vector $s$')
plt.show()

# Visualize the compressed-sensing reconstruction signal
figw, figh = 7.0, 5.0 # figure width and height
plt.figure(figsize=(figw, figh))
plt.plot(t, x,    'b', label='Original signal')
plt.plot(t, xrec, 'r', label='Reconstructed signal')
plt.xlim(0.4, 0.5)
legend = plt.legend(loc='upper center', shadow=True, fontsize='x-large')
# Put a nicer background color on the legend.
legend.get_frame().set_facecolor('C0')
plt.show()
