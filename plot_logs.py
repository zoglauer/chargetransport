import numpy as np
import glob

# Pull phi values from the generated log files. One value of phi per file
# Then, sort s.t. the first entry is the first iteration, i.e. phi_0
phi_paths = glob.glob('/home/jacqueline/chargetransport/CorysCode/planDet/logs' + '/phi_*')
phi_paths = sorted(phi_paths, key=lambda x: float(x.split('_')[1]))


phi = []
for phipath in phi_paths:
    f = open(phipath,'r')
    data = np.genfromtxt(f)
    phi.append(data)


# Pull phi values from the generated log files. One value of phi residual per file
# Then, sort s.t. the first entry is the first iteration, i.e. phiFinalRes_0
phires_paths = glob.glob('/home/jacqueline/chargetransport/CorysCode/planDet/logs' + '/phiFinalRes_*')
phires_paths = sorted(phires_paths, key=lambda x: float(x.split('_')[1]))

phi_res = []
for phirespath in phires_paths:
    f = open(phirespath,'r')
    data = np.genfromtxt(f)
    phi_res.append(data)

import matplotlib.pyplot as plt
fig = plt.figure()
ax = fig.add_subplot(211)  
ax.plot(range(len(phi)),phi)
ax.set_ylabel('phi')
#ax.set_xlabel('iteration')

ax2 = fig.add_subplot(212)
ax2.plot(range(len(phi_res)),phi_res)
ax2.set_ylabel('phi_res')
ax2.set_xlabel('iteration')

plt.subplots_adjust(bottom=0.1)
plt.show()
