import numpy as np

symbols=[]
try:
    coor=np.genfromtxt('coor.xyz', comments="H")[:,1:]
    n_atoms=coor.shape[0]
    for i in range(n_atoms):
        symbols.append(np.genfromtxt('coor.xyz', comments="H",dtype="|U5")[i][0])    
except:
    coor=np.vstack([np.genfromtxt('coor.xyz', comments="H")[1:],np.zeros(3)])
    n_atoms=1
    symbols.append(np.genfromtxt('coor.xyz', comments="H",dtype="|U5")[0])


x=np.loadtxt("x.esp")
y=np.loadtxt("y.esp")
z=np.loadtxt("z.esp")
grid=np.loadtxt("grid.dat")
no=np.loadtxt("0.esp")
n_gridpoints=grid.shape[0]
r=np.zeros((n_atoms,n_gridpoints,4))
units=1.6022*100/(4*np.pi*8.854188*2.7211385)
X=np.zeros((n_atoms,n_gridpoints,3))
field=0.0008

for k in range(n_atoms):
    for i in range(n_gridpoints):
        for j in range(3):
            r[k,i,j]=coor[k,j]-grid[i,j]
        r[k,i,3]=np.sqrt(r[k,i,0]**2+r[k,i,1]**2+r[k,i,2]**2)
for k in range(n_atoms):
    for i in range(n_gridpoints):
        for j in range(3):
            X[k,i,j]=units*r[k,i,j]/(r[k,i,3]**3)

phi_ind_x=x-no
phi_ind_y=y-no
phi_ind_z=z-no



print('%6s %7s %7s %7s %10s' %("Name","a_xx", "a_yy", "a_zz", "a_iso"))
print('-----------------------------------------')
if n_atoms>1:
    mu_ind=np.zeros((n_atoms,3))
    mu_ind[:,0]=np.linalg.lstsq(X[:,:,0].T,phi_ind_x)[0]
    mu_ind[:,1]=np.linalg.lstsq(X[:,:,1].T,phi_ind_y)[0]
    mu_ind[:,2]=np.linalg.lstsq(X[:,:,2].T,phi_ind_z)[0]
    for i in range(n_atoms):
        a_xx=mu_ind[i,0]/field*0.529177249**2
        a_yy=mu_ind[i,1]/field*0.529177249**2
        a_zz=mu_ind[i,2]/field*0.529177249**2
        print('%6s %7.3f %7.3f %7.3f %10.3f' %(symbols[i],a_xx,a_yy,a_zz,(a_xx+a_yy+a_zz)/3))
else:
    mu_ind=np.zeros(3)
    mu_ind[0]=np.linalg.lstsq(np.vstack([X[0,:,0],np.zeros(n_gridpoints)]).T,phi_ind_x)[0][0]
    mu_ind[1]=np.linalg.lstsq(np.vstack([X[0,:,1],np.zeros(n_gridpoints)]).T,phi_ind_y)[0][0]
    mu_ind[2]=np.linalg.lstsq(np.vstack([X[0,:,2],np.zeros(n_gridpoints)]).T,phi_ind_z)[0][0]
    
    a_xx=mu_ind[0]/field*0.529177249**2
    a_yy=mu_ind[1]/field*0.529177249**2
    a_zz=mu_ind[2]/field*0.529177249**2

    print('%6s %7.3f %7.3f %7.3f %10.3f' %(symbols[0],a_xx,a_yy,a_zz,(a_xx+a_yy+a_zz)/3))



