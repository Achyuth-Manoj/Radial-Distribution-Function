import numpy as np
import math
from matplotlib import pyplot as plt
#A POSCAR file input is required.
#THIS PROGRAM ONLY ACCEPTS BINARY COMPOUNDS
f = open("AgN_mp-1091417_Cmce_mmm.vasp","r")
#Input desired resolution as number of bins
bins=int(input("Enter desired number of bins: "))
#number of bins is an input from user
unit=float(10/bins)                        
coordinates_atom1=np.zeros(3)
coordinates_atom2=np.zeros(3)
lattice_parameters=[]
lattice_multipliers=[]
lattice_vectors=np.zeros(3)

#this loop for reading the file line by line
#it reads the coordinates of each atom and makes a 2D list from that 
#takes into account whether coordinates given are direct or cartesian
#if direct program converts to cartesian
#also finds out how many atoms of each type there are
line=f.readline()
l=1
while line:
    l=l+1
    line=f.readline()
    if l<=2:
        continue
    elif 2<l<=5:
        z=line.split()
        lattice_vect=[]
        for w in z:
            lattice_vect.append(float(w))
        parameter=math.sqrt((lattice_vect[0]**2)+(lattice_vect[1]**2)+(lattice_vect[2]**2)) #this line calculates each lattice parameter
        lattice_parameters.append(float(parameter))
        lattice_vectors=np.vstack((lattice_vectors,lattice_vect))
        if parameter<1:
            m=20
        elif 1<=parameter<2:
            m=10
        elif 2<=parameter<3:
            m=5
        elif 3<parameter<=4:
            m=4
        elif 4<parameter<=5:
            m=3
        elif 5<=parameter<10:
            m=2
        else:
            m=1
        lattice_multipliers.append(float(m))
    elif l==6:
        z=line.split()
    elif l==7:
        a=line.split()
        atoms=[]
        for w in a:
            atoms.append(float(w))
    elif l==8:
        a=line.split()
        s=a[0]
        C=["Cartesian"]
        D=["Direct"]
    elif 8<l<=int(8+atoms[0]):
        x=line.split()
        if s in C:
             #cartesian coordinates directly added to list
            floats=[]
            for y in x:
                m=float(y)
                floats.append(float(m))
        elif s in D:
            #direct coordinates converted to cartesian
            r=1
            floats=np.zeros(3)
            for y in x:
                m=np.multiply(lattice_vectors[int(r)],float(y))
                floats=np.add(floats,m)
                r=r+1
        coordinates_atom1=np.vstack((coordinates_atom1,floats))
    elif int(8+atoms[0])<l<=int(8+sum(atoms)):
        x=line.split()
        if s in C:
            #cartesian coordinates directly added to list
            floats=[]
            for y in x:
                m=float(y)
                floats.append(float(m))
        elif s in D:
            #direct coordinates converted to cartesian
            r=0 
            floats=np.zeros(3)
            for y in x:
                m=np.multiply(lattice_vectors[int(r+1)],float(y))
                floats=np.add(floats,m)
                r=r+1
        coordinates_atom2=np.vstack((coordinates_atom2,floats))


f.close()
#removes the first row of zeros which is not part of the coordinates
coordinates_atom1=np.delete(coordinates_atom1,0,axis=0)
coordinates_atom2=np.delete(coordinates_atom2,0,axis=0)
lattice_vectors=np.delete(lattice_vectors,0,axis=0)
#The program prints the basis vectors of the lattice and their magnitudes
#Also prints the multiplication factors for creating lesser lattice
print(lattice_vectors)
print(lattice_parameters)
print(lattice_multipliers)

#unit cell multiplier, creating lesser lattice
for p in [0,1,2]:
    m=lattice_multipliers[p]
    #the atoms whose coorinates are in the folllowing two lists are the ones
    #for which rdf is calculated
    coordinates_atom1_lesser=np.zeros(3)
    coordinates_atom2_lesser=np.zeros(3)
    # multiplying the coordinates of the first atom
    for x in coordinates_atom1:
        n=0
        a=x
        while n<m:
            if n==0:
                coordinates_atom1_lesser=np.vstack((coordinates_atom1_lesser,a))
            else:
                #creating a new particle position
                a=np.add(a,lattice_vectors[p])
                coordinates_atom1_lesser=np.vstack((coordinates_atom1_lesser,a))
            n=n+1
    
    # multiplying the coordinates of the second atom
    for y in coordinates_atom2:
        n=0
        a=y
        while n<m:
            if n==0:
                coordinates_atom2_lesser=np.vstack((coordinates_atom2_lesser,a))
            else:
                #creating a new particle position
                a=np.add(a,lattice_vectors[p])
                coordinates_atom2_lesser=np.vstack((coordinates_atom2_lesser,a))
            n=n+1
    coordinates_atom1_lesser=np.delete(coordinates_atom1_lesser,0,axis=0)
    coordinates_atom2_lesser=np.delete(coordinates_atom2_lesser,0,axis=0)
    coordinates_atom1=coordinates_atom1_lesser
    coordinates_atom2=coordinates_atom2_lesser

#creating a list containing the number of atoms of each type for lesser lattice
#and stacking the coordinates for atom1 and atom2
atoms_lesser=[0,0]
atoms_lesser[0]=np.shape(coordinates_atom1_lesser)[0]
atoms_lesser[1]=np.shape(coordinates_atom2_lesser)[0]
coordinates_lesser=np.vstack((coordinates_atom1_lesser,coordinates_atom2_lesser))

#updating basis vectors of the lesser lattice
lattice_parameters_lesser=[0,0,0]
lattice_vectors_lesser=np.zeros(9).reshape((3,3))
for p in [0,1,2]:
    lattice_parameters_lesser[p]=lattice_parameters[p]*lattice_multipliers[p]
    lattice_vectors_lesser[p]=lattice_vectors[p]*lattice_multipliers[p]
print(lattice_vectors_lesser)

#lesser lattice multiplier, creating the greater lattice
for p in [0,1,2]:
    coordinates_atom1_greater=np.zeros(3)
    coordinates_atom2_greater=np.zeros(3)
    # multiplying the coordinates of the first atom
    for x in coordinates_atom1:
        n=0
        a=x
        while n<3:
            if n==0:
               coordinates_atom1_greater=np.vstack((coordinates_atom1_greater,a))
            elif n==1:
                 #creating a new particle position in the negative direction
                a=np.subtract(x,lattice_vectors_lesser[p])
                coordinates_atom1_greater=np.vstack((coordinates_atom1_greater,a))
            else:
                #creating a new particle position in positive direction
                a=np.add(x,lattice_vectors_lesser[p])
                coordinates_atom1_greater=np.vstack((coordinates_atom1_greater,a))
            n=n+1
    
    # multiplying the coordinates of the second atom
    for y in coordinates_atom2:
        n=0
        a=y
        while n<3:
            if n==0:
               coordinates_atom2_greater=np.vstack((coordinates_atom2_greater,a))
            elif n==1:
                 #creating a new particle position in the negative direction
                a=np.subtract(y,lattice_vectors_lesser[p])
                coordinates_atom2_greater=np.vstack((coordinates_atom2_greater,a))
            else:
                #creating a new particle position in positive direction
                a=np.add(y,lattice_vectors_lesser[p])
                coordinates_atom2_greater=np.vstack((coordinates_atom2_greater,a))
            n=n+1
    
    coordinates_atom1_greater=np.delete(coordinates_atom1_greater,0,axis=0)
    coordinates_atom2_greater=np.delete(coordinates_atom2_greater,0,axis=0)
    coordinates_atom1=coordinates_atom1_greater
    coordinates_atom2=coordinates_atom2_greater

#creating a list containing the number of atoms of each type for greater lattice
#and stacking the coordinates for atom1 and atom2
atoms_greater=[0,0]
atoms_greater[0]=np.shape(coordinates_atom1_greater)[0]
atoms_greater[1]=np.shape(coordinates_atom2_greater)[0]
coordinates_greater=np.vstack((coordinates_atom1_greater,coordinates_atom2_greater))

#updating basis vectors of the greater lattice, not a necessary step
lattice_parameters_greater=[0,0,0]
lattice_vectors_greater=np.zeros(9).reshape((3,3))
for p in [0,1,2]:
    lattice_parameters_greater[p]=lattice_parameters_lesser[p]*3
    lattice_vectors_greater[p]=lattice_vectors_lesser[p]*3

#saving the coordinates of the lesser and greater lattice in separate files
f=open("coordinates_lesser.txt","w")
np.savetxt(f,coordinates_lesser,delimiter=' ')
f.close()
f=open("coordinates_greater.txt","w")
np.savetxt(f,coordinates_greater,delimiter=' ')
f.close()

#Finding global density of structure
#global volume is the volume of the lesser lattice
#number of atoms considered is that of the lesser lattice
N_global=sum(atoms_lesser)
u=np.cross(lattice_vectors_lesser[1],lattice_vectors_lesser[2])
V_global=np.dot(lattice_vectors_lesser[0],u)
V_global_mag=np.linalg.norm([V_global])
Density_global=N_global/V_global_mag
print("Volume of lesser lattice: " ,V_global_mag)
print("Number of atoms in lesser lattice: ",N_global)
print("Number of atoms in greater lattice: ",sum(atoms_greater))

#Calculation of Radial Distribution function for all atoms in lesser lattice
n=0
#rdf_sum will be the added up rdf for this calculation
rdf_sum=np.zeros((bins,1))
rdf_cum_total=np.zeros((bins,1))
count_sum=np.zeros((bins,1))

for a in coordinates_lesser:
    
    a=coordinates_lesser[int(n)]
    calculations=[]
    #The following loop calculates the distances of all other atoms in the greater lattice
    #from the selected lattice
    for c in coordinates_greater:
        b=np.subtract(a,c)
        d=math.sqrt((b[0]**2)+(b[1]**2)+(b[2]**2))
        calculations.append(float(d))
    #removing the distance of selected atom to itself (which is zero)
    calculations.remove(calculations[27*int(n)])
    n=n+1
    
    radius=0
    rdf_data=np.zeros(1)
    rdf_cum=np.zeros(1)
    #first column of our results will be the radius, which is why radii is 
    #declared only in this rdf calculation and not the following ones, since
    #it is only needed once for recording the radii
    radii=np.zeros(1)
    #j is cumulative atom count
    j=0
    count=np.zeros(1)
    
    #This loop calculates rdf for the single atom
    for f in np.arange(0,10,unit):
        #l is the number of atoms in the shell
        l=0
        for h in calculations:
            if (radius-(unit/2))<=h<=(radius+(unit/2)):
                #this condition is used to count the atoms in a given shell
                l=l+1
                j=j+1
            else:
                continue
        count=np.vstack((count,[l]))
        V_local=(4/3)*(np.pi)*(((radius)+(unit/2))**3-((radius)-(unit/2))**3)
        Density_local=(l)/(V_local)
        rdf_calc=(Density_local)/(Density_global)
        rdf_cum=np.vstack((rdf_cum,[j]))
         #stacking each count
        rdf_data=np.vstack((rdf_data,[rdf_calc]))
        radii=np.vstack((radii,[float(radius)]))
        #increasing radius by unit
        radius=radius + unit
    
    #these two lines remove first row of zeros
    radii=np.delete(radii,0,axis=0)
    count=np.delete(count,0,axis=0)
    rdf_cum=np.delete(rdf_cum,0,axis=0)
    rdf_data=np.delete(rdf_data,0,axis=0)
    rdf_cum_total=np.add(rdf_cum_total,rdf_cum)
    #adding to rdf
    rdf_sum=np.add(rdf_sum,rdf_data)
    count_sum=np.add(count_sum,count)
    
result=radii
#finding rdf from sum_rdf by dividing number of atoms
rdf=np.divide(rdf_sum,float(sum(atoms_lesser)))
rdf_cumulative=np.divide(rdf_cum_total,float(sum(atoms_lesser))) 
count_average=np.divide(count_sum,float(sum(atoms_lesser)))
#this command stacks the columns together
result=np.hstack((result,rdf,rdf_cumulative,count_average))
#plotting the graphs with these two lines, haven't given them titles yet, 
#so the graphs are printed in order of the columns of the result 
plt.plot(radii,rdf)                               
#plt.plot(radii,rdf_cumulative)
plt.show()

#Calculation of Radial Distribution function for cations-cations
n=0
#rdf_sum will be the added up rdf for this calculation
rdf_sum=np.zeros((bins,1))
rdf_cum_total=np.zeros((bins,1))
count_sum=np.zeros((bins,1))
for a in coordinates_atom1_lesser:
    
    a=coordinates_atom1_lesser[int(n)]
    calculations=[]
    #The following loop calculates the distances of all other atoms in the greater lattice
    #from the selected lattice
    for c in coordinates_atom1_greater:
        b=np.subtract(a,c)
        d=math.sqrt((b[0]**2)+(b[1]**2)+(b[2]**2))
        calculations.append(float(d))
    #removing the distance of selected atom to itself (which is zero)
    calculations.remove(calculations[27*int(n)])
    n=n+1
    
    radius=0
    rdf_data=np.zeros(1)
    rdf_cum=np.zeros(1)
    #first column of our results will be the radius, which is why radii is 
    #declared only in this rdf calculation and not the following ones, since
    #it is only needed once for recording the radii
    radii=np.zeros(1)
    #j is cumulative atom count
    j=0
    count=np.zeros(1)
    
    #This loop calculates rdf for the single atom
    for f in np.arange(0,10,unit):
        #l is the number of atoms in the shell
        l=0
        for h in calculations:
            if (radius-(unit/2))<=h<=(radius+(unit/2)):
                #this condition is used to count the atoms in a given shell
                l=l+1
                j=j+1
            else:
                continue
        count=np.vstack((count,[l]))
        V_local=(4/3)*(np.pi)*(((radius)+(unit/2))**3-((radius)-(unit/2))**3)
        Density_local=(l)/(V_local)
        rdf_calc=(Density_local)/(Density_global)
        rdf_cum=np.vstack((rdf_cum,[j]))
         #stacking each count
        rdf_data=np.vstack((rdf_data,[rdf_calc]))
        radii=np.vstack((radii,[float(radius)]))
        #increasing radius by unit
        radius=radius + unit
    
    #these two lines remove first row of zeros
    radii=np.delete(radii,0,axis=0)
    count=np.delete(count,0,axis=0)
    rdf_cum=np.delete(rdf_cum,0,axis=0)
    rdf_data=np.delete(rdf_data,0,axis=0)
    rdf_cum_total=np.add(rdf_cum_total,rdf_cum)
    #adding to rdf
    rdf_sum=np.add(rdf_sum,rdf_data)
    count_sum=np.add(count_sum,count)

#finding rdf from sum_rdf by dividing number of atoms
rdf=np.divide(rdf_sum,float(atoms_lesser[0]))
rdf_cumulative=np.divide(rdf_cum_total,float(atoms_lesser[0])) 
count_average=np.divide(count_sum,float(atoms_lesser[0]))
#this command stacks the columns together
result=np.hstack((result,rdf,rdf_cumulative,count_average))
#plotting the graphs with these two lines, haven't given them titles yet, 
#so the graphs are printed in order of the columns of the result 
plt.plot(radii,rdf)                               
#plt.plot(radii,rdf_cumulative)
plt.show()

#Calculation of Radial Distribution Function for cation-anion
n=0
#rdf_sum will be the added up rdf for this calculation
rdf_sum=np.zeros((bins,1))
rdf_cum_total=np.zeros((bins,1))
count_sum=np.zeros((bins,1))
for a in coordinates_atom1_lesser:
    
    a=coordinates_atom1_lesser[int(n)]
    calculations=[]
    #The following loop calculates the distances of all other atoms in the greater lattice
    #from the selected lattice
    for c in coordinates_atom2_greater:
        b=np.subtract(a,c)
        d=math.sqrt((b[0]**2)+(b[1]**2)+(b[2]**2))
        calculations.append(float(d))
    n=n+1
    
    radius=0
    rdf_data=np.zeros(1)
    rdf_cum=np.zeros(1)
    #first column of our results will be the radius, which is why radii is 
    #declared only in this rdf calculation and not the following ones, since
    #it is only needed once for recording the radii
    radii=np.zeros(1)
    #j is cumulative atom count
    j=0
    count=np.zeros(1)
    
    #This loop calculates rdf for the single atom
    for f in np.arange(0,10,unit):
        #l is the number of atoms in the shell
        l=0
        for h in calculations:
            if (radius-(unit/2))<=h<=(radius+(unit/2)):
                #this condition is used to count the atoms in a given shell
                l=l+1
                j=j+1
            else:
                continue
        count=np.vstack((count,[l]))
        V_local=(4/3)*(np.pi)*(((radius)+(unit/2))**3-((radius)-(unit/2))**3)
        Density_local=(l)/(V_local)
        rdf_calc=(Density_local)/(Density_global)
        rdf_cum=np.vstack((rdf_cum,[j]))
         #stacking each count
        rdf_data=np.vstack((rdf_data,[rdf_calc]))
        radii=np.vstack((radii,[float(radius)]))
        #increasing radius by unit
        radius=radius + unit
    
    #these two lines remove first row of zeros
    radii=np.delete(radii,0,axis=0)
    count=np.delete(count,0,axis=0)
    rdf_cum=np.delete(rdf_cum,0,axis=0)
    rdf_data=np.delete(rdf_data,0,axis=0)
    rdf_cum_total=np.add(rdf_cum_total,rdf_cum)
    #adding to rdf
    rdf_sum=np.add(rdf_sum,rdf_data)
    count_sum=np.add(count_sum,count)

#finding rdf from sum_rdf by dividing number of atoms
rdf=np.divide(rdf_sum,float(atoms_lesser[0]))
rdf_cumulative=np.divide(rdf_cum_total,float(atoms_lesser[0])) 
count_average=np.divide(count_sum,float(atoms_lesser[0]))
#this command stacks the columns together
result=np.hstack((result,rdf,rdf_cumulative,count_average))
#plotting the graphs with these two lines, haven't given them titles yet, 
#so the graphs are printed in order of the columns of the result 
plt.plot(radii,rdf)                               
#plt.plot(radii,rdf_cumulative)
plt.show()

#Calculation of Radial Distribution Function for Anion-anion
n=0
#rdf_sum will be the added up rdf for this calculation
rdf_sum=np.zeros((bins,1))
rdf_cum_total=np.zeros((bins,1))
count_sum=np.zeros((bins,1))
for a in coordinates_atom2_lesser:
    
    a=coordinates_atom2_lesser[int(n)]
    calculations=[]
    #The following loop calculates the distances of all other atoms in the greater lattice
    #from the selected lattice
    for c in coordinates_atom2_greater:
        b=np.subtract(a,c)
        d=math.sqrt((b[0]**2)+(b[1]**2)+(b[2]**2))
        calculations.append(float(d))
    #removing the distance of selected atom to itself (which is zero)
    calculations.remove(calculations[27*int(n)])
    n=n+1
    
    radius=0
    rdf_data=np.zeros(1)
    rdf_cum=np.zeros(1)
    #first column of our results will be the radius, which is why radii is 
    #declared only in this rdf calculation and not the following ones, since
    #it is only needed once for recording the radii
    radii=np.zeros(1)
    #j is cumulative atom count
    j=0
    count=np.zeros(1)
    
    #This loop calculates rdf for the single atom
    for f in np.arange(0,10,unit):
        #l is the number of atoms in the shell
        l=0
        for h in calculations:
            if (radius-(unit/2))<=h<=(radius+(unit/2)):
                #this condition is used to count the atoms in a given shell
                l=l+1
                j=j+1
            else:
                continue
        count=np.vstack((count,[l]))
        V_local=(4/3)*(np.pi)*(((radius)+(unit/2))**3-((radius)-(unit/2))**3)
        Density_local=(l)/(V_local)
        rdf_calc=(Density_local)/(Density_global)
        rdf_cum=np.vstack((rdf_cum,[j]))
         #stacking each count
        rdf_data=np.vstack((rdf_data,[rdf_calc]))
        radii=np.vstack((radii,[float(radius)]))
        #increasing radius by unit
        radius=radius + unit
    
    #these two lines remove first row of zeros
    radii=np.delete(radii,0,axis=0)
    count=np.delete(count,0,axis=0)
    rdf_cum=np.delete(rdf_cum,0,axis=0)
    rdf_data=np.delete(rdf_data,0,axis=0)
    rdf_cum_total=np.add(rdf_cum_total,rdf_cum)
    #adding to rdf
    rdf_sum=np.add(rdf_sum,rdf_data)
    count_sum=np.add(count_sum,count)

#finding rdf from sum_rdf by dividing number of atoms
rdf=np.divide(rdf_sum,float(atoms_lesser[1]))
rdf_cumulative=np.divide(rdf_cum_total,float(atoms_lesser[1])) 
count_average=np.divide(count_sum,float(atoms_lesser[1]))
#this command stacks the columns together
result=np.hstack((result,rdf,rdf_cumulative,count_average))
#plotting the graphs with these two lines, haven't given them titles yet, 
#so the graphs are printed in order of the columns of the result 
plt.plot(radii,rdf)                               
#plt.plot(radii,rdf_cumulative)
plt.show()

print(result)
f=open("rdf.txt","w")
np.savetxt(f,result,delimiter=' ')
column_names=["radius", "rdf for all atoms", "cumulative rdf for all atoms", "rdf for cation-cation", "cumulative rdf for cation-cation", "rdf for cation-anion", "cumulative rdf for cation-anion","rdf for anion-cation","cumulative rdf for anion-cation", "rdf for anion-anion", "cumulative rdf for anion-anion"]
#result=np.vstack((column_names,result))
#g=open("rdf_with_headers.txt","w")
#np.savetxt(g,result,delimiter=' ')
f.close()
#g.close()

