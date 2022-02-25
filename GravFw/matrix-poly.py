#Este codigo lee una matriz con la estructura de densidades de un archivo.txt y genera un archivo de poligonos en el formato 
#que requiere el programa gravscript.m para utilizar el metodo de Talwani(1959) basado en el calculo de la 
#anomalia gravitacional causada por poligonos. Tambien genera un archivo con el arreglo en X de distancia en el 
#perfil en el formato que requiere gravscript.m

#De haber datos de gravimetria medidos en campo estos se pueden ingresar en la columna 2 del archivo Model_X.txt que se genera.

import numpy as np 
import matplotlib.pyplot as plt

#Lee matriz en pixeles para convertir a poligonos. Debe ingresar un archivo de texto que contiene la matriz de pixeles.
Datos=input("Ingrese el nombre del archivo con la estructura de densidades:")
X=input("Ingrese la distancia total del perfil en metros:")
Z=input("Ingrese la profundidad total en metros:")

#Numero de lineas de pixeles que no se quiere considerar en el calculo de la anomalia. 
#Por ejemplo si no se quiere considerar la corteza, se debe ingrasar la cantidad de filas que corresponden a corteza.
crustLines=int(float(input("Ingerse el numero de filas de la superficie que no se quieren considerar para el calculo de la anomalía:")))

#Lee los datos ingresados por input()
A=np.loadtxt(str(Datos))


#Imprime dimensiones de la estructura de densidades.
print('Distancia en x=', X, 'm', 'Distancia en z=', Z, 'm')


#Genera los arreglos de las distancias.
x=np.linspace(0,np.max(int(float(X))),len(A[0,:])+1)
z=np.linspace(0,np.max(int(float(Z))),len(A[:,0])+1)


#Grafica la estructura de densidades leida.
plt.pcolormesh(x,z,A)
plt.xlim(0,max(x))
plt.ylim(0,max(z))
plt.title('Estructura de densidades ingresada', fontsize=15)
plt.xlabel('$x(m)$', fontsize=15)
plt.ylabel('$z(m)$', fontsize=15)
clb = plt.colorbar()
clb.ax.set_title(r'$\rho$ $kg/m^{3}$')
plt.show()
plt.close()

#Calcula dx y dz. 
dx=x[1]-x[0]
dz=z[1]-z[0]

#Densidad de referencia para realizar el calculo de la anomalia gravitacional en gravscript.m
ref_rho=int(float(input("Ingrese la densidad de referencia para el calculo de la anomalía en kg/m^3:")))

#Numero de vertices por poligono mas uno (En gravscript.m siempre en cada poligono el primer vertice se ingresa dos veces, 
#una al inicio y una al final para cerrar el poligono). Aca los poligonos son pixeles de 4 vertices. 
vtx=5 

#Cambia unidades de distancia a kilometros. gravscript.m requiere distancias en kilometros.
x=x/1000
z=z/1000

#Estribe archivo con la informacion de lo poligonos que intresa a gravscript.m para el calculo de la anomalia
#gravitacional y lo guarda como "Model_mod.txt"
file=open('Model_mod.txt', 'w')
file.write(str(ref_rho)+'\n' )
file.write(str(np.size(A))+'\n' )
for i in range(len(A[:,0])):
    for j in range(len(A[0,:])):
        if(i>(len(A[:,0])-crustLines)):
            file.write(str(ref_rho)+'\n' )
        else:
            file.write(str(A[i,j])+'\n' )
        file.write(str(vtx)+'\n' )
        file.write('\t'+str(x[j])+' ' )
        file.write('\t'+str(max(z)-z[i])+'\n' )
        file.write('\t'+str(x[j+1])+' ' )
        file.write('\t'+str(max(z)-z[i])+'\n' )
        file.write('\t'+str(x[j+1])+' ' )
        file.write('\t'+str(max(z)-z[i+1])+'\n' )
        file.write('\t'+str(x[j])+' ' )
        file.write('\t'+str(max(z)-z[i+1])+'\n' )
        file.write('\t'+str(x[j])+' ' )
        file.write('\t'+str(max(z)-z[i])+'\n' )
file.close()
    

#Guarda arreglo de distancia x como Model_X.txt. gravscript.m recibe un archivo de 3 columnas donde
#la primera es el arreglo de distancia en x en metros. Las otras dos columnas se llenan con unos 
#porque no se van a usar.
MATRIX=np.ones((len(x),3))
MATRIX[:,0]=x*1000 #Se cambia unidades de x a metros 
np.savetxt('Model_X.txt', MATRIX, fmt='%f', delimiter='\t', newline='\r\n')