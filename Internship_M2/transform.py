# Este es un ejemplo de como se transforman datos del codigo Axisimetrico: Ollin-axis a coordenadas tridimensionales 
# Utilizando geometria cartesiana y las propiedades simetricas del sistema
import numpy as np #Se importa la libreria numpy
import matplotlib.pyplot as plt # Se importa la libreria matplotlib.pyplot para checar que las transformaciones sean correctas
import h5py # Se importa la libreria h5py para guardar los datos en un archivo .h5 (Tipo de archivos usados por Simflorny)

Data_phi = np.loadtxt("phi_b0l0.2D") # Leemos el archivo y lo guardamos en un arreglo llamado Data_phi
# phi_b0l0.2D es el Nombre del archivo a importar, se tiene que generar primero con Ollin Axis.
# En este ejemplo, este archivo corresponde a los datos iniciales del logaritmo del factor conforme para datos de tamano: 50 puntos
# en la coordenada x y 50 puntos en la coordenada z. Tambien se tienen tres puntos adicionales en cada coordenada por los puntos fantasma
# Hay que precisar que hay que conocer bien la estructura de los datos para poder leerlos.
# En este caso, el archivo consta de 53x53 lineas y tres columnas. Las columnas corresponde respectivamente a la coordenada en r, la coordenada en z y el valor del
# Logaritmo del factor conforme.
Phi0 = Data_phi[:,2] # Guardamos los datos de Phi (el logaritmo del factor conforme) en un arreglo unidimencional
Phi0dr = np.reshape(Phi0,(53,53)) # Conociendo la estructura de los datos, convertimos el arreglo de los valores de Phi en una matriz bidimensional
#Donde cada columna representa una coordenada en r y cada linea en z. De esta forma se tiene un campo de datos bidimensionales.
Phi0dr = np.delete(Phi0dr,slice(3),axis=0) # Eliminamos las tres primeras lineas (Datos de los puntos fantasma)
Phi0dr = np.delete(Phi0dr,slice(3),axis=1) # Eliminamos las tres primeras columnas (Datos de los puntos fantasma)
# Phi0dr sera el arreglo que crooespondera al cuadrante inferior derecho


##################### Simetria de revolucion ######################## 
# ESte procedimiento utiliza simetrias de refleccion en geometria tridimensional por lo que es un procedimiento que utiliza muchos recursos numericos
# Para optimizar el codigo, seria necesario generar primero todo el plano (y=0) en dos dimensiones y luego aplicarle la simetria de revolucion
Phi3Dur1 = np.zeros((50,50,50)) #Generamos un arreglo tridimensional de con 50 0' en cada eje
for i in range(50): # Nos desplazamos sobre z en el arreglo dimensional
    for l in range(50): # Hacemos un ciclo sobre los valores de r del arreglo bidimensional
        for j in range(i+1): # Hacemos un ciclo sobre los valores de y y x del arreglo tridimensional
            for k in range(i+1):
                if k==i or j==i: # Por analisis, si la posicion de x o la de y del arreglo tridimensional es igual a la posicion de x del arreglo 
                    #bidimensional, adjuntamos el valor, y repetimos esto para cada z creado revoluciones artificiales por niveles            
                    Phi3Dur1[j,k,l]=Phi0dr[l,i]
# Este procedimiento genera los datos para el octante superior izquierdo de frente
#plt.imshow(Phi3Dur1[:,:,0])
#plt.show()
Phi3Dur2 = np.flip(Phi3Dur1,axis=1) # Invertimos el arreglo para genera los datos restantes del octante superior derecho
Phi3Dur=np.concatenate((Phi3Dur2,Phi3Dur1),axis=1) # Concatenamos para unir las dos partes
Phi3Ddr=np.flip(Phi3Dur,axis=2) # Invertimos para obtener los octantes inferiores derechos
Phi3Dr=np.concatenate((Phi3Ddr,Phi3Dur), axis=2) # Concatenamos para obtener toda la parte inferior del espacio
Phi3Dl=np.flip(Phi3Dr,axis=0) # Invertimos para obtener la parte superior
Phi3D=np.concatenate((Phi3Dl,Phi3Dr),axis=0) # Concatenamos para obtener el espacio entero

"""
Este es un test para ver como Simflowny leia los datos. Se genero un arreglo conocido de datos, se hizo el primer paso de simulacion
y se analizo como se inicializaban los datos de Phi.
Phitest=np.zeros((100,100,100))
for k in range(100):
    for j in range(100):
        for i in range(100):
            Phitest[i,j,k] = i+100*(j+100*k)
"""
# Arreglamos el arreglo para que tenga la forma correcta puesto que algunos octantes estan invertidos
# Tambien arreglamos el arreglo de tal manera que el primer indice corresponda a x, el segundo a y y el tercero a z 
for i in range(100): 
    Phi3D[:,:,i]= Phi3D[:,:,i].T
Phi3D=np.swapaxes(Phi3D,1,2)
Phi3D=np.swapaxes(Phi3D,0,1)
hf = h5py.File('phi.h5', 'w') # Creamos el archivo donde guardar los datos
hf.create_dataset("Phi", data=Phi3D) # Guardamos los datos en el archivo
hf.close() # Cerramos el archivo
