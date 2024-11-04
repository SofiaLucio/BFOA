from bacteria import bacteria
from chemiotaxis import chemiotaxis
import numpy

poblacion = []
path = "C:\\secuenciasBFOA\\multiFasta.fasta"  # Asegúrate de usar dobles barras invertidas
numeroDeBacterias = 5
numRandomBacteria = 1
iteraciones = 30
tumbo = 1  # Número de gaps a insertar 
nado = 3
chemio = chemiotaxis()
veryBest = bacteria(path)  # Mejor bacteria   
tempBacteria = bacteria(path)  # Bacteria temporal para validaciones
original = bacteria(path)  # Bacteria original sin gaps
globalNFE = 0  # Número de evaluaciones de la función objetivo

dAttr = 0.1  # 0.1
wAttr = 0.2  # 0.2
hRep = dAttr
wRep = 10    # 10


def clonaBest(veryBest, best):
    veryBest.matrix.seqs = numpy.array(best.matrix.seqs)
    veryBest.blosumScore = best.blosumScore
    veryBest.fitness = best.fitness
    veryBest.interaction = best.interaction

def validaSecuencias(path, veryBest):
    # Clona a veryBest en tempBacteria   
    tempBacteria.matrix.seqs = numpy.array(veryBest.matrix.seqs)
    # Descartar los gaps de cada secuencia
    for i in range(len(tempBacteria.matrix.seqs)):
        tempBacteria.matrix.seqs[i] = tempBacteria.matrix.seqs[i].replace("-", "")
    
    # Valida que las secuencias originales sean iguales a las secuencias de tempBacteria
    for i in range(len(tempBacteria.matrix.seqs)):
        if tempBacteria.matrix.seqs[i] != original.matrix.seqs[i]:
            print("*****************Secuencias no coinciden********************")
            return


# Inicializar población
for i in range(numeroDeBacterias):  # Población inicial
    poblacion.append(bacteria(path))

# Bucle de iteraciones
for iteracion in range(iteraciones):
    # Ajustar parámetros de interacción
    d_attr, w_attr, h_rep, w_rep = chemio.adjust_interaction_parameters(iteracion, iteraciones)
    
    chemio.doChemioTaxis(poblacion, d_attr, w_attr, h_rep, w_rep)
    
    globalNFE += chemio.parcialNFE  # Actualiza NFE
    best = max(poblacion, key=lambda x: x.fitness)  # Encuentra la mejor bacteria

    # Clonación de la mejor bacteria
    if (veryBest is None) or (best.fitness > veryBest.fitness):
        clonaBest(veryBest, best)

    # Imprimir resultados de la iteración
    print(f"Iteracion: {iteracion + 1}, Mejor fitness: {veryBest.fitness}, NFE: {globalNFE}")
    
    # Eliminar y clonar bacterias
    chemio.eliminarClonar(path, poblacion)

    # Insertar bacterias aleatorias cada 5 iteraciones
    if iteracion % 5 == 0:  
        chemio.insertRamdomBacterias(path, numRandomBacteria, poblacion)


veryBest.showGenome()
validaSecuencias(path, veryBest)

