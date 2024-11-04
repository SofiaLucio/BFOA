import math
import random
from bacteria import bacteria

class chemiotaxis():
    def __init__(self):
        self.parcialNFE = 0  
        self.base_d_attr = 0.1
        self.base_w_attr = 0.2
        self.base_h_rep = 0.1
        self.base_w_rep = 10
        self.initial_population_size = 10  # Aumentar el tamaño de la población inicial

    def compute_cell_interaction(self, bacteria, poblacion, d, w):
        total = 0.0
        for other in poblacion:
            diff = (bacteria.blosumScore - other.blosumScore) ** 2.0
            total += d * math.exp(w * diff)
        return total

    def attract_repel(self, bacteria, poblacion, d_attr, w_attr, h_rep, w_rep):
        attract = self.compute_cell_interaction(bacteria, poblacion, -d_attr, -w_attr)
        repel = self.compute_cell_interaction(bacteria, poblacion, h_rep, -w_rep)
        return attract + repel

    def chemio(self, bacteria, poblacion, d_attr, w_attr, h_rep, w_rep):
        bacteria.interaction = self.attract_repel(bacteria, poblacion, d_attr, w_attr, h_rep, w_rep)
        bacteria.fitness = bacteria.blosumScore + bacteria.interaction

    def doChemioTaxis(self, poblacion, d_attr, w_attr, h_rep, w_rep):
        self.parcialNFE = 0
        for bacteria in poblacion:
            self.chemio(bacteria, poblacion, d_attr, w_attr, h_rep, w_rep)
            self.parcialNFE += bacteria.NFE
            bacteria.NFE = 0

    def eliminarClonar(self, path, poblacion):
        poblacion.sort(key=lambda x: x.fitness)
        del poblacion[:int(len(poblacion) * 0.6)]  # Elimina el 60% de las bacterias
        clones = self.clonacion(path, poblacion)
        poblacion.extend(clones)

    def clonacion(self, path, poblacion):
        poblacionClones = []
        best = max(poblacion, key=lambda x: x.fitness)
        for bacteria in poblacion:
            newBacteria = bacteria.clonar(path)
            mutacion = int((best.fitness - bacteria.fitness) / 10)  # Mutación adaptativa
            newBacteria.tumboNado(mutacion)
            newBacteria.autoEvalua()
            poblacionClones.append(newBacteria)
        return poblacionClones

    def randomBacteria(self, path):
        bact = bacteria(path)
        bact.tumboNado(random.randint(1, 10))
        bact.autoEvalua()
        return bact 

    def insertRamdomBacterias(self, path, num, poblacion):
        for _ in range(num):
            new_bact = self.randomBacteria(path)
            poblacion.append(new_bact)
            poblacion.sort(key=lambda x: x.fitness, reverse=True)
            poblacion.pop()  # Elimina la bacteria con menor fitness

    def adjust_interaction_parameters(self, iteracion_actual, max_iteraciones):
        self.d_attr = self.base_d_attr + (0.05 * iteracion_actual / max_iteraciones)
        self.w_attr = self.base_w_attr + (0.05 * iteracion_actual / max_iteraciones)
        self.h_rep = self.base_h_rep - (0.05 * iteracion_actual / max_iteraciones)
        self.w_rep = self.base_w_rep - (0.5 * iteracion_actual / max_iteraciones)
        return self.d_attr, self.w_attr, self.h_rep, self.w_rep

