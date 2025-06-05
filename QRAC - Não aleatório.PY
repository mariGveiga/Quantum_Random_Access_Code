# --------------- QRAC sem a aleatoriedade --- somente os estados da base computacionl ---------------

from qiskit import QuantumCircuit, transpile
from qiskit_aer import AerSimulator
import numpy as np
import math as math

simulator = AerSimulator()
d = 2           # words size 

# function to create the superposition of the qubit so Alice could send it to Bob
def create_states(lista):
    x0 = lista[0]
    x1 = lista[1]
       
    amplitudes = []

    # creating amplitudes according to the formula
    ampl_1 = 1/math.sqrt(2)     # phi = 0 -- variação ao longo do eixo longitudinal
    amplitudes.append(ampl_1)

    ampl_2 = ((-1)**x0 +1j*(-1)**x1)/2      # phi = 45 (maximizado) 
    amplitudes.append(ampl_2)

    return(amplitudes)


def initialize_circuit(ampl_list):
    qc = QuantumCircuit(1,1)
    qc.initialize(ampl_list,0)      # initialize the circuit with the qubit Alice sent 
    qc.save_statevector()

    # Transpile for simulator
    qc_aer = transpile(qc, backend=simulator)

    result = simulator.run(qc_aer).result()
    #psi = result.get_statevector()

    return qc

# create psi states 
states = []
for i in range(d):
    linha = []
    for j in range(d):
        linha.append([i,j])
    states.append(linha)

print("Os estados gerados para a função psi são:", states)

# gera uma amplitude para cada estado de função criado em 'states'
ampl_states = []
for i in range(d):
    linha = []
    for j in range(d):
        state = create_states(states[i][j])
        linha.append(state)
    ampl_states.append(linha)

all_circuits = []           # matriz que vai armazenar todas as informações de todos os processamentos de cada estado da função psi
for i in range(d):
    linha = []
    for j in range(d):
        circ_ampl = ampl_states[i][j]
        # print("Amplitudes:", circ_ampl)
        circ_simulation = initialize_circuit(circ_ampl)         # inicializa o circuito para cada amplitude dos estados

        indiv = []              # cria uma linha individual para cada estado 

        indiv.append(states[i][j])
        indiv.append(circ_ampl)
        
        for indice in range(d):     # gera a medição tanto para o índice 0 quanto para o índice 1 para todas as amplitudes
            if indice == 0:         # primeira letra
                indiv.append("O qubit foi medido na base sigma x")
                circ_simulation.h(0)
            else:                   # segunda letra 
                indiv.append("O qubit foi medido na base sigma y")
                circ_simulation.sdg(0)
                circ_simulation.h(0)

            # realização da medição e da simulação de cada circuito 
            circ_simulation.measure(0,0)
            result = simulator.run(circ_simulation, shots = 1024).result()
            counts = result.get_counts()
            indiv.append(f"{counts}")

            # mapear o resultado (0 or 1) para autovalores (+1 and -1)
            eigenvalue_mapping = {'0': +1, '1': -1}
            eigenvalues = {eigenvalue_mapping[bitstring]: count for bitstring, count in counts.items()}

            eigenvalue_counts = {+1: 0, -1: 0}
            for eigenvalue, count in eigenvalues.items():
                eigenvalue_counts[eigenvalue] += count

            # colocar os autovalores individuais de cada amplitude
            indiv.append(eigenvalue_counts)

            # encontrar qual é o máximo do autovalor 
            max_eigenvalue = max(eigenvalue_counts, key=eigenvalue_counts.get)
            indiv.append(f"{max_eigenvalue}")

            
        linha.append(indiv)
    all_circuits.append(linha)

# esse for é apenas para conseguir checar as informações para ver se estão coerentes 
# for i in range(d):
#     for j in range(d):        
#         item = all_circuits[i][j]
#         print("O estado que está sendo analisado da função psi é:", item[0])
#         print("Sua respectiva amplitude é:", item[1])

#         print("---- Primeira Simulação ----")
#         print(item[2])
#         print("Resultado da simulação:", item[3])
#         print("Cálculo dos autovalores:", item[4])
#         print("Maior autovalor:", item[5])

#         print("---- Segunda Simulação ----")
#         print(item[6])
#         print("Resultado da simulação:", item[7])
#         print("Cálculo dos autovalores:", item[8])
#         print("Maior autovalor:", item[9])

#         print()
 
correct_results = []            # matriz que vai armazenar apenas os resultados em que Bob conesegue acertar a letra da Alice 
for i in range(d):
    for j in range(d):
        item = all_circuits[i][j]
        s1 = int(item[5])        # autovalor da simulação 1 -- índice = 0
        s2 = int(item[9])        # autovalor da simulação 2 -- índice = 1

        # calculando os valores respectivos em binário dos autovalores: -1 = 1, 1 = 0
        if s1 == -1:
            s1_bin = 1
        else:
            s1_bin = 0        
        if s2 == -1:
            s2_bin = 1
        else:
            s2_bin = 0

        for indice in range(d):            
            letra_alice = states[i][j][indice]      # pega letra da Alice, índice por índice 
            linha = []
            if indice == 0:
                if s1_bin == letra_alice:
                    # print("Certo: %d %d, indice 0" %(s1_bin, letra_alice))
                    linha.append(indice)                    # qual é o índice que foi analisado
                    linha.append(item[2])                   # em qual base foi medida
                    linha.append(letra_alice)               # qual foi a letra enviada por Alice 
                    linha.append(item[5])                   # autovalor
                else:
                    # print("Incorreto: %d %d, indice 0" %(s1_bin, letra_alice))
                    pass
            else:
                if s2_bin == letra_alice:
                    # print("Certo: %d %d, indice 1" %(s2_bin, letra_alice))
                    linha.append(indice)                    # qual é o índice que foi analisado                 
                    linha.append(item[6])                   # em qual base foi medida
                    linha.append(letra_alice)               # qual foi a letra enviada por Alice    
                    linha.append(item[9])                   # autovalor
                else:
                    # print("Incorreto: %d %d, indice 1" %(s2_bin, letra_alice))
                    pass

            if len(linha) != 0:
                correct_results.append(linha)

print("As respostas em que Bob conseguiu acertar a letra de Alice são:")
for i in range(len(correct_results)):
    print(correct_results[i])           

# Dicionário para armazenar probabilidades de acerto
probabilidades = {}
media_acerto = []
for i in range(d):
    for j in range(d):
        estado = states[i][j]  # Ex: [0, 0]
        item = all_circuits[i][j]
        
        # Resultados para sigma-x (índice 0)
        counts_x = eval(item[3])  # Converte a string de counts para dicionário
        acertos_x = counts_x.get('0', 0) if estado[0] == 0 else counts_x.get('1', 0)
        prob_x = acertos_x / 1024
        
        # Resultados para sigma-y (índice 1)
        counts_y = eval(item[7])
        acertos_y = counts_y.get('0', 0) if estado[1] == 0 else counts_y.get('1', 0)
        prob_y = acertos_y / 1024
        
        # Armazena no dicionário
        probabilidades[f"Estado {estado} - Medição x (bit 0)"] = prob_x
        probabilidades[f"Estado {estado} - Medição y (bit 1)"] = prob_y

        if(prob_x > 0.8):
            media_acerto.append(prob_x)
        elif(prob_y > 0.8):
            media_acerto.append(prob_y)

# Imprime as probabilidades
# print("\nProbabilidades de acerto para cada caso:")
# for caso, prob in probabilidades.items():
#     print(f"{caso}: {prob:.4f} ({prob * 100:.2f}%)")

print("Probabilidade média de acerto: %.2f" %((sum(media_acerto)/len(media_acerto))*100))