from qiskit import QuantumCircuit, transpile
from qiskit_aer import AerSimulator
from qiskit.quantum_info import Operator
import numpy as np
import math as math
import cmath as cmath
from random import uniform
import picos as pc
import qutip as qt

simulator = AerSimulator()
d = 2           # words size 

def is_unitary(matrix):
    identity_check = np.dot(np.conjugate(matrix.T), matrix)
    return np.allclose(identity_check, np.eye(2))

def initialize_circuit(ampl_list):
    qc = QuantumCircuit(1,1)
    qc.initialize(ampl_list,0)      # initialize the circuit with the qubit Alice sent 
    qc.save_statevector()

    # Transpile for simulator
    qc_aer = transpile(qc, backend=simulator)

    result = simulator.run(qc_aer).result()
    #psi = result.get_statevector()
    return qc

def Pi_random(): # Return a random projector for dimension d
    # Definir theta e phi para o exemplo
    d=2
    theta = np.random.rand(1)[0] * np.pi # Theta da esfera de Bloch
    phi = np.random.rand(1)[0] * 2 * np.pi # Phi da esfera de Bloch
    # criação da unitária U
    out = np.cos(theta)*np.sin(theta)*np.exp(1j*phi)
   
    
    return np.array([[np.cos(theta)**2,out],[out.conjugate(),np.sin(theta)**2]])

# create psi states 
states = []
for i in range(d):
    linha = []
    for j in range(d):
        linha.append([i,j])
    states.append(linha)

# Matriz Densidade Aleatória
rho=np.zeros([d,2]).tolist() 
for x0 in range(d): 
    for x1 in range(d):
        # Passo 1: Criar uma matriz aleatória complexa
        A = np.random.rand(d,d) + 1j * np.random.rand(d,d)  
        # Passo 2: Transformar em uma matriz hermitiana -- matriz A é = a sua transposta conjugada 
        H = (A + A.conj().T) / 2
        # Passo 3: Calcular os autovalores e autovetores
        autovalores, autovetores = np.linalg.eig(H)
        # Ajustar os autovalores para que sejam positivos
        autovalores = np.abs(autovalores)
        # Normalizar os autovalores para que o traço seja 1
        autovalores = autovalores / np.sum(autovalores)
        # Reconstituir a matriz hermitiana com os novos autovalores
        H = autovetores @ np.diag(autovalores) @ autovetores.conj().T
        #salva a matriz H em rho[x0][x1]
        rho[x0][x1]=np.round(H,4)

phi_a = 0                   # na esfera de Bloch, os estados estão cortando no eixo longitudinal 
phi_b = []                  # agrupamento de ângulos aleatórios para ter a otimização depois

# criando ângulos aleatórios para phi_b
for i in range(d):
    linha = []
    for j in range(d):
        linha.append(uniform(0, 2 * 3.141592653589793))
    phi_b.append(linha)

# amplitudes individuais de cada estado
ampl = []
for i in range(d):
    linha = []
    for j in range(d):
        # pegar elementos da diagonal principal de rho --- |a|^2 e |b|^2 
        alpha = np.sqrt(rho[i][j][0][0])
        beta = np.sqrt(rho[i][j][1][1])

        # adicionar a parte imaginária a cada amplitude
        img_1 = alpha*cmath.exp(1j*phi_a)
        img_2 = beta* cmath.exp(1j*phi_b[i][j])
        estado = [img_1, img_2]         # cria uma lista interna de amplitudes para cada estado (+ fácil a separação)
        linha.append(estado)
    ampl.append(linha)

'''
ampl é uma matriz tridimensional, onde cada elemento é uma lista [img_1, img_2] que representa as amplitudes de dois 
estados quânticos associados aos índices i e j
ampl[i][j]: Retorna o estado [img_1, img_2] correspondente aos índices i e j.
ampl[i][j][k]: Retorna o valor específico do estado, onde:
    k = 0 para img_1.
    k = 1 para img_2
'''

# Definir theta e phi para o exemplo
theta = np.random.rand(d*2) * np.pi # Theta da esfera de Bloch
phi = np.random.rand(d*2) * 2 * np.pi # Phi da esfera de Bloch

# M = matriz de decodificação ------------ NÃO ESTÁ SENDO USADO 
# M=[]
# for x in range(d):
#     m1=np.array([np.cos(theta[x]/2),np.exp(1j*phi[x])*np.sin(theta[x]/2)])  # alpha
#     m2=np.array([-np.sin(theta[x]/2),np.exp(1j*phi[x])*np.cos(theta[x]/2)])     # beta 
#     # print("m1:", m1)
#     # print("m2:", m2)
#     # m1 e m2 formam uma base de estados ortogonais quaisquer.
#     # Aqui estou definindo como se fosse os parâmetros iniciais, pois é necessário iniciar com alguma coisa a otimização do sistema.
#     M.append([np.outer(m1,np.conjugate(m1).T),np.outer(m2,np.conjugate(m2).T)])
#     # M[x][y] contém a base que vai decodificar x=x_i onde y=0 (=1) temos x_i=0 (1)

'''
se j = 0 -- M[0][x_0]
se j = 1 -- M[1][X_1]
00 -- rho[0][0] -- se j = 0: M[0][0]; se j = 1: M[1][0]
01 -- rho[0][1] -- se j = 0: M[0][0]; se j = 1: M[1][1]
10 -- rho[1][0] -- se j = 0: M[0][1]; se j = 1: M[1][0]
11 -- rho[1][1] -- se j = 0: M[0][1]; se j = 1: M[1][1]
'''

# criação da unitária
U = []
for x in range(d):
    n = np.array([np.cos(theta[x]/2), -np.sin(theta[x]/2)])
    j = np.array([np.sin(theta[x]/2)*np.exp(1j*phi[x]), np.cos(theta[x]/2)*np.exp(1j*phi[x])])
    linha = []
    linha.append(n)
    linha.append(j)
    U.append(linha)

U = np.array(U)     # transforma U em matriz numpy 

'''
U -- matrix 2x2 (cada uma): U[0] -- primeira letra
                            U[1] -- segunda letra
'''

all_circuits = []           # matriz que vai armazenar todas as informações de todos os processamentos de cada estado da função psi
for i in range(d):
    linha = []
    
    for j in range(len(ampl[i])):
        estado = ampl[i][j]
        img_1, img_2 = estado
        circ_ampl = [img_1, img_2]
        circ_simulation = initialize_circuit(circ_ampl)         # inicializa o circuito para cada amplitude dos estados

        indiv = []              # cria uma linha individual para cada estado 

        indiv.append(states[i][j])
        indiv.append(circ_ampl)        
        
        for indice in range(d):     # gera a medição tanto para o índice 0 quanto para o índice 1 para todas as amplitudes
            if indice == 0 :         # primeira letra
                indiv.append("O qubit foi medido na base sigma x")
            else:                    # segunda letra                     
                indiv.append("O qubit foi medido na base sigma y")
            
            # Selecionar a matriz de decodificação
            decod_matrix = U[indice]

            # Verificar unitariedade
            if not is_unitary(decod_matrix):
                print(f"Erro: A matriz U[{indice}] não é unitária.")
                print("Matriz U[indice]:\n", decod_matrix)
                exit()

            # Continuar apenas se for unitária
            decod_operator = Operator(decod_matrix)
            circ_simulation.unitary(decod_operator, [0], label="Decodificação")

            # print(circ_simulation.data)
            # print(circ_simulation)
            
            # Realização da medição e da simulação de cada circuito 
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

            # Liberação do circuito para o circuito do j = 1
            circ_simulation.data.pop(-1)
            circ_simulation.data.pop(-1)

        linha.append(indiv)
    all_circuits.append(linha)

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
            letra_alice = states[i][j][indice]      # inicializa a letra da Alice, índice por índice 
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

        if(prob_x > 0.5):
            media_acerto.append(prob_x)
        elif(prob_y > 0.5):
            media_acerto.append(prob_y)

# Imprime as probabilidades
print("\nProbabilidades de acerto para cada caso:")
for caso, prob in probabilidades.items():
    print(f"{caso}: {prob:.4f} ({prob * 100:.2f}%)")

print("Probabilidade média de acerto: %.2f" %((sum(media_acerto)/len(media_acerto))*100))